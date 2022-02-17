//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : QuantoCDSAlgorithm.cpp
//
//   Description : How to apply quanto adjustment to CDS Par Spread Curves
//
//   Author      : Mark A Robson
//
//   Date        : November 30, 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif
#include "edginc/QuantoCDSAlgorithm.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CUPSAnalytics.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/SRMFXUtil.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/FindRootBrent.hpp"
#include "edginc/IRGridPointCache.hpp"
#include "edginc/Hashtable.hpp"

DRLIB_BEGIN_NAMESPACE
static const double SPREAD_LOWER_BOUND = -1.e4;
static const double SPREAD_UPPER_BOUND =  1.e10;
class QuantoCDSAlgorithm::Imp{
public:

    //// returns list of forward rates and forward discount factors
    static void ratesAndDiscountFactors(
        IYieldCurveConstSP   curve,           // in
        const DateTime&      today,           // in
        const DateTimeArray& timeline,        // in
        DoubleArray&         timeOffset,      // out
        DoubleArray&         fwdDiscFactor,   // out
        DoubleArray&         rate){           // out
        int numDates = timeline.size();
        timeOffset.resize(numDates+1); // awful but safer to keep it as it was
        timeOffset[0] = 0.0;
        fwdDiscFactor.resize(numDates);
        rate.resize(numDates);
        // fast indexing off yield curve
        auto_ptr<IYieldCurve::IKey> key(curve->logOfDiscFactorKey());
        for (int i = 0; i < numDates; i++){
            const DateTime& start = i == 0? today: timeline[i-1];
            fwdDiscFactor[i] = exp(key->calc(start, timeline[i]));
            // Use continuosly compounded rates
            timeOffset[i+1] = Actual365F::yearFraction(start, timeline[i]);
            // will need to revisit this as this could happen
            if (Maths::isZero(timeOffset[i+1])){
                throw ModelException("QuantoCDSAlgorithm::Imp::"
                                     "ratesAndDiscountFactors",
                                     "Zero year fraction on timeline between "+
                                     start.toString()+" and "+
                                     timeline[i].toString());
            }
            rate[i] = RateConversion::discountToRateYearFrac(
                fwdDiscFactor[i], timeOffset[i+1], CompoundBasis::CONTINUOUS);
            if (i > 0){
                // timeOffset in the original code is yearFrac from today.
                // Although it seems that you then have to go round subtracting
                // the previous value off to get the interval...
                timeOffset[i+1] += timeOffset[i];
            }
        }
    }
    
    /** this function returns product of spread spot vols and spread
        forward rates
        lBpVol - spread "basis point" vol */
    static void spreadBasisPointVol(
        const DateTime&      today,    // the first point on the timeline
        const DateTimeArray& allDates, /* all points on the timeline
                                          but the first one */
        DefaultRatesConstSP  spreadCurve,    // Domestic spread curve
        const DoubleArray&   spreadVolCurve, // Domestic spread vol curve
        DoubleArray&         lBpVol){    //out
        int numDates = allDates.size();
        lBpVol.resize(numDates);
        Actual365F act365F;
        for (int i = 0; i < numDates; i++) {
            const DateTime& startDate = i == 0? today: allDates[i-1];
            const DateTime& endDate = allDates[i];
            // probably want more methods on DefaultRates
            double defaultPV = spreadCurve->calcDefaultPV(startDate, endDate);
            double defaultFwdRate = RateConversion::discountToRate(
                defaultPV, startDate, endDate, &act365F, 
                CompoundBasis::ANNUAL); // why ANNUAL?
            lBpVol[i] = spreadVolCurve[i] * defaultFwdRate;
        }
    }

    /** Adjust clean spreads for adjustment due to IR correct. Algorithm is 
        clean spread -> CDS pricing -> continuous fee par spreads ->
        -> adjusted CDS bootstrapping -> clean spread -> 
        -> subtract rates covariance adjustment -> result.
        From SpreadCurveAlgorithms.cpp: SpreadCurveIRCorrelationAdjustment.
        NB Here domestic means original currency of cds par spreads */
    static DefaultRatesConstSP irCorrelationAdjustment(
        const DateTime&        today,       // not included in timeline
        const DateTimeArray&   timeline,    // when quantities below are defined
        DefaultRatesConstSP    spreadCurve,    // Domestic spread curve
        const DoubleArray&     spreadVolCurve, // Domestic spread vol curve
        double                 spreadMeanReversion, // for spread vol
        SRMRatesUtil&          ir,            // Domestic IR stuff
        double                 spreadIRCorr,    // Correlation of spread and IR
        bool                   floorIntensity){  // floored intensity
        const string method("QuantoCDSAlgorithm::Imp::irCorrelationAdjustment");
        if (Maths::isZero(spreadIRCorr)){
            return spreadCurve;
        }
        // The integration is correct for flat forwards and constant spot
        // vol.  We do not enforce these conventions, however.
        
        // Set up a datelist of all "critical dates" corresponding to dates in
        // vol curves and zero curves.  By making the timeline as fine as all
        // the dates in all curves, we are ensuring that we have constant rates
        // and vols (for flat forwards and constant spot vol integration).
        
        double betaL = spreadMeanReversion;
        double betaR = ir.getBeta(0 /* only one factor here */);

        DateTimeArraySP spreadCurveDates(spreadCurve->getDefaultDates());
        // fudge the time of day
        DateTime::setTimeOfDay(*spreadCurveDates, today.getTime());
        
        DoubleArray    fwdRlDF; // forward discount factor
        DoubleArray    rRate;  // interest rate ( r )
        DoubleArray    tOffset;
        vector<int> maturityIndices(DateTime::getIndexes(timeline,
                                                         *spreadCurveDates));
        ratesAndDiscountFactors(ir.getDiscYC(), today, timeline,
                                tOffset, fwdRlDF, rRate);
        // default rate ( lambda )
        // to do: put appropriate methods on DefaultRates
        CashFlowArraySP fwdRates = spreadCurve->getCleanSpreadCurve();
        // basis point vols
        const DoubleArray& rBpVol = ir.basisPointVol(timeline);
        DoubleArray lBpVol;
        spreadBasisPointVol(today, timeline, 
                            spreadCurve, spreadVolCurve,
                            lBpVol);

        UnadjustedCDSFunction unadjCDS(tOffset, fwdRlDF, rRate);
        auto_ptr<CUPSAnalyticsBase> adjCDS(
            CUPSAnalyticsBase::create(betaR, betaL));
        adjCDS->initialize( 
            tOffset, fwdRlDF, rRate,
            rBpVol, lBpVol, spreadIRCorr);

        DoubleArray result(spreadCurveDates->size());
        double covAdjusment = 0; // = CRK(tn,tn)
        double startOffset = 0; // first interval starts at t0
        Actual365F dcc;
        for (int benchmarkIndex = 0; 
             benchmarkIndex < spreadCurveDates->size(); benchmarkIndex++) {
            
            int maturityIndex = maturityIndices[benchmarkIndex];
            // calculate the probability of not defaulting
            double unadjustedLambda = (*fwdRates)[benchmarkIndex].amount;
            unadjCDS.nextBenchmark(maturityIndex);
            double continuousParSpread = unadjCDS.parSpread(unadjustedLambda);
            
            adjCDS->nextBenchmark(maturityIndex, continuousParSpread);
            double adjustedLambda;
            try{
                adjustedLambda = FindRootBrent::solve(*adjCDS,
                                                      unadjustedLambda,
                                                      SPREAD_LOWER_BOUND,
                                                      SPREAD_UPPER_BOUND).x;
            } catch (exception& e){
                throw ModelException(e, method, "Failed to solve on "+
                                     (*spreadCurveDates)[benchmarkIndex].
                                     toString());
            }

            double newCovAdjusment = adjCDS->SpreadIrCov(maturityIndex);
            
            adjustedLambda += 
                (newCovAdjusment - covAdjusment) / 
                (tOffset[maturityIndex+1] - startOffset);
            
            if (!floorIntensity && adjustedLambda < 0){
                string m("Adjusted default probability that corresponds "
                         "to benchmark with maturity " + 
                         (*spreadCurveDates)[benchmarkIndex].toString() + 
                         " is negative.");
                throw ModelException(method, m);
            }
            
            result[benchmarkIndex] = Maths::max(adjustedLambda, 0.0);
            
            covAdjusment = newCovAdjusment;
            startOffset = tOffset[maturityIndex+1];
        }
        return DefaultRatesConstSP(
            new CDSHelper::CParSpreadDefaultRates(today, *spreadCurveDates, 
                                                  result));
    }

    /** algorithm is 
     clean spread adjusted for domestic IR spread -> 
     add rates and FX covariance adjustment -> 
     clean spread clean spread adjusted in foreign ccy-> 
     adjusted CDS pricing -> 
     continuous fee par spreads -> CDS bootstrapping -> result
     NB Here foreign means instrument ccy, domestic means original
     currency of cds par spreads */
    static DefaultRatesConstSP irAdjustedToCCYAdjusted(
        const DateTime&       today,       // not included in timeline
        const DateTimeArray&  timeline,    // when quantities below are defined
        DefaultRatesConstSP   spreadCurve, /* the clean spread curve with
                                              removed effect of correlation
                                              between spreads and domestic
                                              interest rates. */
        DefaultRatesConstSP   domesticSpreadCurve,  // Domestic spread curve
        const DoubleArray&    spreadVolCurve,       // Domestic spread vol curve

        double                spreadMeanReversion, // for spread vol
        SRMRatesUtil&         ir,            // Foreign IR stuff

        const DoubleArray&    fxVolCurve,    // FX vol curve

        double                spreadIRCorr,  /* Correlation of spread and
                                                foreign IR */
        double                spreadFxCorr,  // Correlation of spread and FX
        bool                  floorIntensity) {
        static const string method("QuantoCDSAlgorithm::Imp::"
                                   "irAdjustedToCCYAdjusted");
        if (Maths::isZero(spreadIRCorr) && Maths::isZero(spreadFxCorr)){
            return spreadCurve;
        }

        // The integration is correct for flat forwards and constant spot
        // vol.  We do not enforce these conventions, however.

        // Set up a datelist of all "critical dates" corresponding to dates in
        // vol curves and zero curves.  By making the timeline as fine as all
        // the dates in all curves, we are ensuring that we have constant rates
        // and vols (for flat forwards and constant spot vol integration).

        double betaL = spreadMeanReversion;
        double betaR = ir.getBeta(0 /* only one factor here */);
        DateTimeArraySP spreadCurveDates(spreadCurve->getDefaultDates());
        // fudge the time of day
        DateTime::setTimeOfDay(*spreadCurveDates, today.getTime());
        DoubleArray    fwdRlDF; // forward discount factor
        DoubleArray    rRate;  // interest rate ( r )
        DoubleArray    tOffset; // NB from 0, ... to timeline.size()+1
        vector<int> maturityIndices(DateTime::getIndexes(timeline,
                                                         *spreadCurveDates));
        ratesAndDiscountFactors(ir.getDiscYC(), today, timeline,
                                tOffset, fwdRlDF, rRate);
        // default rate ( lambda )
        // to do: put appropriate methods on DefaultRates
        CashFlowArraySP fwdRates = spreadCurve->getCleanSpreadCurve();
        DoubleArray    lRate(fwdRates->size());
        for (int i = 0; i < fwdRates->size(); i++){
            lRate[i] = (*fwdRates)[i].amount;
        }
        // basis point vols
        const DoubleArray& rBpVol = ir.basisPointVol(timeline);
        DoubleArray    lBpVol;
        spreadBasisPointVol(today, timeline,
                            domesticSpreadCurve, spreadVolCurve,
                            lBpVol );
        DoubleArray result(spreadCurveDates->size());
        UnadjustedCDSFunction unadjCDS(tOffset, fwdRlDF, rRate );
        auto_ptr<CUPSAnalyticsBase> adjCDS(CUPSAnalyticsBase::create(betaR,
                                                                     betaL));
        adjCDS->initialize(tOffset, fwdRlDF, rRate,
                           rBpVol, lBpVol, spreadIRCorr );
    
        double covAdjusment = 0; // = CRK(tn,tn) - CXK(tn)
        double startOffset = 0; // first interval starts at t0
        for(int benchmarkIndex = 0; 
            benchmarkIndex < spreadCurveDates->size(); benchmarkIndex++ ) {

            int maturityIndex = maturityIndices[benchmarkIndex];
            double adjustedLambda = lRate[benchmarkIndex];

            double newCovAdjusment = adjCDS->SpreadIrCov(maturityIndex) - 
                adjCDS->SpreadFxCov(fxVolCurve, spreadFxCorr, maturityIndex);

            adjustedLambda -= 
                (newCovAdjusment - covAdjusment) / 
                (tOffset[maturityIndex+1] - startOffset);

            if (!floorIntensity && adjustedLambda < 0){
                string m("Adjusted default probability that corresponds "
                         "to benchmark with maturity " + 
                         (*spreadCurveDates)[benchmarkIndex].toString() + 
                         " is negative.");
                throw ModelException(method, m);
            }

            adjCDS->nextBenchmark(maturityIndex);
            double continuousParSpread = adjCDS->parSpread(adjustedLambda); 

            unadjCDS.nextBenchmark(maturityIndex, continuousParSpread );
            double unadjustedLambda;
            try{
                unadjustedLambda = 
                    FindRootBrent::solve(unadjCDS,
                                         adjustedLambda,
                                         SPREAD_LOWER_BOUND,
                                         SPREAD_UPPER_BOUND,
                                         FindRoot::X_TOLERANCE,
                                         FindRoot::F_TOLERANCE,
                                         FindRoot::MAX_NUM_OF_ITERATIONS,
                                         fabs(adjustedLambda)*1.e-4 + 1.e-8).x;
            } catch (exception& e){
                throw ModelException(e, method, "Failed to solve on "+
                                     (*spreadCurveDates)[benchmarkIndex].
                                     toString());
            }
            result[benchmarkIndex] = Maths::max(unadjustedLambda, 0.0);
            covAdjusment = newCovAdjusment;
            startOffset = tOffset[maturityIndex+1];
        }
        return DefaultRatesConstSP(
            new CDSHelper::CParSpreadDefaultRates(today, *spreadCurveDates, 
                                                  result));
    }
};
/** Just a class for returning debug results */
class QuantoCDSAlgorithm::Debug: public CObject{
public:
    static CClassConstSP const TYPE;

    DateTime        today;
    DateTimeArray   timeline;
    DoubleArray     domIRVolCurve;
    DoubleArray     forIRVolCurve;
    DoubleArray     fxVolCurve;
    DoubleArray     spreadVolCurve;
    CashFlowArraySP adjustedCleanSpreads;
    CashFlowArraySP quantoCleanSpreads;

    Debug(): CObject(TYPE){}

private:
    static IObject* defaultConstructor(){
        return new Debug();
    }

    static void load(CClassSP& clazz){
        REGISTER(Debug, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(today, "today");
        FIELD(timeline, "timeline used to compute quantoed curve");
        FIELD(domIRVolCurve, "domestic IR spot vol");
        FIELD(forIRVolCurve, "foreign IR spot vol");
        FIELD(fxVolCurve, "fx spot vol");
        FIELD(spreadVolCurve, "spread spot vols");
        FIELD(adjustedCleanSpreads, "orig spreads adjusted for IR correlation");
        FIELD(quantoCleanSpreads, "clean spreads in new [domestic] ccy");
    }
};

CClassConstSP const QuantoCDSAlgorithm::Debug::TYPE =
CClass::registerClassLoadMethod("QuantoCDSAlgorithm::Debug", typeid(Debug),
                                load);

QuantoCDSAlgorithm::~QuantoCDSAlgorithm(){}
    
/** Return clean spreads in domestic currency from clean spreads in foreign
    currency. NB Here domestic means instrument ccy, foreign means original
    currency of cds par spreads */
DefaultRatesConstSP QuantoCDSAlgorithm::currencyAdjust(
    IYieldCurveConstSP      domYC,   // set to Projection curve
    IYieldCurveConstSP      forYC,   // set to Projection curve
    ICDSParSpreadsConstSP   parSpreads,
    FXAssetConstSP          fx,
    CorrelationConstSP      corrFxCDS,      // fx and cds par spreads
    CorrelationConstSP      corrFxDom,      // fx and dom IR
    CorrelationConstSP      corrFxFor,      // fx and for IR
    CorrelationConstSP      corrCDSDom,     // cds and domestic iR
    CorrelationConstSP      corrCDSFor,     // cds and foreign IR
    CorrelationConstSP      corrDomFor) const { // dom and for IR
    static const string method("QuantoCDSAlgorithm::currencyAdjust");
    try{
        const DateTime& today = fx->getToday();
        // get hold of the cds clean spreads and the relevant dates
        DefaultRatesConstSP   spreadCurve(parSpreads->defaultRates());
        DateTimeArraySP spreadCurveDates(spreadCurve->getDefaultDates());
        // fudge the time of day (the approach is to force all relevant dates
        // to have the same time of day as today)
        DateTime::setTimeOfDay(*spreadCurveDates, today.getTime());

        // get hold of processed vol, domestic
        VolProcessedBSIRSP domProcessedVol = getIRProcessedVol(domYC);
        // and then foreign
        VolProcessedBSIRSP forProcessedVol = getIRProcessedVol(forYC);

        // 1. Calculate spot vols for both YCs
        /* Do this essentially by building a SRMRatesUtil and calling getSpotVols */
        // here foreign means the currency of the original CDS Par Spreads
        SRMRatesHJMUtilSP forIR(new SRMRatesHJMUtil(today, fx,
                                    modelParamsKey, forProcessedVol, 
                                    IYieldCurveConstSP::dynamicCast(fx->getRiskCcy()),
                                    forYC, skipIRBadVols, corrSwapStart,
                                    corrSwapMat, corrSwapDCC, corrSwapFreq));
        const DateTimeArray& forTimeline = forIR->getExtendedTimeLine();
        // here domestic means the currency of the instrument (which wants to
        // use the CDS Par Spreads curve)
        SRMRatesHJMUtilSP domIR(new SRMRatesHJMUtil(today, fx,
                                    modelParamsKey, domProcessedVol,
                                    IYieldCurveConstSP::dynamicCast(fx->getBaseCcy()),
                                    domYC, skipIRBadVols, corrSwapStart,
                                    corrSwapMat, corrSwapDCC, corrSwapFreq));
        const DateTimeArray& domTimeline = domIR->getExtendedTimeLine();
    
        // 2. Retrieve cds vol dates
        MRSpotVolRequest cdsRequest;
        CVolProcessedSP processed(parSpreads->getProcessedVol(&cdsRequest));
        MRSpotVolProcessedSP cdsVol(
            MRSpotVolProcessedSP::dynamicCast(processed));
        DateTimeArraySP cdsVolDates(cdsVol->getSpotVolDates());
        // fudge the time of day
        DateTime::setTimeOfDay(*cdsVolDates, today.getTime());
        
        // 3. Calculate spot vols for fx 
        // Do this by building a SRMFXUtil object. 
        SRMFXUtilSP fxUtil(new SRMFXUtil(
            forIR, 
            domIR, 
            fx, 
            corrFxFor->getCorrelation(),
            corrFxDom->getCorrelation(),
            corrDomFor->getCorrelation(),
            fxVolBootstrapMode,
            fxCutOffLevel,
            // we really don't care about the type - need to revisit. to do.
            SRMFXVol::VOL_TYPE->getName(),
            true) );
        // when created fxUtil doesn't assume there are dates in the timeline
        // moreover, it is SRMFXDiffuse asset that informs Util about dates it is interested in via setTimeLine()
        // here we use special case of RatesUtil that knows its dates, so we update fxUtil with the same dates.
        fxUtil->setTimeLine(DateTimeArrayConstSP(new DateTimeArray(domIR->getSimDates())));
        DateTimeArray  fxTimeline(fxUtil->calibratedSpotVolDates());
        // fudge the time of day
        DateTime::setTimeOfDay(fxTimeline, today.getTime());
        const vector<double>& fxSpotVols = fxUtil->calibratedSpotVols();

        // now build our combined timeline
        vector<const DateTimeArray*> dtArrayVector(4);
        dtArrayVector[0] = &forTimeline;
        DateTimeArray forZeroDates(forYC->zeroDates());
        // fudge the time of day
        DateTime::setTimeOfDay(forZeroDates, today.getTime());
        
        dtArrayVector[1] = &forZeroDates;
        dtArrayVector[2] = spreadCurveDates.get();
        dtArrayVector[3] = cdsVolDates.get();
        DateTimeArray timeline1(DateTime::merge(dtArrayVector));
        // remove today from timeline
        if (timeline1.front() != today){
            throw ModelException(method, "Internal error. First date in SRM "
                                 "timeline is not today");
        }
        timeline1.erase(timeline1.begin());
        // Retrieve cds vols - ask for them on timeline
        DoubleArray spreadVols1(timeline1.size());
        cdsVol->spotVol(today, timeline1, spreadVols1);
        double spreadMeanReversion = cdsVol->meanReversion();
        // 4. call method to adjust clean spreads for 
        // adjustment to IR correlation
        DefaultRatesConstSP adjustedCleanSpreads(
            Imp::irCorrelationAdjustment(today,
                                         timeline1,
                                         spreadCurve,
                                         spreadVols1,
                                         spreadMeanReversion,
                                         *forIR,
                                         corrCDSFor->getCorrelation(),
                                         floorIntensity));
        // build time line for next call
        dtArrayVector.resize(5);
        dtArrayVector[0] = &domTimeline;
        DateTimeArray domZeroDates(domYC->zeroDates());
        // fudge the time of day
        DateTime::setTimeOfDay(domZeroDates, today.getTime());
        dtArrayVector[1] = &domZeroDates;
        dtArrayVector[2] = &fxTimeline;
        dtArrayVector[3] = spreadCurveDates.get();
        dtArrayVector[4] = cdsVolDates.get();
        DateTimeArray timeline2(DateTime::merge(dtArrayVector));
        // and then extend our spot vols along this timeline
        vector<double> fxVolCurve(SRMUtil::extendVol(fxTimeline, fxSpotVols, timeline2));
        // remove today from timeline
        if (timeline2.front() != today){
            throw ModelException(method, "Internal error. First date in SRM "
                                 "timeline is not today");
        }
        timeline2.erase(timeline2.begin());
        if (timeline2.size() != (int) fxVolCurve.size()){
            throw ModelException(method, "Internal error. Arrays wrong length");
        }
        DoubleArray spreadVols2(timeline2.size());
        cdsVol->spotVol(today, timeline2, spreadVols2);

        // 5. call method to adjust clean spreads due to switch in ccy
        DefaultRatesConstSP quantoSpreads(
            // note switch in terminology regarding domestic and foreign. Here
            // the domestic currency is what we quanto into. Whilst in CMLib
            // it's the other way round here. Therefore the sign of the 
            // correlation needs to be reversed since the CMLib code is for
            // where the fx is quoted as number of domestic per unit of foreign
            // using CMLib terminology for domestic and foreign.
            Imp::irAdjustedToCCYAdjusted(today,
                                         timeline2,
                                         adjustedCleanSpreads,
                                         spreadCurve,
                                         spreadVols2,
                                         spreadMeanReversion,
                                         *domIR,
                                         DoubleArray(fxVolCurve.begin(),
                                                     fxVolCurve.end()),
                                         corrCDSDom->getCorrelation(),
                                         -corrFxCDS->getCorrelation(),
                                         floorIntensity));
        if (debug){
            debugOutput.reset(new Debug());
            debugOutput->today = today;
            debugOutput->timeline = fxTimeline;
            debugOutput->domIRVolCurve = domIR->getCalibratedSpotVols();
            debugOutput->forIRVolCurve = forIR->getCalibratedSpotVols();
            debugOutput->fxVolCurve = DoubleArray(fxSpotVols.begin(),
                                                  fxSpotVols.end());
            debugOutput->spreadVolCurve = spreadVols2;
            // the whole interface needs reviewing ...
            const CDSHelper::CParSpreadDefaultRates* defRates = 
                &dynamic_cast<const CDSHelper::CParSpreadDefaultRates&>(
                    *adjustedCleanSpreads);
            debugOutput->adjustedCleanSpreads = 
                defRates->annualiseDefaultRates()->getCleanSpreadCurve();
            defRates = &dynamic_cast<const CDSHelper::CParSpreadDefaultRates&>(
                *quantoSpreads);
            debugOutput->quantoCleanSpreads = 
                defRates->annualiseDefaultRates()->getCleanSpreadCurve();
        }
        return quantoSpreads;
    } catch (exception& e){
        throw ModelException(e, method,
                             "When applying quanto conversion for "+
                             parSpreads->getName()+" from "+forYC->getCcy()+
                             " to "+domYC->getCcy());
    }
}

/** Returns debug info - may be null eg if the algorithm has not been 
    used. Object must have been initialised with debug flag on */
IObjectSP QuantoCDSAlgorithm::getDebugInfo() const{
    return debugOutput;
}

/** Constructor */
QuantoCDSAlgorithm::QuantoCDSAlgorithm(
    const string&      modelParamsKey, // for domestic and foreign IR
    const string&      calibrationStyle,// for domestic and foreign IR
    const string&      calibrationMaturity,// for domestic and foreign IR
    bool               skipIRBadVols,// for domestic and foreign IR
    const string&      fxVolBootstrapMode,
    double             fxCutOffLevel,
    bool               floorIntensity, // true: floor negative vols at 0
    string             corrSwapStart, // eg 1Y  (offset to today)
    string             corrSwapMat,   // eg 10Y (offset to start)
    string             corrSwapDCC,   // eg Act/365F
    string             corrSwapFreq, // eg 6M
    bool               debug): // true: store extra info
    modelParamsKey(modelParamsKey),
    calibrationStyle(calibrationStyle),
    calibrationMaturity(calibrationMaturity), skipIRBadVols(skipIRBadVols),
    fxVolBootstrapMode(fxVolBootstrapMode), fxCutOffLevel(fxCutOffLevel),
    floorIntensity(floorIntensity), corrSwapStart(corrSwapStart),
    corrSwapMat(corrSwapMat), corrSwapDCC(corrSwapDCC), 
    corrSwapFreq(corrSwapFreq), debug(debug){}

/** The supplied IRGridPointCache will be updated as IRVols are used.
    This replaces any previous IRGridPointCache. Note no clone of the
    supplied parameter is taken */
void QuantoCDSAlgorithm::setIRGridPointCache(IRGridPointCacheSP irGridPtsCache){
    this->irGridPtsCache = irGridPtsCache;
}


//------------------------------
// Cache-related methods
//------------------------------
/** Hash code function - the cache needs improved performance compared
 * to the default "hashCode" function in CObjet: only the required
 * components are hashed */
int QuantoCDSAlgorithm::hashCodeOpt() const {
    int hcode = hash_string(modelParamsKey); // for domestic and foreign IR
    hcode ^= hash_string(calibrationStyle);// for domestic and foreign IR
    hcode ^= hash_string(calibrationMaturity);// for domestic and foreign IR
    hcode ^= hash_string(fxVolBootstrapMode);
    hcode ^= hash_string(corrSwapStart); // eg 1Y  (offset to today)
    hcode ^= hash_string(corrSwapMat);   // eg 10Y (offset to start)
    hcode ^= hash_string(corrSwapDCC);   // eg Act/365F
    hcode ^= hash_string(corrSwapFreq); // eg 6M
    hcode ^= CDouble::hashCode(fxCutOffLevel);
    hcode ^= CBool::hashCode(skipIRBadVols);// for domestic and foreign IR
    hcode ^= CBool::hashCode(floorIntensity); // true: floor negative vols at 0
    return hcode;
}


/** Comparison function - the cache needs improved performance compared 
 * to the default "equalTo" function in CObjet: only the required
 * components are compared */
bool QuantoCDSAlgorithm::equalToOpt(
    const QuantoCDSParSpreads::IAlgorithm* algorithm2) const 
{
    if (this == algorithm2) { // Obvious first test
        return true;
    }
    if (!algorithm2) {
        return false;
    }

    const QuantoCDSAlgorithm* algo2 = 
        dynamic_cast<const QuantoCDSAlgorithm*>(algorithm2);

    if (!algo2) {
        return false;
    }

    if ((modelParamsKey != algo2->modelParamsKey)           ||
        (calibrationStyle != algo2->calibrationStyle)       ||
        (calibrationMaturity != algo2->calibrationMaturity) ||
        (fxVolBootstrapMode != algo2->fxVolBootstrapMode)   ||
        (corrSwapStart != algo2->corrSwapStart)             ||
        (corrSwapMat != algo2->corrSwapMat)                 ||
        (corrSwapDCC != algo2->corrSwapDCC)                 ||
        (corrSwapFreq != algo2->corrSwapFreq)               ||
        (fxCutOffLevel != algo2->fxCutOffLevel)             ||
        (skipIRBadVols != algo2->skipIRBadVols)             ||
        (floorIntensity != algo2->floorIntensity))
    {
        return false;
    }
    return true;
}

//-------------------------------------
// IR Grid point cache related methods
//-------------------------------------

/** Factor out the caching of IR vol points */
void QuantoCDSAlgorithm::cacheGridPoints(const IYieldCurveConstSP domYC,
                                         const IYieldCurveConstSP forYC) const
{
    if (irGridPtsCache.get())
    {
        //be paranoid
        if (!!domYC)
        {
            VolProcessedBSIRSP domProcessedVol = getIRProcessedVol(domYC);
            irGridPtsCache->cacheGridPoints(domYC, domProcessedVol);
        }

        if (!!forYC)
        {
            VolProcessedBSIRSP forProcessedVol = getIRProcessedVol(forYC);
            irGridPtsCache->cacheGridPoints(forYC, forProcessedVol);
        }
    }
}

/** Factor out the retrieval of IR vols */
VolProcessedBSIRSP QuantoCDSAlgorithm::getIRProcessedVol(const IYieldCurveConstSP irYC) const
{
    ExpirySP calibExpiry(new MaturityPeriod(calibrationMaturity));
    CVolRequestSP request(new SwapMaturityVolRequest(calibExpiry.get(), 
                                                     calibrationStyle));

    // get hold of processed vol
    CVolProcessedSP volProcessed(irYC->getProcessedVol(request.get()));
    VolProcessedBSIRSP processedVol(VolProcessedBSIRSP::dynamicCast(volProcessed));

    return processedVol;
}

DRLIB_END_NAMESPACE
