//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditMetricsModel.cpp
//
//   Description : Simple closed form model for Credit Metrics
//
//   Author      : Antoine Gregoire
//
//   Date        : April 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CreditMetricsModel.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/ConvolutionEngine.hpp"
#include "edginc/CCMLossUnit.hpp"
#include "edginc/CreditMetricsLossCalculator.hpp"
#include "edginc/CCMLossCalculator.hpp"
#include "edginc/CreditMetricsFastLossCalculator.hpp"
#include "edginc/CCMFastLossCalculator.hpp"
#include "edginc/FixedTrancheLossCalculator.hpp"
#include "edginc/MarketDataFetcherCIS.hpp"
#include "edginc/CreditIndexBasis.hpp"
#include "edginc/Results.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/NToDefaultLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

/** CreditMetricsModel has been rewritten to automatically convert itself into
 *  a ConvolutionEngine (so for all intents and purposes it will look like a 
 *  ConvolutionEngine, even in instanceOf calls). Previously this has been 
 *  done by implementing RunMulti and delegating that call to a ConvultionEngine. 
 *  However this did not only delegate the pricing, but also the flow of control, 
 *  which is a bad thing.
 *  Because of the the automatic conversion CreditMetricsModel is no longer
 *  implementing IModel, but CObject.
 *  The conversion is done in the "convert" method.
 *  Linus Thand, July 2006
*/


//// Automatically convert to an instance of ConvolutionEngine.
void CreditMetricsModel::convert(IObjectSP&    object,
                                 CClassConstSP requiredType) const {
    static const string method = "CreditMetricsModel::convert";
    try {
        // the condition ConvolutionEngine::TYPE != requiredType has been added 
        // for internal use of that method in CDOQuotesBCGenerator
        // TO REVIEW
        if (ConvolutionEngine::TYPE != requiredType && ConvolutionEngine::TYPE->isAssignableFrom(requiredType)) {
            throw ModelException(method,
                                 "Cannot convert a CreditMetricsModel into "
                                 "object of type "+requiredType->getName());
        }
        object = IObjectSP(
            new ConvolutionEngine(IConvolutionModelSP::dynamicCast(object),
                                  DateTime(), //Null
                                  creditChargeViewType,
                                  pvToSpot, cfCutOffDate,
                                  calibrationStyle,
                                  calibrationMaturity));        
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CreditMetricsModel::~CreditMetricsModel(){}

CreditMetricsModel::CreditMetricsModel(const CClassConstSP& clazz): 
    CObject(clazz),
    pvToSpot(false),
    maxNbSlice(5000), gcdDivisor(20),
    lossUnitOverride(-1), fastConvThreshold(90), quickGreeks("DEFAULT"),
    // sensible default for ir vol interpolation (only used for
    // quanto adjustment)
    calibrationStyle("CMS"), calibrationMaturity("10Y"),
    calculateIndexBasis(0),creditChargeViewType(ConvolutionProduct::CC_AGGRESSIVE),
    effSpreadFreq("1Y"),lossInterpolation(EffectiveCurve::FLAT_FORWARD)
{
    //this functionality is bypassed for the timebeing
    //convCache = ConvolutionCacheSP(new ConvolutionCache());
}

/** Invoked by the containing model before the instrument data is fetched
    ie before CInstrument::GetMarket is invoked. Uses for indexMap data */
void CreditMetricsModel::preFetchMarketData(IModel*           model,
                                            MarketDataConstSP market)
{
    CModelSP newModel = CModelSP( dynamic_cast<CModel*>(model) );
    if (!indexMap.isEmpty()) {
        indexMap.getData(model, market.get());
        //and then clear out the fetcher since the credit index
        //it contains will not be setup, and of no use
        newModel->releaseMDF();
    }   
}


/** Invoked by the containing model after the instrument data is fetched
    ie after CInstrument::GetMarket is invoked. This implementation 
    sorts out index data */
void CreditMetricsModel::postFetchMarketData(IModel*            model,
                                             MarketDataConstSP market)
{
    if (!indexMap.isEmpty()) {
        //only get if it hasnt already been
        if (!indexMap.getSP()) {
            indexMap.getData(model, market.get());
        }
    }
    
}

/** Create a MarketDataFetcher which will be used for retrieving
    market data etc */
MarketDataFetcherSP CreditMetricsModel::createMDF() const {
    if (indexMap.isEmpty()) {
        return MarketDataFetcherSP(new MarketDataFetcherCDS(true));
    } 
    else {
        return MarketDataFetcherSP(
            new MarketDataFetcherCIS(true,
                                    indexMap.getSP(),
                                    indexMapTreatment,
                                    calculateIndexBasis->boolValue()));
    }
}

IModel::WantsRiskMapping CreditMetricsModel::wantsRiskMapping() const {
    return IModel::riskMappingIrrelevant;
}

/** Returns true if the 'fast' convolution should be used */
bool CreditMetricsModel::useFastConvolution(
    CreditTrancheLossConfigConstSP tranche) const
{
    return (fastConvThreshold < tranche->numInnerLossConfigs());
}


void CreditMetricsModel::validatePop2Object() {
    static const string method = "CreditMetricsModel::validatePop2Object";

    // called immediately after object constructed
    try {
        effSpreadFreqObj.reset(new MaturityPeriod(effSpreadFreq));

        if (creditChargeViewType != ConvolutionProduct::CC_CONSERVATIVE &&
            creditChargeViewType != ConvolutionProduct::CC_AGGRESSIVE)
        {
            throw ModelException(method, "Counterparty: creditCharge View "
                                 "type flag must be "+
                                 ConvolutionProduct::CC_CONSERVATIVE+" or "+
                                 ConvolutionProduct::CC_AGGRESSIVE);
        }
        
        /* check model param */
        if (!Maths::equals(lossUnitOverride, -1.0) && 
            !Maths::isPositive(lossUnitOverride))
        {
            throw ModelException(method, "Loss unit override "+
                                 Format::toString(lossUnitOverride)+
                                 " invalid");
        }
        if (gcdDivisor <=0) {
            throw ModelException(method, "Loss unit gcdDivisor "+
                                 Format::toString(gcdDivisor)+"invalid");
        }
        /* check model param */
        if (maxNbSlice < 0) {
            throw ModelException(method, "Loss unit max slice "+
                                 Format::toString(maxNbSlice)+" invalid");
        }
        if (lossInterpolation != EffectiveCurve::FLAT_FORWARD && 
            lossInterpolation != EffectiveCurve::LINEAR)
        {
            throw ModelException("Interpolation method flag invalid");
        }

        if (fastConvThreshold < 0) {
            throw ModelException(method, 
                                 "Fast convolution threshold must be positive");
        }

        //check quick greeks input
        if ((quickGreeks != QG_DEFAULT) &&
            (quickGreeks != QG_YES) &&
            (quickGreeks != QG_NO))
        {
            throw ModelException(method,
                "quickGreeks must be " +
                QG_DEFAULT + ", " +
                QG_YES + " or " +
                QG_NO + ".");
        }

        //check indexMapTreatment input
        if (!indexMap.isEmpty()) {
            if ((indexMapTreatment != CreditIndexBasis::APPLY_TO_NONE) &&
                (indexMapTreatment != CreditIndexBasis::APPLY_TO_MAPPED_NAMES_ONLY) &&
                (indexMapTreatment != CreditIndexBasis::APPLY_TO_ALL_NAMES))
            {
                throw ModelException(method,
                    "indexMapTreatment must be " +
                    CreditIndexBasis::APPLY_TO_NONE + ", " +
                    CreditIndexBasis::APPLY_TO_MAPPED_NAMES_ONLY + " or " +
                    CreditIndexBasis::APPLY_TO_ALL_NAMES + ".");
            }

            //check if calculate flag has been passed
            if (!calculateIndexBasis.get()) {
                //default to true since the index map has been provided
                calculateIndexBasis.reset(CBool::create(true));
            }
        }
        else {
            //check that calculate flag has been passed
            if (!!calculateIndexBasis.get()) {
                throw ModelException(method,
                    "calculateIndexBasis has been supplied without an indexMap.");
            }
        }
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::1validatePop2Object");
    }
}


/** Computes the loss unit (bin size for loss distribution discretization)
    as dictated by this model (NB model can override it) */
double CreditMetricsModel::calculateLossUnit(
    CreditTrancheLossConfigConstSP tranche) const
{
    if (lossUnitOverride > 0.0) {
        return lossUnitOverride;
    }
    // this is a bit weak too - copying all the name notionals to an array
    int numNames = tranche->numInnerLossConfigs();
    DoubleArray notionals(numNames);
    for (int i = 0; i < numNames; ++i) {
        notionals[i] = tranche->nameNotional(i);
    }
    return CCMLossUnit::calcLossUnit(notionals, maxNbSlice, gcdDivisor);
}



/** Generate the timeline on which the effective curve will be calculated */
DateTimeArraySP CreditMetricsModel::generateTimeline(
    const DateTime& today,
    const DateTime& lastObservationDate) const
{
    return CreditMetricsModel::generateTimelineStatic(today, 
                                                      lastObservationDate,
                                                      effSpreadFreqObj);
}


/** Static method to generate the timeline on which the effective curve 
    will be calculated */
DateTimeArraySP CreditMetricsModel::generateTimelineStatic(
    const DateTime&  today,
    const DateTime&  lastObservationDate,
    MaturityPeriodSP effSpreadFreqObj)
{
    static const string method("CreditMetricsModel::generateTimelineStatic");
    try {
        // timeline offset from today      
        /* create a timeline ending strictly larger than lastObservationDate */
        DateTimeArraySP timeline(new DateTimeArray(1, today));
        DateTime tempDate(today.rollDate(1));
        do {
            tempDate = effSpreadFreqObj->toDate(tempDate);
            timeline->push_back(tempDate);
            // use getDate() to compare here to avoid having eg 1/1/05 SOD
            // AND 1/1/05 EOD in the timeline when 'effSpreadFreq' = 1d
        } while (tempDate.getDate() < lastObservationDate.getDate());
        
        /* maxObsEndDate replaces the last date in the timeline */
        if (lastObservationDate > today && lastObservationDate != tempDate) {
            timeline->back() = lastObservationDate;
        }
        return timeline;
    } 
    catch (exception& e) {
        throw ModelException(e, method, "Failed to build timeline");
    }
}


/** Returns an IFixedTrancheLossCalculator which is capable of
    returning expected tranche losses along the specified timeline
    (with the strikes as specified in product). See ConvolutionModel
    for further details. The implementation here uses
    createLossCalculator() fed into 
    ITrancheLossCalculator::createFixedTrancheLossCalculator. The
    conditionalLossCalc parameter is not set */
IFixedTrancheLossCalculator* CreditMetricsModel::createFixedLossCalculator(
    const DateTimeArray&                timeline,           /* (I) */
    CreditTrancheLossConfigConstSP      tranche,            /* (I) */
    CounterPartyCreditConstSP           cpty,               /* (I) */
    YieldCurveConstSP                   discount,           /* (I) */
    Control*                            control,            /* (I) */
    Results*                            results,            /* (I) */
    IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const // (O)
{
    try {
        // look up the strikes
        double lowStrike, highStrike;
        tranche->getTrancheStrikes(lowStrike, highStrike);

        // create loss calculator
        ITrancheLossCalculatorConstSP calculator(
            createLossCalculator(timeline, tranche, cpty));

        // and then into a fixed strikes one
        return ITrancheLossCalculator::createFixedTrancheLossCalculator(
            lowStrike, highStrike, calculator);
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::"
                             "createFixedTrancheLossCalculator(tranche)");
    }
}


/** Analogous to the loss calculator creation, but this allows us
    to model recovered notional, by pricing an inverted 1-k2, 1-k1 tranche
    instead. */
IFixedTrancheLossCalculator* CreditMetricsModel::createFixedRecoveredNotionalCalculator(
    const DateTimeArray&                timeline,           /* (I) */
    CreditTrancheLossConfigConstSP      tranche,            /* (I) */
    CounterPartyCreditConstSP           cpty,               /* (I) */
    YieldCurveConstSP                   discount,           /* (I) */
    Control*                            control,            /* (I) */
    Results*                            results,            /* (I) */
    IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const // O
{
    try {
        // look up the strikes
        double lowStrike, highStrike;
        tranche->getTrancheStrikes(lowStrike, highStrike);
        // and the total notional
        double prtfNtnl = tranche->portfolioNotional();
        // create recovered notional calculator
        ITrancheLossCalculatorConstSP calculator(
            createRecoveredNotionalCalculator(timeline, tranche, cpty));
        // and then into a fixed strikes one
        // NB we invert the strikes
        return ITrancheLossCalculator::createFixedTrancheLossCalculator(
            prtfNtnl-highStrike, prtfNtnl-lowStrike, calculator);
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::"
                             "createFixedRecoveredNotionalCalculator");
    }
}


void CreditMetricsModel::createLossCalculators(
    const DateTimeArray&                timeline,           /* (I) */
    CreditTrancheLossConfigConstSP      tranche,            /* (I) */
    CounterPartyCreditConstSP           cpty,               /* (I) */
    const DateTime&                     maturity,           /* (I) */
    Control*                            control,            /* (I) */
    Results*                            results,            /* (I) */
    bool                                recoverNotional,    /* (I) */
    YieldCurveConstSP                   discount,           /* (I) */
    IFixedTrancheLossCalculatorConstSP& lossCalculator,               /* (O) */
    IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,  /* (O) */
    IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,          /* (O) */
    IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const /* (O) */
{
    lossCalculator.reset(createFixedLossCalculator(timeline, tranche, cpty,
                                                   discount, control, results, 
                                                   conditionalLossCalc));
    if (recoverNotional) {
        recoveredNotionalCalculator.reset(createFixedRecoveredNotionalCalculator(
            timeline, tranche, cpty, discount,
            control, results, conditionalRecNtnlCalc));
    }
}

/** Creates an IEffectiveCurveGen for an 'NtD' */
ICreditLossGenSP CreditMetricsModel::createEffCurveGenerator(
    NToDefaultLossConfigConstSP ntdLossCfg,
    CounterPartyCreditConstSP   cpty,
    const bool                  recoverNotional) const
{
    static const string method("CreditMetricsModel::createEffCurveGenerator");
    throw ModelException (method, "This convolution model can NOT be used "
                          "to price NtDs.");
}


/** Indicates whether this model supports stochastic recovery rates
    AND there are any engine parameters for any names specifying so.
    By default returns false - derived types may override this method */
const bool CreditMetricsModel::hasStochasticRecoveries(
    CreditTrancheLossConfigConstSP tranche) const
{
    return false;
}


/** Recovery of notional from the top of the portfolio requires
    an additional call to the convolution.
    It only has effect if the upper strike is greater than the
    sum of (name notional * name recovery).
    Therefore we can avoid making this additional call with the
    models assistance if the product displays the correct
    characteristics. Typically this means the model must be of
    a fixed recovery type */
const bool CreditMetricsModel::modelRecoveredNotional(
    CreditTrancheLossConfigConstSP tranche) const /* (I) */
{
    // here we are interested in
    // i)  do any names call for stochastic recovery
    //     - the information of interest is really on the portfolio name
    //     - however, you could argue that the name has to have the correct
    //       engine parameters for the model...
    //     - but for simplicity we can ignore that for now, and just
    //       ask the names 
    // ii) if not, does the total recovered notional bite into the
    //     top of the tranche
    
    bool stochasticRecoveries = hasStochasticRecoveries(tranche);
    if (stochasticRecoveries) {
        return true; //we have to do the recovered notional calculation
    }
    else {
        //we only need to do the recovered notional calc
        //if sum(Ni) - sum(Ri*Ni) < upperStrike
        double sumNi = tranche->notional();

        double sumNiRi = 0.0; //this we need to calculate
        int numNames = tranche->numInnerLossConfigs();
        for (int i=0; i < numNames; ++i) {
            //we define recovered notional as notional - loss
            //where loss takes LGD parameters into account
            ICreditLossConfigConstSP name = tranche->getInnerLossConfig(i);
            sumNiRi += name->notional() - name->maxPossibleLoss();
        }

        double lowStrike;
        double highStrike;
        tranche->getTrancheStrikes(lowStrike, highStrike);

        //now compare
        return ((sumNi - sumNiRi) < highStrike);
    }
}



/** Returns an ITrancheLossCalculator which is capable of
    returning expected tranche losses along the specified timeline
    using the methodology appropriate to this model. The computeCondCurve
    indicates whether losses conditional on the counterparty
    surviving are required. The implementation here just calls 
    createLossCalculatorBase. The indirection/complexity is due to the
    rather dubious inheritance hierarchy employed.  */
ITrancheLossCalculator* CreditMetricsModel::createLossCalculator(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    return createLossCalculatorBase(timeline, tranche, cpty);
}


/** Analagous to the loss calculator above, but here we are modelling
    recovered notional. The implementation here just calls 
    createLossCalculatorBase. The indirection/complexity is due to the
    rather dubious inheritance hierarchy employed.  */
ITrancheLossCalculator* CreditMetricsModel::createRecoveredNotionalCalculator(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    return createRecoveredNotionalCalculatorBase(timeline, tranche, cpty);
}

/** Same as above createLossCalculator() method but returns one which
    is derived from CreditMetricsLossCalculatorBase. This method won't
    make sense for certain derived types - or rather its meaning changes.
    The implementation here creates an appropriate CreditMetrics loss
    calculator */
CreditMetricsLossCalculatorBase* CreditMetricsModel::createLossCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CreditMetricsFastLossCalculator(timeline, 
                                                       tranche,
                                                       false, //recoverNotional
                                                       cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CreditMetricsLossCalculator(timeline, 
                                                   tranche, 
                                                   lossUnitToUse,
                                                   false, //recoverNotional
                                                   cpty);
        }
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::createLossCalculatorBase");
    }
}


/** Analagous to createLossCalculator() above but for recovered notional */
CreditMetricsLossCalculatorBase* CreditMetricsModel::createRecoveredNotionalCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CreditMetricsFastLossCalculator(timeline, 
                                                       tranche,
                                                       true, //recoverNotional
                                                       cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CreditMetricsLossCalculator(timeline, 
                                                   tranche,
                                                   lossUnitToUse,
                                                   true, //recoverNotional
                                                   cpty);
        }
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::createLossCalculatorBase");
    }
}


/** A utility method (as inheritance hiearchy is questionable).
    Returns an ITrancheLossCalculator which is capable of returning
    expected tranche losses along the specified timeline using CCM
    algorithm.  The computeCondCurve indicates whether losses
    conditional on the counterparty surviving are required. */
CreditMetricsLossCalculatorBase* CreditMetricsModel::createCCMLossCalculator(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CCMFastLossCalculator(timeline, 
                                             tranche, 
                                             false, //recoverNotional
                                             cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CCMLossCalculator(timeline, 
                                         tranche, 
                                         lossUnitToUse,
                                         false, //recoverNotional
                                         cpty);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "CreditMetricsModel::createCCMLossCalculator");
    }
}

/** A utility method (as inheritance hiearchy is questionable).
    Returns an ITrancheLossCalculator which is capable of returning
    expected tranche losses along the specified timeline using CCM
    algorithm.  The computeCondCurve indicates whether losses
    conditional on the counterparty surviving are required. */
CreditMetricsLossCalculatorBase* CreditMetricsModel::createCCMRecoveredNotionalCalculator(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CCMFastLossCalculator(timeline, 
                                             tranche, 
                                             true, //recoverNotional
                                             cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CCMLossCalculator(timeline, 
                                         tranche, 
                                         lossUnitToUse,
                                         true, //recoverNotional
                                         cpty);
        }
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsModel::"
                             "createCCMRecoveredNotionalCalculator");
    }
}

/** Returns how the effectiveCurve should be interpolated. */
const string& CreditMetricsModel::getLossInterpolation() const{
    return lossInterpolation;
}

/** Returns all the points on the skew surface (of given name) to which this
    product is sensitive. This implementation returns null */
BetaSkewGridPointArrayConstSP CreditMetricsModel::getSensitiveBetaSkewPoints(
    OutputNameConstSP              outputName,
    const DateTime&                lastObservationDate, 
    CreditTrancheLossConfigConstSP tranche,  
    YieldCurveConstSP              discount,
    bool                           tweakAll, 
    const int                      numberOfNeighbours) const
{
    return BetaSkewGridPointArrayConstSP();
}

//------------------------------
// IHasForwardRatePricer methods
//------------------------------

/** Key method providing access to the pricer */
IForwardRatePricerSP CreditMetricsModel::getForwardRatePricer() const
{
    if (!forwardRateModel)
    {
        //not provided by the user, so supply the default
        return getDefaultForwardRatePricer();
    }
    else
    {
        //use the supplied pricer
        return forwardRateModel;
    }
}

class CreditMetricsModelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditMetricsModel, clazz);
                SUPERCLASS(CObject);
        IMPLEMENTS(IConvolutionModel);
                IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultCreditMetricsModel);
        
        FIELD        (pvToSpot, "flag to compute pv (default is FALSE)");
        FIELD_MAKE_OPTIONAL (pvToSpot);
        FIELD        (cfCutOffDate, "ignore CF before and on this date");
        FIELD_MAKE_OPTIONAL (cfCutOffDate);
        FIELD        (maxNbSlice, "max nb of slices for loss distribution");
        FIELD_MAKE_OPTIONAL (maxNbSlice);
        FIELD        (gcdDivisor, "nb by which to divide the GCD of notionals");
        FIELD_MAKE_OPTIONAL (gcdDivisor);
        FIELD        (lossUnitOverride, "Loss unit override");
        FIELD_MAKE_OPTIONAL (lossUnitOverride);
        FIELD        (lossInterpolation, "LINEAR, FLAT_FORWARD");
        FIELD        (effSpreadFreq, " date template frequency for "
                                            "effective spread curve");
        FIELD        (creditChargeViewType, "credit charge view type either "
                                                   "CONSERVATIVE or AGGRESSIVE");
        FIELD        (fastConvThreshold, "if nbName >= threshold, use fast "
                                                "convolution");
        FIELD_MAKE_OPTIONAL (fastConvThreshold);
        FIELD        (quickGreeks, "allow use of reconvolution for greeks, "
                                          "DEFAULT, YES, NO");
        FIELD_MAKE_OPTIONAL (quickGreeks);
        FIELD        (calibrationStyle, "calibration style (eg CMS) for quanto "
                                               "adjustment");
        FIELD_MAKE_OPTIONAL (calibrationStyle);
        FIELD        (calibrationMaturity, "calibration Maturity (eg 10Y) for "
                                                  "quanto adjustment");
        FIELD_MAKE_OPTIONAL (calibrationMaturity);
        FIELD_NO_DESC       (effSpreadFreqObj);
        FIELD_MAKE_TRANSIENT(effSpreadFreqObj);
        FIELD        (indexMap, "adjust single names by index basis according to the "
                                       "relationship defined in this map");
        FIELD_MAKE_OPTIONAL (indexMap);
        FIELD        (indexMapTreatment, "additional control for applying index basis, "
                                                "APPLY_TO_NONE, APPLY_TO_MAPPED_NAMES_ONLY, "
                                                "APPLY_TO_ALL_NAMES");
        FIELD_MAKE_OPTIONAL (indexMapTreatment);
        FIELD               (calculateIndexBasis, "TRUE  : calculate the basis from the index definition; "
                                                  "FALSE : retrieve named basis from the cache");
        FIELD_MAKE_OPTIONAL (calculateIndexBasis);

        FIELD               (forwardRateModel, "A model capable of pricing all fees");
        FIELD_MAKE_OPTIONAL (forwardRateModel);
    }
    static IObject* defaultCreditMetricsModel(){
        return new CreditMetricsModel(CreditMetricsModel::TYPE);
    }
};

CClassConstSP const CreditMetricsModel::TYPE = CClass::registerClassLoadMethod(
    "CreditMetricsModel", typeid(CreditMetricsModel), CreditMetricsModelHelper::load);

const string CreditMetricsModel::QG_DEFAULT="DEFAULT";
const string CreditMetricsModel::QG_YES="YES";
const string CreditMetricsModel::QG_NO="NO";

bool CreditMetricsModelLoad() {
    return (CreditMetricsModel::TYPE != 0);
}

DRLIB_END_NAMESPACE



