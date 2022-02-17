//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RainbowRangeKO.cpp
//
//   Description : Payoff is same as RangeNote, but sampling is done on worst/best performer of basket.
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/IRebate.hpp"
#include "edginc/RangeAccrue.hpp"
#include "edginc/KOStubRule.hpp"
#include "edginc/KOFixedLeg.hpp"
#include "edginc/KOLiborLeg.hpp"
#include "edginc/Events.hpp"

DRLIB_BEGIN_NAMESPACE
class RainbowRangeKO: public GenericNFBase, 
                 virtual public IMCIntoProduct,
                 virtual public LegalTerms::Shift,
                 virtual public BarrierBreach::IEventHandler{
public:
    static CClassConstSP const TYPE;
    friend class RainbowRangeKOMC;
    
    // validation
    void validatePop2Object(){
        static const string method = "RainbowRangeKO::validatePop2Object";
        GenericNFBase::validatePop2Object();

        try
        {
            // barrier can be only "N".
            if (barrier->getInterp() != "N") {
                throw ModelException(method,
                                     "Linear or Step barrier schedule is not allowed. Only 'N' is available.");
            }

            if (hasFloater && !floater.get()){
                throw ModelException(method,
                                     "hasFloater = true, but floater is not found");
            }

            // Make monitoringDate.
            DateTimeArray monitoringDates = DateTime::merge(rangeAccrueMaker->getMonitorDates(), 
                                                            barrier->getDateArray());             

            // KO date should be subset of range accrue date.
            if (!DateTime::isSubset(rangeAccrueMaker->getMonitorDates(), barrier->getDateArray())) {
                throw ModelException(method,
                                     "Require barrier schedule to be a subset of range monitor dates");
            }

            // request always averageOutDates
            if (averageOutDates.size() <=0)
                throw ModelException(method, "Require at least one date for averageOutDates");

            // painful but too much code implicitly depending on this in Barrier class
            // require all simulation dates are monitoring dates - so ...
            if (!DateTime::isSubset(monitoringDates, averageOutDates)) {
                throw ModelException(method,
                                     "Require averageOutDates to be a subset of barrier and range monitor dates");
            }

            //now validate the economic barrier if it exists
            if (economicLevels.get()) {
                // barrier can be only "N".
                if (economicLevels->getInterp() != "N") {
                    throw ModelException(method,
                                        "Linear or Step economic barrier schedule is not allowed. Only 'N' is available.");
                }

                // KO date should be subset of range accrue date.
                if (!DateTime::isSubset(rangeAccrueMaker->getMonitorDates(), economicLevels->getDateArray())) {
                    throw ModelException(method,
                                        "Require economic barrier schedule to be a subset of range monitor dates");
                }
            }
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        DateTimeArray monitoringDates = DateTime::merge(rangeAccrueMaker->getMonitorDates(), 
                                                        barrier->getDateArray());
        return DateTime::merge(averageOutDates, monitoringDates);
    }


    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    // BarrierBreach::IEventHandler interface
    virtual void getEvents(const BarrierBreach* breach, IModel* model, 
        const DateTime& eventDate, EventResults* events) const;

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
        // parent
        GenericNFBase::GetMarket(model, market);
        // relevant elements of self
        rebate->getMarket(model, market.get());
        if(hasFloater)
        {
            floater->setCouponCurve(discount); // in case floater gets coupon curve from inst
            floater->getMarket(model, market.get());
        }
    }

    // set barriers to be the economic (legal) ones
    bool sensShift(LegalTerms* shift) {
        // just replace all barriers with the corresponding economic one
        economicLevels = barrier;

        return true; // continue shifting
    }

protected:
    /// fields 
    IDoubleArrayModifierMakerSP  matPerformance;    // option at maturity
    IAggregateMakerSP            assetBasket;       // pct/rainbow/product basket etc of ...
    DateTimeArray                averageOutDates;   // purely for getting average perf of each asset, XXX what about extremes?
    string                       assetMeasurement;  // "Avg", "Min", "Max"
    IRebateMakerSP               rebate;            // the "opposite" side of the payoff from 'assetBasket'
    bool                         hasFloater;    
    
    RangeAccrueMakerSP           rangeAccrueMaker;  //interface
    
    BarrierUnionSP               barrierUnion;      // Barrier Union Class, only for internal use.  Not public interface. $unregistered
    
    ScheduleSP                   barrier;           // KO barrier.  Sampling on (plain/rainbow) basket level.
    ScheduleSP                   economicLevels;    // economical barrier 
    bool                         isUp;              // isUp barrier?    
    bool                         isOut;             // "True=>Knock Out; False=> In");
    bool                         addCpnOnKODate;    // (false) sum up Coupon before KO Date. (true) add Coupon on KO Date.
    DateTime                     breachDate;        // KO date if it's already breached.

    double                       redemption;        // redemption amount at maturity (only maturity.  If it's KO, it is not paid).
    LiborLegSP                   floater;
    string                       koStubRule;        // ko Stub Rule for Libor

private:
    RainbowRangeKO(): GenericNFBase(TYPE), assetMeasurement("Avg") {
    koStubRule = "N";
    } // for reflection
    RainbowRangeKO(const RainbowRangeKO& rhs);     // not implemented
    RainbowRangeKO& operator=(const RainbowRangeKO& rhs); // not implemented

    static IObject* defaultRainbowRangeKO(){
        return new RainbowRangeKO();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RainbowRangeKO, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultRainbowRangeKO);
        FIELD(matPerformance, "Overall Perf Modifier");
        FIELD(assetBasket, "How to aggregate asset perfs");
        FIELD(averageOutDates, "Averaging out dates for basket");
        FIELD(assetMeasurement, "Avg/Min/Max");
        FIELD_MAKE_OPTIONAL(assetMeasurement);
        FIELD(rangeAccrueMaker, "Range accrual");
        FIELD(barrier, "KO barrier.  Sampling on (plain/rainbow) basket level");
        FIELD(redemption, "redemption amount at maturity. (only maturity.  If it's KO, it is not paid.)");
        FIELD(isUp, "is Up barrier?");    
        FIELD(isOut, "True=>Knock Out; False=> In");
        FIELD(addCpnOnKODate, "[false] sum up Coupon before KO Date. [true] add Coupon on KO Date.");
        FIELD(economicLevels, "Economic barrier schedule");
        FIELD_MAKE_OPTIONAL(economicLevels);
        FIELD(breachDate, "breached date is necessary when it's breached.");
        FIELD(rebate, "Opposite side of payoff to assetBasket");
        FIELD_MAKE_OPTIONAL(rebate);
        FIELD(hasFloater, "has Libor Leg?")
        FIELD(floater, "libor leg");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(koStubRule, "KO Stub rule for swap leg stream.  Active only whan hasFloater is true.");
        FIELD_MAKE_OPTIONAL(koStubRule);
    }
};

/* MC product class for RainbowRangeKO */
class RainbowRangeKOMC: public IMCProduct,
                   virtual public IMCProductLN,
                   virtual public IMCProductImplied{
private:
    const RainbowRangeKO*          inst;
    BarrierScheduleSP           barSched;       // non-const note
    DateTimeArray               monitoringDates;    
    DateTimeArray               payDates;       // range pay dates

    int                       nbAssets;         // convenient

    SimpleDoubleArray         assetComps;
    SimpleDoubleArray         assetCompsForRange;
    IAggregateSP              assetBasket;
    IAggregateSP              assetBasketForRange;
    TrivialDoubleArray        basket;  // just a double ... being the aggregation of per-asset "perfs"
    IRebateSP                 rebate;
    IDoubleArrayModifierSP    overallOption;

    int                       assetMeasure; // 0 - Avg, 1 - Min, 2 - Max
    IntArray                  avgMap;    // to track averaging dates
    DoubleArray               sum;       // [nbAssets], saves alloc later
    DoubleArray               fvFactors; // if rebate at hit the this FVs to mat. 
                                         // Past rebates drop out via fvFactor=0.0
    double                    df;    // pvFromPaymentDate

    // for past
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 

    RangeAccrueSP            rangeAccrue;
    double                   redemption;    
    LiborLegSP               floater;
    KOLiborLegSP             koFlt;         // for internal usage
    // the following member is necessary to be member of RRKMC member to be "alive" in entire object.
    // formally, I temporary made this variable at constructor,  but koFlt's info could be deleted.
    KOStubRuleSP             koRule;        
    
    DoubleArray              liborValues;   // libor value when it's knocked out at each time step.
    double                   plainLibor;    // plain libor value
    bool                     hasFloater;    
    string                   koStubRule;    

    CashFlowSP               knownRebate;   // if already touched, for known cash flow

public:

    /** equivalent to  InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to RainbowRangeKO) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    RainbowRangeKOMC(const RainbowRangeKO*          inst,
                const SimSeriesSP&        simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        assetComps(nbAssets, 0.0),
        assetCompsForRange(nbAssets, 0.0),
        basket(0.0),
        rebate(inst->rebate->getRebate(simSeries->getAllDates(), discount)),
        sum(nbAssets, 0.0),
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0){

        // -----------  Create Monitoring Array ----------//
        // make the monitoringDates by merging barrier and rangeAccrue's monitoring daes
        // However, souceBarDates are assumed (validate) to be subset of rangeMonitorDates.
        // Thus, for time being, there is not much meaning to make monitoringDate by merge.
        DateTimeArray sourceBarDates = inst->barrier->getDates();
        DateTimeArray rangeMonitorDates = inst->rangeAccrueMaker->getMonitorDates();
        DoubleArray sourceBarValues = inst->barrier->getValues();
        monitoringDates = DateTime::merge(rangeMonitorDates, sourceBarDates);
        int nbSteps = monitoringDates.size();

        // -----------  Set up internal BarrierSchedule Class ----------//
        // make internal Barrier Schedule class from one factor barrier interface
        // Need add barrier schedule date and dummmy level, 
        // so as to avoid interpolation problem at any monitoring date.
        int i=0,j=0;
        double extremeValue = inst->isUp ? 1.0/1.0e-10 : -1.0e-10;
        DoubleArray barLevels = DoubleArray(nbSteps, extremeValue);
        for (i=0;i<nbSteps;i++){
            if(monitoringDates[i] == sourceBarDates[j]){
                barLevels[i] = sourceBarValues[j];
                j++;
                if (j>=sourceBarDates.size())
                    break;  
            }
        }
        ScheduleSP bar_along_monDates = ScheduleSP(new Schedule(monitoringDates,barLevels,inst->barrier->getInterp()));

        // do the same for the economic barrier otherwise it might fail in the past
        // for the economic barrier I'm assuming dates are subset of monitoring dates 
        // and I'm not merging otherwise this will open a can of worms
        // XXX this needs revisiting XXX
        bool hasEconomicBarrier = inst->economicLevels.get() && 
                                    !(inst->economicLevels->getDates().empty()) ? 
                                        true : false;
        DateTimeArray sourceEconBarDates = hasEconomicBarrier ? 
                                            inst->economicLevels->getDates() :
                                            inst->barrier->getDates();
        DoubleArray sourceEconBarValues = hasEconomicBarrier ? 
                                            inst->economicLevels->getValues() :
                                            inst->barrier->getValues();
        DoubleArray econLevels = DoubleArray(nbSteps, extremeValue);
        j=0;
        for (i=0;i<nbSteps;i++){
            if(monitoringDates[i] == sourceEconBarDates[j]){
                econLevels[i] = sourceEconBarValues[j];
                j++;
                if (j>=sourceEconBarDates.size())
                    break;  
            }
        }
        ScheduleSP econBar_along_monDates = ScheduleSP(new Schedule(monitoringDates,econLevels,inst->barrier->getInterp()));

        barSched = BarrierScheduleSP(new BarrierSchedule(bar_along_monDates,
                                 inst->isOut,
                                 inst->isUp,
                                 econBar_along_monDates,
                                 nbAssets,                 // numHits is ignored.
                                 IntArraySP(new IntArray(nbAssets, 0)),
                                 inst->breachDate,
                                 "E",                      // only European monitoring allowed.
                                 false,                    // no smoothing allowed
                                 0.0,                   
                                 0.0));
        barSched->validate(nbAssets);
        // create an interpolated barrier to match the monitoring dates
        barSched->createInterpBarrier(inst->valueDate, monitoringDates);


        // -----------  Set Up Basket for Final Option / Range Sampling ----------//
        // Tie the pieces together :
        // First creating the per-asset performances, then aggregating them into a single number and
        // finally turning that into an overall perf
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
        assetBasketForRange = IAggregateSP(inst->assetBasket->getAggregate(&assetCompsForRange));
        overallOption = IDoubleArrayModifierSP(inst->matPerformance->getModifier(&basket));

        // -----------  Make Range Accrue class ----------//
        // 1. prepare the PVed amount array
        // 2. Let RangeClass know which step is range monitoring step.
        rangeAccrue = RangeAccrueSP(inst->rangeAccrueMaker->getDiscounted(inst->valueDate, discount));
        rangeAccrue->setIsMonitorStep(simSeries->getAllDates());
        
        
        // -----------  Set up Product.... ----------//        
        df =  pvFromPaymentDate();                    
        redemption = inst->redemption;

        // Account for inst settlement, but forbid fixed settlement date
        if (CashSettleDate::TYPE->isInstance(settlement.get())) {
            throw ModelException("RainbowRangeKOMC","Rebate pay at hit does not allow instrument settlement of CashSettleDate!");
        }

        
        // -----------  Set up pvFactors for Rebate ----------//
        const DateTime& today = getToday();
        payDates = rangeAccrue->getPayDates();
        fvFactors = DoubleArray(nbSteps, 1.0); // default is pay at mat
        for(i=0; i<fvFactors.size(); i++) {
            // payDate is not detemined by settlement but by RangeAccrue Paydates
            // Flat and Schedule Rebate has no problem by this treatment, 
            // but I don't know it's correct for streamRebate.
            if (payDates[i] <= today) {
                // past so already paid
                fvFactors[i] = 0.0;
            } else {
                // infrastructure will pv from payment date, so override for this case
                fvFactors[i] = discount->pv(payDates[i]) / df; 
            }
        }

        // -----------  Set up Libor Leg ----------//
        // I'd like to move many parts of those to KOFixedLeg.
        liborValues = DoubleArray(nbSteps,0.0);     //initialize by zero value.
        koStubRule = inst->koStubRule;    
        hasFloater = inst->hasFloater;
        if (hasFloater){
            floater = inst->floater;
            DateTimeArray liborPayDates = floater->PayDates;

            // set KO settlement info
            // "N"one :  all future payment coupon will be cancelled. (Nothing to do)
            // "S"wap :  coupon is accrued up to next payment date in range coupon schedule.
            // "B"ond :  Already accrue started coupon is fully paid at next (libor's) payment date.  
            bool isAccrueUpToSettle = true;        // 
            DateTimeArray settleDates = payDates;   // give the Range Accrue Pay Dates.
            DateTimeArray lastMonDates = DateTimeArray(0);
            DateTimeArray AccrueDates = DateTimeArray(0);
            if (koStubRule == "S") {
                AccrueDates.push_back(rangeMonitorDates[0]);
                for (i=1;i<settleDates.size();i++){
                    if (settleDates[i-1] < settleDates[i]){                                   
                        lastMonDates.push_back(rangeMonitorDates[i-1]);
                        AccrueDates.push_back(rangeMonitorDates[i]);
                    }                
                }
                lastMonDates.push_back(rangeMonitorDates[i-1]);        // And maturity
                if (rangeMonitorDates[i-1] == AccrueDates[AccrueDates.size()-1])
                    AccrueDates.push_back(rangeMonitorDates[i-1].rollDate(1));  // dummy Accrue End.
                else
                    AccrueDates.push_back(rangeMonitorDates[i-1]);
                DateTime::removeDuplicates(settleDates, false);     //remove the duplicate.
                if (lastMonDates.size() != settleDates.size())
                    throw ModelException("RainbowRangeKOMC","libor pay date size is not equal to lastMonDates size.");
            }
            else{
                settleDates = liborPayDates;    // Thus, those schedule are not used.
                lastMonDates = liborPayDates;   // Thus, those schedule are not used.
                AccrueDates = floater->AccrualDates;
                isAccrueUpToSettle = false;
            }

            string setlleTiming = "AsSchedule";
            koRule = KOStubRuleSP(new KOStubRule(koStubRule,isAccrueUpToSettle,setlleTiming));    //koStubRule,isAccrueUpToSettle,payTimingRule
            koFlt = KOLiborLegSP(new KOLiborLeg(floater,
                                                koRule,
                                                inst->instSettle,
                                                settleDates,
                                                AccrueDates,
                                                inst->valueDate,
                                                inst->discount.get()));
            
            // get plain libor leg value and value at each time line
            // values are converted as of last settlement date.
            plainLibor = koFlt->getPV(inst->valueDate, discount)/df;
            liborValues = koFlt->getPVsAlongKOTimeStep(monitoringDates,1.0/df);      

            // prepare KNOWN_CASH_FLOW.  plain libor leg at first level.
            koFlt->makeKnownCashFlows(inst->valueDate, false);
        }


        // -----------  prepare for Average-In ----------//
        bool isTrivial;
        avgMap = DateTime::createMapping(simSeries->getAllDates(),
                                         inst->averageOutDates,
                                         isTrivial);

#define DO_AVG 0
#define DO_MIN 1
#define DO_MAX 2
        if (inst->assetMeasurement == "Min") {
            assetMeasure = DO_MIN;
            for(int j=0; j<nbAssets; j++) {
                sumSoFar[j] = 1e30;
            }
        } else if (inst->assetMeasurement == "Max") {
            assetMeasure = DO_MAX;
        } else {
            assetMeasure = DO_AVG;
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        static const string method = "RainbowRangeKO::payoff";
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset, iStep;

        // 1. determine raw performances for each asset
        // 2. establish KO multiplier ("prob") for the basket (to be)
        // 3. (new). calculate Range Acrrue Payment..
        // 4. aggregate into a basket
        // 5. apply any modifier (call/put/etc) to basket
        // 6. form "expected payoff" combining overall option and rebate
        sum = sumSoFar;
        
        // 1a. determine raw performances for each asset
        for (iStep=beginIdx+avgMap[beginIdx]; iStep<endIdx; 
             iStep++, iStep += avgMap[iStep]) {
            if (assetMeasure==DO_MIN) {
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    sum[iAsset] = Maths::min(pathGen->Path(iAsset, 0)[iStep], sum[iAsset]);
                }
            } else if (assetMeasure==DO_MAX) {
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    sum[iAsset] = Maths::max(pathGen->Path(iAsset, 0)[iStep], sum[iAsset]);
                }
            } else {
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];
                }
            }
        }

        // 2. establish KO multiplier
        // Note this is outside the "if past" clause, which allows
        // the barrier to handle its own past.
        double koFactor;
        bool isConditionMet;
        int metAtStep;
        barSched->hitValueAndTimeAsBasket(pathGen,
                                 inst->assetBasket,
                                 koFactor,
                                 isConditionMet,
                                 metAtStep);

        int atStep = isConditionMet ? metAtStep : endIdx;
        atStep += inst->addCpnOnKODate ? 1 : 0;
        atStep = atStep > endIdx ? endIdx : atStep; //capped so as to avoid ABR error.

        if (pathGen->doingPast()){ 
            sumSoFar = sum;
            // This for BARRIER_LEVEL request. See RainbowKO.cpp for more complete comment
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }
            // prepare the KNOWN_CASH_FLOW by looking at the past values.
            rangeAccrue->makeKnownCashFlow(pathGen, atStep, hasFuture(), assetBasketForRange, assetCompsForRange);
            // re-set libor KNOWN_CASH_FLOW when it's already breached.  
            if (hasFloater && isConditionMet)
                koFlt->makeKnownCashFlows(monitoringDates[metAtStep], isConditionMet);
            if (isConditionMet)
                knownRebate = CashFlowSP(new CashFlow(payDates[metAtStep],rebate->getLevel(metAtStep)));
        }


        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().

            // 1b. finished this job
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                double avg = sum[iAsset];
                if (assetMeasure==DO_AVG) {
                    avg /= inst->averageOutDates.size();
                }
                assetComps[iAsset] = avg/pathGen->refLevel(iAsset, 0);
            }

            // 3(new). calculate range coupon.  
            // range accrue class has value as of today.
            double rangePay = rangeAccrue->getValue(pathGen, atStep, assetBasketForRange, assetCompsForRange)/df;
            
            // 4. aggregate into a basket
            basket() = assetBasket->aggregate();

            // 5. apply any modifier (call/put/etc) to basket
            overallOption->apply();

            // 6. form "expected payoff" combining overall option and rebate
            // Note if condition NOT met then rebate level at mat (last step) is used.
            double rebateValue = isConditionMet? 
                (fvFactors[metAtStep] * rebate->getLevel(metAtStep)) : 
                rebate->getLevel(endIdx-1);
            double payoff = koFactor * (basket() + redemption) + (1.0 - koFactor) * rebateValue + rangePay;
            if (hasFloater)
                payoff +=  isConditionMet ? liborValues[metAtStep] : plainLibor ;
            
            // Done...
            prices.add(inst->notional * payoff); 
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        barSched->adjustLN(pathGen,
                          this,
                          getMultiFactors(),
                          0);
    }

    /** Use this opportunity to do any Implied driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustments */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{
        barSched->adjustImplied(pathGen,
                               this,
                               getMultiFactors(),
                               inst->discount.get(),
                               inst->instSettle.get(),
                               0);
    }

    // for the LogNormal path generator
    // XXX since there could be several interp levels needed we should
    // XXX defer construction of reqarr to a method on the barrier. TBD
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "RainbowRangeKO::getVolInterp";

        try
        {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime&    today = getToday();
            const DateTime&    startDate = getRefLevel()->getAllDates().front();
            bool               fwdStarting = startDate.isGreater(today);
            const DateTime&    lastSimDate = getSimSeries()->getLastDate();
            double interpLevel = (fwdStarting? 1.0 : pathGen->refLevel(iAsset, 0));  // to do review

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
            return reqarr;
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    /** For now am using this instead of implementing the IHandlePaymentEvents interface
        since I'd like to reuse the PAYMENT_DATES and KNOWN_CASHFLOWS features from
        the IMCProduct, but need to add BARRIER_LEVEL support. Should review. XXX */
    void recordExtraOutput(CControl*     control,
                           Results*      results,
                           const IMCPrices& prices) const {
        static const string method = "RainbowRangeKO::recordExtraOutput";
        try{
            CashFlowArraySP knownCfl = CashFlowArraySP(new CashFlowArray(0));

            // KNOWN CASH FLOW
            // not yed implemented after Maturity...
            OutputRequest* request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished())
            {   // add plain swap as known cash flow.  (Although, could be cancelled....) 
                CashFlowArraySP rng = rangeAccrue->getKnownCashFlows();
                if (hasFloater)
                    knownCfl = CashFlow::merge(rng,koFlt->getKnownCashFlows());
                else 
                    knownCfl = rng;
                if (!!knownRebate)
                    knownCfl->append(knownRebate);                
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        knownCfl.get());   
            }

            // PAYMENT_DATES  
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                // in the known cash flow (i.e. libor + known Range Coupon)
                DateTimeArray cashPayDates = DateTimeArray(0);
                for (int i=0;i<(*knownCfl).size();i++){
                    cashPayDates.push_back((*knownCfl)[i].date);
                }
                DateTimeArray listPayDates = rangeAccrue->getPayDates();
                DateTime::removeDuplicates(listPayDates,false);
                DateTimeArray payDateArray = DateTime::merge(listPayDates,cashPayDates);
                OutputRequestUtil::recordPaymentDates(control,results,&payDateArray);
            }

            // BARRIER_LEVEL ...
            // cannot report the level of basket.....
            // currently, calc the barrier level by assuming one asset.  Not good as
            // no necessary Barrier Breach report would be generated.
            request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
            const DateTime& today = getToday();
            const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
            // Only try to satisfy this request when have some past (so easy to get ref levels).
            if (request && (today >= lastRefDate)) {
                int iAsset;
                // This operates on a per-asset basis
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    // Mostly delegate to Barrier class...
                    // histRefLevels are to allow absolute barriers to be reported.
                    BarrierLevelArraySP levels = 
                        barSched->reportLevels(today, histRefLevels[iAsset], iAsset);
                    if (!levels->empty()) {
                        OutputRequestUtil::recordBarrierLevels(
                            control, results,
                            getMultiFactors()->assetGetTrueName(iAsset),
                            levels.get());
                    }
                }
            }

        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    // barrier breach event based on this path
    virtual void getEvents(const IMCPathGenerator*  pathGen, 
        EventResults* events,
        const DateTime& eventDate) {

        static const string method = "RainbowRangeKOMC::getEvents";
        try {
            barSched->getEvents(pathGen,events,"Rainbow Range KO Barrier", inst->assetBasket);
        }catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RainbowRangeKO::createProduct(const MonteCarlo* model) const {    
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(samplingDates());

    return new RainbowRangeKOMC(this, simSeries);
}

// implementation of Events interface
void RainbowRangeKO::getEvents(const BarrierBreach* breach,
    IModel* model, 
    const DateTime& eventDate,
    EventResults* events) const {
    static const string method = "RainbowRangeKO::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));
            MCPathGeneratorSP pastPathGenerator(
                mc->getPathConfig()->pastPathGenerator(prod.get()));
            prod->getEvents(pastPathGenerator.get(), events, eventDate);
        } else {
            throw ModelException(method, 
                "Internal error - expected Monte Carlo model for RRK pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


CClassConstSP const RainbowRangeKO::TYPE = CClass::registerClassLoadMethod(
    "RainbowRangeKO", typeid(RainbowRangeKO), RainbowRangeKO::load);

// * for class loading (avoid having header file) */
bool RainbowRangeKOLoad() {
    return (RainbowRangeKO::TYPE != 0);
}

DRLIB_END_NAMESPACE
