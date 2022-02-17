//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RainbowDKO.cpp
//
//   Description : RainbowDKO except with DoubleBarrier
//
//   Date        : Oct 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/IRebate.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/Events.hpp"

DRLIB_BEGIN_NAMESPACE
class RainbowDKO: public GenericNFBase, 
                  virtual public IMCIntoProduct,
                  virtual public BarrierBreach::IEventHandler{
protected:
    /// fields 
    IDoubleArrayModifierMakerSP  overallOption; // call/put/digital etc on sum of ...
    IAggregateMakerSP            assetBasket;   // pct/rainbow/product basket etc of ...
    IDoubleArrayModifierMakerSP  assetBasketComponents;    // call/put/digital etc on each asset
    DateTimeArray                averageOutDates; // purely for getting average perf of each asset, XXX what about extremes?
    string                       assetMeasurement;  // "Avg", "Min", "Max"
    IRebateMakerSP               rebate;        // the "opposite" side of the payoff from 'assetBasket'

    BarrierUnionSP               barrierUnion;  // determines KO "prob" -> hitValue()

    // floating/fix leg linked to barrierUnion
    bool                         streamIsReverse;    // true if use barrierUnion to KO stream, false to KI (default)
    LiborLegSP                   floater;
    string                       floatStubRule;
    FixedLegSP                   fixedLeg;
    string                       fixedStubRule;

    // allow for 1KI and 1KO
    bool                         useKiTrigger;  // true if use triggerUnion as KI and barrierUnion has to be KO
    BarrierUnionSP               triggerUnion; 

public:
    static CClassConstSP const TYPE;
    friend class RainbowDKOMC;

    // validation
    void validatePop2Object(){
        static const string method = "RainbowDKO::validatePop2Object";
        GenericNFBase::validatePop2Object();

        try
        {
            (barrierUnion->getBarrier())->validate(assets->NbAssets());

            // some processing for 1KI/1KO barrier 
            if( useKiTrigger )
            {
                (triggerUnion->getBarrier())->validate(assets->NbAssets());

                // must be barrierSchedule and correct KO/KI field
                BarrierSchedule *ko = dynamic_cast<BarrierSchedule*>(barrierUnion->getBarrier());
                BarrierSchedule *ki = dynamic_cast<BarrierSchedule*>(triggerUnion->getBarrier());
                if( !ko || !ki || !ko->isKO() || ki->isKO() )
                    throw ModelException("RainbowKOMC", "If useKiTrigger, triggerUnion/barrierUnion must be KI/KO barrierSchedule respectively");
            }

            // if stream rebate, does not allow payAtMaturity nor non-empty paydates
            if( StreamRebateMaker::TYPE->isInstance(rebate.get()) )
                    throw ModelException("RainbowDKOMC", "Can not be stream rebate");

        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
        // parent
        GenericNFBase::GetMarket(model, market);
        // relevant elements of self
        if( !!rebate )
            rebate->getMarket(model, market.get());
        if( !!floater && floater->getSize() > 0 )
        {
//            floater->setCouponCurve(discount); // in case floater gets coupon curve from inst
            floater->getMarket(model, market.get());
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        // stripped from createProduct - all rather messy and code is duplicated
        // there must be a better way. 
        BarrierSP realBarrier = BarrierSP(copy(barrierUnion->getBarrier()));
        if( useKiTrigger ) {
            realBarrier = BarrierSP(new DoubleBarrier(true/*kiKeepKO*/, false/*simultaneous*/, realBarrier, 
                                                    BarrierSP(copy(triggerUnion->getBarrier()))));
        }

        realBarrier->createBarrierUtil(valueDate);
        bool hasContinuous;
        DateTimeArray monitoringDates;
        realBarrier->getBarrierUtil()->relevantDates(monitoringDates, hasContinuous);

        return DateTime::merge(averageOutDates, monitoringDates);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below
    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach, IModel* model, 
        const DateTime& eventDate, EventResults* events) const;

private:
    RainbowDKO(): GenericNFBase(TYPE), assetMeasurement("Avg"), streamIsReverse(false),  
                  floatStubRule("N"), fixedStubRule("N"), useKiTrigger(false) {} // for reflection
    RainbowDKO(const RainbowDKO& rhs);     // not implemented
    RainbowDKO& operator=(const RainbowDKO& rhs); // not implemented

    static IObject* defaultRainbowDKO(){
        return new RainbowDKO();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RainbowDKO, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultRainbowDKO);
        FIELD(overallOption, "Overall Perf Modifier");
        FIELD(assetBasket, "How to aggregate asset perfs");
        FIELD(assetBasketComponents, "Asset Perf Modifier");
        FIELD(averageOutDates, "Averaging out dates for basket");
        FIELD(assetMeasurement, "Avg/Min/Max");
        FIELD_MAKE_OPTIONAL(assetMeasurement);
        FIELD(rebate, "Opposite side of payoff to assetBasket");
        FIELD_MAKE_OPTIONAL(rebate);
        FIELD(barrierUnion, "Barrier Union");
        FIELD(useKiTrigger,  "true if use triggerUnion for KI and barrierUnion has to be KO");
        FIELD_MAKE_OPTIONAL(useKiTrigger);
        FIELD(triggerUnion,         "ignored if useTriggerKI is false");
        FIELD_MAKE_OPTIONAL(triggerUnion);
        FIELD(streamIsReverse, "true if KO/KI stream when barrierUnion KO/KI, false to KI/KO stream when barrierUnion KO/KI (default)");
        FIELD_MAKE_OPTIONAL(streamIsReverse);
        FIELD(floater,              "floater linked to barrierUnion");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(floatStubRule, "N(one), B(ond) or S(wap)");
        FIELD_MAKE_OPTIONAL(floatStubRule);
        FIELD(fixedLeg,             "fixedLeg linked to barrierUnion");
        FIELD_MAKE_OPTIONAL(fixedLeg);
        FIELD(fixedStubRule, "N(one), B(ond) or S(wap)");
        FIELD_MAKE_OPTIONAL(fixedStubRule);
    }
};

/* MC product class for RainbowDKO */
class RainbowDKOMC: public IMCProduct,
                   virtual public IMCProductLN,
                   virtual public IMCProductImplied{
private:
    const RainbowDKO*          inst;
    BarrierSP                 barrier;   // non-const note
    int                       nbAssets;  // convenient

    SimpleDoubleArray         assetComps;
    IDoubleArrayModifierSP    basketComponents;
    IAggregateSP              assetBasket;
    TrivialDoubleArray        basket;  // just a double ... being the aggregation of per-asset "perfs"
    IDoubleArrayModifierSP    overallOption;

    int                       assetMeasure; // 0 - Avg, 1 - Min, 2 - Max
    IntArray                  avgMap;    // to track averaging dates
    DoubleArray               sum;       // [nbAssets], saves alloc later

    // for past
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 

    // IBarrier elements
    mutable IBarrierMgrSP       barrierStat;
    mutable IBarrierMgrSP       rebateMgr;
    mutable IBarrierMgrSP       floatMgr;
    mutable IBarrierMgrSP       fixedMgr;
    double                      matFV;
    double                      payoffHist; // value from historical observation

    CashFlowArraySP             knownCFs;

public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to RainbowDKO) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    RainbowDKOMC(const RainbowDKO*          inst,
                const SimSeriesSP&        simSeries,
                BarrierSP                 realBarrier,
                BarrierSP                 streamBarrier):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst),
        barrier(realBarrier),
        nbAssets(getNumAssets()),
        assetComps(nbAssets, 0.0),
        basket(0.0),
        sum(nbAssets, 0.0),
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0),
        payoffHist(0.0),
        knownCFs(CashFlowArraySP(new CashFlowArray(0))) {

        static const string method("RainbowDKOMC");

        // overwrite paymentDate, because FixedLeg or EqLink could have paydates which
        // is longer then instSettle.
        if( !!inst->floater && inst->floater->getSize() > 0 )
        {
            DateTime lastPayDate = inst->floater->PayDates.back();
            if (lastPayDate > paymentDate)
                paymentDate = lastPayDate;
        }
        if( !!inst->fixedLeg && inst->fixedLeg->getSize() > 0 )
        {
            DateTime lastPayDate = inst->fixedLeg->PaymentDatesArray.back();
            if (lastPayDate > paymentDate)
                paymentDate = lastPayDate;
        }

        
        // Tie the pieces together :
        // First creating the per-asset performances, then aggregating them into a single number and
        // finally turning that into an overall perf
        basketComponents = IDoubleArrayModifierSP(inst->assetBasketComponents->getModifier(&assetComps));
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
        overallOption = IDoubleArrayModifierSP(inst->overallOption->getModifier(&basket));

        bool isTrivial;
        avgMap = DateTime::createMapping(simSeries->getAllDates(),
                                         inst->averageOutDates,
                                         isTrivial);

        matFV = (paymentDate<=getToday())?0.0:(1.0/pvFromPaymentDate());

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

        // create objects for BarrierUtil
        DateTimeArraySP simDates(new DateTimeArray(simSeries->getAllDates()));
        DateTimeArray obsDates(1, simDates->back());
        try {
            barrierStat = BarrierFactory::makeManager(barrier->getBarrierUtil(), BarrierFactory::makeTrivialPay(obsDates));
            barrierStat->preprocess(inst->valueDate, simDates, discount);
        } catch ( exception& e ) {
            throw ModelException(e, method, "Can not create maturity option barrier pay");
        }

        // create objects for rebate
        if( !!inst->rebate )
        {
            try {
                BarrierSchedule *barSchd = dynamic_cast<BarrierSchedule*>(streamBarrier.get());
                if( !barSchd )
                    throw ModelException(method, "Stream barrier must be BarrierSchedule in order to create rebate barrier pay");
                
                rebateMgr = BarrierFactory::makeManager(streamBarrier->getBarrierUtil(), 
                    IRebateMaker::createBarrierPay(inst->rebate.get(), barSchd->isKO(), settlement.get()));
                rebateMgr->preprocess(inst->valueDate, simDates, discount);
            } catch ( exception& e ) {
                throw ModelException(e, method, "Can not create rebate barrier pay");
            }
        }

        // create floater/fixed leg
        if( !!inst->floater && inst->floater->getSize() > 0 )
        {
            try {
                floatMgr = BarrierFactory::makeManager(streamBarrier->getBarrierUtil(), 
                                inst->floater->createBarrierPay(inst->streamIsReverse, inst->floatStubRule));
                floatMgr->preprocess(inst->valueDate, simDates, discount);
            } catch ( exception& e ) {
                throw ModelException(e, method, "Can not create float leg barrier pay");
            }
        }
        if( !!inst->fixedLeg && inst->fixedLeg->getSize() > 0 )
        {
            try {
                fixedMgr = BarrierFactory::makeManager(streamBarrier->getBarrierUtil(), 
                                inst->fixedLeg->createBarrierPay(inst->streamIsReverse, inst->fixedStubRule));
                fixedMgr->preprocess(inst->valueDate, simDates, discount);
            } catch ( exception& e ) {
                throw ModelException(e, method, "Can not create fixed leg barrier pay");
            }
        }

    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        static const string method = "RainbowDKO::payoff";
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset;

        // 1. determine raw performances for each asset
        // 2. establish KO multiplier ("prob") for the basket (to be)
        // 3. apply any performance modifiers (call/put/etc)
        // 4. aggregate into a basket
        // 5. apply any modifier (call/put/etc) to basket
        // 6. form "expected payoff" combining overall option and rebate
        sum = sumSoFar;
        
        // 1a. determine raw performances for each asset
        for (int iStep=beginIdx+avgMap[beginIdx]; iStep<endIdx; 
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
        barrier->pathUpdated(pathGen);

        if (pathGen->doingPast()){ 
            sumSoFar = sum;

            // In order to satisfy the BARRIER_LEVEL request we need
            // to know the ref level for each asset. This is only
            // needed for already started case. By conditioning via this
            // "if(doingPast)" we further restrict it so the BARRIER_LEVEL
            // will only be reported once a future sample is past. It is
            // far simpler to implement - otherwise getting access to the
            // ref level is non-trivial, and I'm not inclined to change
            // the design purely to satisfy this output request.
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }
            
            // for ibarrier case, get and store rebate value from historical hit
            if( hasFuture() )
            {
                payoffHist = rebateMgr->calculate(inst->valueDate, pathGen)*matFV;
                if( !!floatMgr )
                    payoffHist += floatMgr->calculate(inst->valueDate, pathGen)*matFV;
                if( !!fixedMgr )
                    payoffHist += fixedMgr->calculate(inst->valueDate, pathGen)*matFV;
            }

            // get known cfs
            PhysicalDeliveryArray phyDs; // dummy
            if( !!rebateMgr )
                rebateMgr->getKnowCFPhyD(inst->valueDate, pathGen, inst->notional, *knownCFs, phyDs);
            if( !!floatMgr ) {
                CashFlowArraySP fltCFs(new CashFlowArray(0));
                floatMgr->getKnowCFPhyD(inst->valueDate, pathGen, inst->notional, *fltCFs, phyDs);
                knownCFs = CashFlow::merge(knownCFs, fltCFs);
            }
            if( !!fixedMgr ) {
                CashFlowArraySP fixCFs(new CashFlowArray(0));
                fixedMgr->getKnowCFPhyD(inst->valueDate, pathGen, inst->notional, *fixCFs, phyDs);
                knownCFs = CashFlow::merge(knownCFs, fixCFs);
            }

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

            // 3. apply any modifiers (call/put/etc)
            basketComponents->apply();

            // 4. aggregate into a basket
            basket() = assetBasket->aggregate();

            // 5. apply any modifier (call/put/etc) to basket
            overallOption->apply();

            // 6. form "expected payoff" combining overall option and rebate
            // Note if condition NOT met then rebate level at mat (last step) is used.
            double koFactor = barrierStat->calculate(inst->valueDate, pathGen);
            double payoff = payoffHist + koFactor * basket();
            if( !!rebateMgr )
                payoff += rebateMgr->calculate(inst->valueDate, pathGen)*matFV;
            if( !!floatMgr )
                payoff += floatMgr->calculate(inst->valueDate, pathGen)*matFV;
            if( !!fixedMgr )
                payoff += fixedMgr->calculate(inst->valueDate, pathGen)*matFV;

            // Done...
            prices.add(inst->notional * payoff); 
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        barrier->adjustLN(pathGen,
                          this,
                          getMultiFactors(),
                          0);                          
    }

    /** Use this opportunity to do any Implied driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustments */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{
        barrier->adjustImplied(pathGen,
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
        static const string method = "RainbowDKO::getVolInterp";

        try
        {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime&    today = getToday();
            const DateTime&    startDate = getRefLevel()->getAllDates().front();
            bool               fwdStarting = startDate.isGreater(today);
            const DateTime&    lastSimDate = getSimSeries()->getLastDate();
            double interpLevel = inst->assetBasketComponents->getInterpLevel(iAsset) * 
                (fwdStarting? 1.0 : pathGen->refLevel(iAsset, 0));

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
        // BARRIER_LEVEL ...
        OutputRequest*  request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        const DateTime& today = getToday();
        const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
        // Only try to satisfy this request when have some past (so easy to get ref levels).
        if (request && (today >= lastRefDate)) {
            // This operates on a per-asset basis
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                // Mostly delegate to Barrier class...
                // histRefLevels are to allow absolute barriers to be reported.
                BarrierLevelArraySP levels = 
                    barrier->reportLevels(today, histRefLevels[iAsset], iAsset);
                if (!levels->empty()) {
                    OutputRequestUtil::recordBarrierLevels(
                        control, results,
                        getMultiFactors()->assetGetTrueName(iAsset),
                        levels.get());
                }
                
            }
        }            
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if ( request && !knownCFs->empty() ) {
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    knownCFs.get());
        }            

        // PAYMENT_DATES  
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            // in the known cash flow (i.e. libor + known Range Coupon)
            DateTimeArray cashPayDates = DateTimeArray(0);
            for (int i=0;i<(*knownCFs).size();i++){
                cashPayDates.push_back((*knownCFs)[i].date);
            }
            OutputRequestUtil::recordPaymentDates(control,results,&cashPayDates);
        }

    }

    // go get the events based on this path
    virtual void getEvents(const IMCPathGenerator*  pathGen, 
                           EventResults* events,
                           const DateTime& eventDate) {

        static const string method = "RainbowDKOMC::getEvents";
        try {
            int numAssets = getNumAssets();
            StringArraySP assetNames(new StringArray(numAssets));
            for (int i = 0; i < numAssets; i++) {
                (*assetNames)[i] = getMultiFactors()->assetGetTrueName(i);
            }
            
            barrier->pathUpdated(pathGen);
            
            barrierStat->getEvents(pathGen, eventDate, events, "Final Option Barriers", 
                                   assetNames.get());
            
            rebateMgr->getEvents(pathGen, eventDate, events, "Rebate Barrier", 
                                 assetNames.get());

            if( !!floatMgr)
                floatMgr->getEvents(pathGen, eventDate, events,"Floating Leg Barrier",
                                    assetNames.get());

            if (!!fixedMgr)
                fixedMgr->getEvents(pathGen, eventDate, events, "Fixed Leg Barrier",
                                    assetNames.get());
        } catch (exception& e) {
            throw ModelException(e,method);
        }

    }

};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RainbowDKO::createProduct(const MonteCarlo* model) const {
            
    // need to copy barriers to be passed to product. we need both barrierUtil here 
    // to get the sim dates and later used in product. do it here to avoid duplicate
    BarrierSP realBarrier = BarrierSP(copy(barrierUnion->getBarrier()));
    BarrierSP streamBarrier = realBarrier;
    if( useKiTrigger ) {
        realBarrier = BarrierSP(new DoubleBarrier(true/*kiKeepKO*/, false/*simultaneous*/, realBarrier, 
                                                  BarrierSP(copy(triggerUnion->getBarrier()))));
    }

    // Validate the barrier based on the model
    realBarrier->validate(model, assets->NbAssets());

    // make a copy of the barrier to give to product
    bool useBarrUtil = realBarrier->createBarrierUtil(valueDate);
    if( !useBarrUtil )
        throw ModelException("RainbowDKO::createProduct", "Can not create BarrierUtil. Check that no smoothing allowed");

    bool hasContinuous;
    DateTimeArray monitoringDates;
    realBarrier->getBarrierUtil()->relevantDates(monitoringDates, hasContinuous);

    // currently we do NOT add dates to monitoring dates if barrier is continuous monitoring 
    // but the barrier schedule does not contain enough monitoring dates
    if( monitoringDates.size() < 1 )
        throw ModelException("RainbowDKO::createProduct", "Barrier schedules must contain at least 1 monitor date");

    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty one */
    simSeries->addDates(averageOutDates);
    simSeries->addDates(monitoringDates);
    return new RainbowDKOMC(this, simSeries, realBarrier, streamBarrier);
}

// implementation of Events interface
void RainbowDKO::getEvents(const BarrierBreach* breach,
                          IModel* model, 
                          const DateTime& eventDate,
                          EventResults* events) const {
    static const string method = "RainbowDKO::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));
            MCPathGeneratorSP pastPathGenerator(
                    mc->getPathConfig()->pastPathGenerator(prod.get()));
            prod->getEvents(pastPathGenerator.get(), events, eventDate);
        } else {
            throw ModelException(method, 
                    "Internal error - expected Monte Carlo model for RDK pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const RainbowDKO::TYPE = CClass::registerClassLoadMethod(
    "RainbowDKO", typeid(RainbowDKO), RainbowDKO::load);

// * for class loading (avoid having header file) */
bool RainbowDKOLoad() {
    return (RainbowDKO::TYPE != 0);
}

DRLIB_END_NAMESPACE
