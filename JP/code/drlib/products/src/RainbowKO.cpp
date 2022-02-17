//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RainbowKO.cpp
//
//   Description : Payoff = Notional * { OptionPayoff (if out & not hit, or in & hit) else Rebate }
//                 allow rebate pay at hit
//                 all breach decisions are NOT simultaneous
//                 OptionPayoff = GenOverallOption( Aggregate( GenPerfsPerAsset ) )
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

DRLIB_BEGIN_NAMESPACE
class RainbowKO: public GenericNFBase, 
                 virtual public IMCIntoProduct, 
                 virtual public BarrierBreach::IEventHandler {
protected:
    /// fields 
    IDoubleArrayModifierMakerSP  overallOption; // call/put/digital etc on sum of ...
    IAggregateMakerSP            assetBasket;   // pct/rainbow/product basket etc of ...
    IDoubleArrayModifierMakerSP  assetBasketComponents;    // call/put/digital etc on each asset
    DateTimeArray                averageOutDates; // purely for getting average perf of each asset, XXX what about extremes?
    string                       assetMeasurement;  // "Avg", "Min", "Max"
    IRebateMakerSP               rebate;        // the "opposite" side of the payoff from 'assetBasket'
    bool                         rebatePayAtHit;  // ?put this inside rebate?

    BarrierUnionSP               barrierUnion;  // determines KO "prob" -> hitValue()
    DateTimeArray                monitoringDates; // feel this should be part of barrier class XXX

    // allow for 1KI and 1KO
    bool                         useKiTrigger;  // true if use triggerUnion as KI and barrierUnion has to be KO
    BarrierUnionSP               triggerUnion; 

    // transient
    BarrierSP                    realBarrier;

public:
    static CClassConstSP const TYPE;
    friend class RainbowKOMC;

    // validation
    void validatePop2Object(){
        static const string method = "RainbowKO::validatePop2Object";
        GenericNFBase::validatePop2Object();

        try
        {
            (barrierUnion->getBarrier())->validate(assets->NbAssets());

            // painful but too much code implicitly depending on this in Barrier class
            // require all simulation dates are monitoring dates - so ...
            if (!DateTime::isSubset(monitoringDates, averageOutDates)) {
                throw ModelException(method,
                                     "Require averageOutDates to be a subset of monitoringDates");
            }

            
            // some processing for 1KI/1KO barrier 
            if( useKiTrigger )
            {
                // must be barrierSchedule and correct KO/KI field
                BarrierSchedule *ko = dynamic_cast<BarrierSchedule*>(barrierUnion->getBarrier());
                BarrierSchedule *ki = dynamic_cast<BarrierSchedule*>(triggerUnion->getBarrier());
                if( !ko || !ki || !ko->isKO() || ki->isKO() )
                    throw ModelException("RainbowKOMC", "If useKiTrigger, triggerUnion/barrierUnion must be KI/KO barrierSchedule respectively");
            
                realBarrier = BarrierSP(new DoubleBarrier(true/*kiKeepKO*/, false/*simultaneous*/, BarrierSP(ki), BarrierSP(ko)));
                realBarrier->validate(assets->NbAssets());
            }
            else
                realBarrier = BarrierSP(copy(barrierUnion->getBarrier()));

            // if stream rebate, does not allow payAtMaturity nor non-empty paydates
            if( StreamRebateMaker::TYPE->isInstance(rebate.get()) )
            {
                if( !rebatePayAtHit )
                    throw ModelException("RainbowKOMC", "Must be payAtHit if has stream rebate");
            }
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
        rebate->getMarket(model, market.get());
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
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
    RainbowKO(): GenericNFBase(TYPE), assetMeasurement("Avg"), rebate(), rebatePayAtHit(true), useKiTrigger(false) {} // for reflection
    RainbowKO(const RainbowKO& rhs);     // not implemented
    RainbowKO& operator=(const RainbowKO& rhs); // not implemented

    static IObject* defaultRainbowKO(){
        return new RainbowKO();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RainbowKO, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultRainbowKO);
        FIELD(overallOption, "Overall Perf Modifier");
        FIELD(assetBasket, "How to aggregate asset perfs");
        FIELD(assetBasketComponents, "Asset Perf Modifier");
        FIELD(averageOutDates, "Averaging out dates for basket");
        FIELD(assetMeasurement, "Avg/Min/Max");
        FIELD_MAKE_OPTIONAL(assetMeasurement);
        FIELD(rebate, "Opposite side of payoff to assetBasket");
        FIELD_MAKE_OPTIONAL(rebate);
        FIELD(rebatePayAtHit, "False=> at maturity");
        FIELD_MAKE_OPTIONAL(rebatePayAtHit);
        FIELD(barrierUnion, "Barrier Union");
        FIELD(monitoringDates, "Monitoring Dates for barrier");
        FIELD(useKiTrigger,  "true if use triggerUnion for KI and barrierUnion has to be KO");
        FIELD_MAKE_OPTIONAL(useKiTrigger);
        FIELD(triggerUnion,         "ignored if useTriggerKI is false");
        FIELD_MAKE_OPTIONAL(triggerUnion);
        FIELD(realBarrier,          "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(realBarrier);
    }
};

/* MC product class for RainbowKO */
class RainbowKOMC: public IMCProduct,
                   virtual public IMCProductLN,
                   virtual public IMCProductImplied{
private:
    const RainbowKO*          inst;
    BarrierSP                 barrier;   // non-const note
    int                       nbAssets;  // convenient

    SimpleDoubleArray         assetComps;
    IDoubleArrayModifierSP    basketComponents;
    IAggregateSP              assetBasket;
    TrivialDoubleArray        basket;  // just a double ... being the aggregation of per-asset "perfs"
    IRebateSP                 rebate;
    IDoubleArrayModifierSP    overallOption;

    int                       assetMeasure; // 0 - Avg, 1 - Min, 2 - Max
    IntArray                  avgMap;    // to track averaging dates
    DoubleArray               sum;       // [nbAssets], saves alloc later
    DoubleArray               fvFactors; // Rebate gives value at hit. Anything else we need to adjust for
                                         // Such as rebate at mat, or  
                                         // past rebates dropping out via fvFactor=0.0

    // for past
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 

    // IBarrier elements
    mutable bool                useBarrUtil;
    mutable IBarrierMgrSP       barrierStat;
    mutable IBarrierMgrSP       rebateMgr;
    double                      matFV;
    double                      payoffHist; // value from historical observation

public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to RainbowKO) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    RainbowKOMC(const RainbowKO*          inst,
                const DateTimeArray&      monitorDates,
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
        barrier(copy(inst->realBarrier.get())),
        nbAssets(getNumAssets()),
        assetComps(nbAssets, 0.0),
        basket(0.0),
        rebate(inst->rebate->getRebate(simSeries->getAllDates(), discount)),
        sum(nbAssets, 0.0),
        fvFactors(monitorDates.size(), 1.0), // default is pay at mat
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0),
        useBarrUtil(false),
        payoffHist(0.0) {

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

        // we do interp barrier even in case of BarrierUtil for backward compatibility
        // in case instrument's monitoring dates are used
        barrier->createInterpBarrier(inst->valueDate, monitorDates);

        if (inst->rebatePayAtHit) {
            // Account for inst settlement, but forbid fixed settlement date
            if (CashSettleDate::TYPE->isInstance(settlement.get())) {
                throw ModelException("RainbowKOMC","Rebate pay at hit does not allow instrument settlement of CashSettleDate!");
            }

            // if inst override rebatePayAtHit to be false, can not use BarrierPay
            createBarrierPay();
        }

        // create an interpolated barrier to match the monitoring dates
        if( !useBarrUtil )
        {
            const DateTime& today = getToday();
            for(int i=0; i<fvFactors.size(); i++) {
                DateTime payDate = inst->rebatePayAtHit?monitorDates[i]:monitorDates.back();
                payDate = settlement->settles(payDate,0);
                if (payDate <= today) {
                    // past so already paid
                    fvFactors[i] = 0.0;
                } else {
                    // infrastructure will pv from payment date, so override for this case
                    fvFactors[i] = discount->pv(today, payDate) / pvFromPaymentDate(); 
                }
            }
        }

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

    }
    
    void createBarrierPay() const
    {
        DateTimeArraySP simDates(new DateTimeArray(simSeries->getAllDates()));
        useBarrUtil = barrier->createBarrierUtil(inst->valueDate);
        if( useBarrUtil )
        {
            DateTimeArray obsDates(1, simDates->back());
            barrierStat = BarrierFactory::makeManager(barrier->getBarrierUtil(), BarrierFactory::makeTrivialPay(obsDates));
            barrierStat->preprocess(inst->valueDate, simDates, discount);
            
            BarrierSchedule *barSchd = dynamic_cast<BarrierSchedule*>(barrier.get());
            if( !barSchd )
                throw ModelException("RainbowKOMC:createBarrierPay", "barrier must be BarrierSchedule in order to create rebate barrier pay");
            rebateMgr = BarrierFactory::makeManager(barrier->getBarrierUtil(), 
                IRebateMaker::createBarrierPay(inst->rebate.get(), barSchd->isKO(), settlement.get()));
            rebateMgr->preprocess(inst->valueDate, simDates, discount);
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        static const string method = "RainbowKO::payoff";
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
        double koFactor;
        bool isConditionMet;
        int metAtStep;
        if( !useBarrUtil )
            barrier->hitValueAndTime(pathGen, koFactor, isConditionMet, metAtStep);
        else
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
            if( useBarrUtil && hasFuture() )
                payoffHist = rebateMgr->calculate(inst->valueDate, pathGen) * matFV;
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
            double rebateValue, payoff = 0;
            if( !useBarrUtil )
            {                
                rebateValue = isConditionMet?
                    (fvFactors[metAtStep] * rebate->getLevel(metAtStep)) :
                    rebate->getLevel(endIdx-1);
                payoff = koFactor * basket() + (1.0 - koFactor) * rebateValue;
            }
            else
            {
                rebateValue = rebateMgr->calculate(inst->valueDate, pathGen) * matFV;
                koFactor = barrierStat->calculate(inst->valueDate, pathGen);
                payoff = payoffHist + koFactor * basket() + rebateValue;
            }

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
        static const string method = "RainbowKO::getVolInterp";

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
    }

    // go get the events based on this path
    virtual void getEvents(const IMCPathGenerator*  pathGen, 
                           EventResults* events,
                           const DateTime& eventDate) {
        int numAssets = getNumAssets();
        StringArraySP assetNames(new StringArray(numAssets));
        for (int i = 0; i < numAssets; i++) {
            (*assetNames)[i] = getMultiFactors()->assetGetTrueName(i);
        }

		// override the barrier today if necessary
		// otherwise barrier today will be the legal one
		barrier->overrideEventBarrier(eventDate);

        barrier->getEvents(pathGen, events, false, "Multi-asset barrier", 
                            assetNames.get());
    }
};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RainbowKO::createProduct(const MonteCarlo* model) const {
            
    // Validate the barrier based on the model
    realBarrier->validate(model, assets->NbAssets());

    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);

    smartPtr<IMCPathConfig> pathConfig = model->getPathConfig();
    DateTimeArray monitorDates = realBarrier->getFutureMonitoringDates(getValueDate(),
                                                                       monitoringDates,
                                                                       pathConfig);
    simSeries->addDates(monitorDates);

    return new RainbowKOMC(this, monitorDates, simSeries);
}

// implementation of Events interface
void RainbowKO::getEvents(const BarrierBreach* breach,
                          IModel* model, 
                          const DateTime& eventDate,
                          EventResults* events) const {
    static const string method = "RainbowKO::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));
            MCPathGeneratorSP pastPathGenerator(
                    mc->getPathConfig()->pastPathGenerator(prod.get()));
            prod->getEvents(pastPathGenerator.get(), events, eventDate);
        } else {
            throw ModelException(method, 
                    "Internal error - expected Monte Carlo model for RKO pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const RainbowKO::TYPE = CClass::registerClassLoadMethod(
    "RainbowKO", typeid(RainbowKO), RainbowKO::load);

// * for class loading (avoid having header file) */
bool RainbowKOLoad() {
    return (RainbowKO::TYPE != 0);
}

DRLIB_END_NAMESPACE
