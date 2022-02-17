//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CalendarDrop.cpp
//
//   Description : Payoff = Notional * SumOverPeriods{ OptionPayoff (if out & not hit, or in & hit) else Rebate }
//                 allow rebate pay at hit
//                 all breach decisions are NOT simultaneous
//                 OptionPayoff = GenOverallOption( Aggregate( GenPerfsPerAsset ) )
//                 Assets may be dropped along the way
//
//   Date        : Feb 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/ITimeAggregate.hpp"
#include "edginc/IRebate.hpp"
#include "edginc/HandlePaymentEvents.hpp"

DRLIB_BEGIN_NAMESPACE

class CalendarDrop: public GenericNFBase, 
                    virtual public IMCIntoProduct, 
                    virtual public BarrierBreach::IEventHandler {
protected:
    /// fields 
    IDoubleArrayModifierMakerSP  overallOption; // call/put/digital etc on periodPieces ...

    // perhaps makes sense for gen time agg to contain dates itself?
    ITimeAggregateMakerSP        timeAggregate;
    IDoubleArrayModifierMakerSP  timeAggregateComponents;    // call/put/digital etc on each period perf

    // determines monitoring across all periods. Subject to review after Barrier redesign
    // This barrier selects between rebate and asset basket (below)
    BarrierUnionSP               barrierUnion;  // determines KO "prob" -> hitValue()
    DateTimeArray                monitoringDates; // feel this should be part of barrier class XXX
    
    // Defined across all periods
    IRebateMakerSP               rebate;        // the "opposite" side of the payoff from 'assetBasket'
    bool                         rebatePayAtHit;  // ?put this inside rebate?

    // From raw asset data across all periods -> basket of perfs (one per period)
    // These mappings are per asset (so apply across all periods)
    IAggregateMakerSP            assetBasket;   // pct/rainbow/product basket etc of ...
    // Since this is where assets come into things - this is where DROP is defined...
    IAssetFilterMakerSP          assetFilter;    
    IDoubleArrayModifierMakerSP  assetBasketComponents;    // call/put/digital etc on each asset

    // Captures raw asset data - across all periods
    // May reset perf-origin at each period (isCliquetStyle)
    bool                         avgFromStart;
    // how about introducing a "Sampler" class which has sampleDates and a variety of ways to get "raw perfs"
    // such as average, max/min, rolling avg -> this is the new MC idea of "state vars"
    // What about having flex vars of types like gen perf etc?
    DateTimeArray                averageOutDates; // purely for getting average perf of each asset, XXX what about extremes?
    bool                         isCliquetStyle;
    bool                         isRefPrevAvgOut;

public:
    static CClassConstSP const TYPE;
    friend class CalendarDropMC;

    // validation
    void validatePop2Object(){
        static const string method = "CalendarDrop::validatePop2Object";
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

            // XXX Would like to offer this but the barrier adjustment is flawed - it can only use the ref level from start
            // XXX I don't see why we shouldn't allow this with European style monitoring though.
/* while testing ...
            if (isCliquetStyle) {
                throw ModelException(method,
                                     "isCliquetStyle still being investigated - not yet suppored.");
                                     } */

        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const;

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return DateTime::merge(averageOutDates, monitoringDates);
    }

private:
    CalendarDrop(): GenericNFBase(TYPE) {} // for reflection
    CalendarDrop(const CalendarDrop& rhs);     // not implemented
    CalendarDrop& operator=(const CalendarDrop& rhs); // not implemented

    static IObject* defaultCalendarDrop(){
        return new CalendarDrop();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CalendarDrop, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultCalendarDrop);
        FIELD(overallOption, "Overall Perf Modifier");
        FIELD(timeAggregate, "How to combine values per time period");
        FIELD(timeAggregateComponents, "Period Perf Modifier");
        FIELD(assetBasket, "How to aggregate asset perfs");
        FIELD(assetFilter, "i.e. How to Drop");
        FIELD(assetBasketComponents, "Asset Perf Modifier");
        FIELD(averageOutDates, "Averaging out dates for basket");
        FIELD(rebate, "Opposite side of payoff to assetBasket");
        FIELD(rebatePayAtHit, "False=> at maturity");
        FIELD(barrierUnion, "Barrier Union");
        FIELD(monitoringDates, "Monitoring Dates for barrier");
        FIELD(avgFromStart, "avgFromStart");
        FIELD(isCliquetStyle, "isCliquetStyle");
        FIELD(isRefPrevAvgOut, "isRefPrevAvgOut");
    }
};

/* MC product class for CalendarDrop */
class CalendarDropMC: public IMCProduct,
                      virtual public IMCProductLN,
                      virtual public IMCProductImplied,
                      virtual public IHandlePaymentEvents,
                      virtual public ITimeAggProductView {
private:
    const CalendarDrop*       inst;
    BarrierSP                 barrier;   // non-const note
    int                       nbAssets;  // convenient
    int                       nbPeriods;  // convenient

    SimpleDoubleArray         assetPerfs; // asset dim at each period
    SimpleDoubleArray         filterMeasures; // spot perfs used to determine asset dropping
    SimpleDoubleArray         periodPerfs; // time dim across periods
    TrivialDoubleArray        overallPerf; // final perf after aggregated across time periods

    IDoubleArrayModifierSP    assetBasketComponents; // per-asset gen perfs
    IAssetFilterSP            assetFilter; // any dropping
    IAggregateSP              assetBasket; // how to turn per-asset gen perfs into a single value at each period
    IDoubleArrayModifierSP    timeAggComponents; // per-period gen perfs
    ITimeAggregateSP          timeAgg; // how to turn per-period gen perfs into a single value
    IRebateSP                 rebate;
    IDoubleArrayModifierSP    overallOption;

    IntArray                  avgMap;    // to track averaging dates
    IntArray                  periodMap; // to track periods
    IntArray                  periodEndIdx; // if rebate not paid at hit then paid at end of each period
    DoubleArray               sum;       // [nbAssets], saves alloc later
    DoubleArray               fvFactors; // if rebate at hit the this FVs to mat. 
                                         // Past rebates drop out via fvFactor=0.0
    // per-asset values (and pairs for past)
    DoubleArray               koFactor;
    DoubleArray               koFactorSoFar;
    DoubleArray               refLevels;
    DoubleArray               refLevelsSoFar;
    BoolArray                 isConditionMet;
    BoolArray                 isConditionMetSoFar;
    IntArray                  metAtStep;
    IntArray                  metAtStepSoFar;

    // XXX not entirely sure about this ...
    IntArray                  activeAssetsAtPeriodStart;
    IntArray                  activeAssetsAtPeriodStartSoFar;

    // per-period values
    DoubleArray               periodKOFactor;
    BoolArray                 periodIsConditionMet;
    IntArray                  periodMetAtStep;
    IntArray                  nbAvgOutPerPeriod; // [nbPeriods] 

    // Working area for N->1 collapse when may have drops (ctg => contiguous)
    DoubleArray               koFactorCtg;
    BoolArray                 isConditionMetCtg;
    IntArray                  metAtStepCtg;
    
    // for past
    SimpleDoubleArray         periodPerfsSoFar;
    int                       iPeriodSoFar;
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 
    
    // to make picking up events easier
    int                       lastHistoricIndex;
    IntArrayArraySP           eventMetAtStep;
    DoubleMatrix              eventRefLevels;
    DoubleMatrix              eventAssetLevels;
public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to CalendarDrop) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    CalendarDropMC(const CalendarDrop*    inst,
                   const DateTimeArray&   monitorDates,
                   const SimSeriesSP&     simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst),
        barrier(copy(inst->barrierUnion->getBarrier())),
        nbAssets(getNumAssets()),
        nbPeriods(inst->timeAggregate->getDates().size()),
        assetPerfs(nbAssets, 0.0),
        filterMeasures(nbAssets, 0.0),
        periodPerfs(nbPeriods, 0.0),
        overallPerf(0.0),
        rebate(inst->rebate->getRebate(simSeries->getAllDates(), discount)),
        sum(nbAssets, 0.0),
        fvFactors(monitorDates.size(), 1.0), // default is to leave to someone else
        koFactor(nbAssets, 0.0),
        koFactorSoFar(nbAssets, 0.0),
        refLevels(nbAssets,0.0),
        refLevelsSoFar(nbAssets,0.0),
        isConditionMet(nbAssets,false),
        isConditionMetSoFar(nbAssets,false),
        metAtStep(nbAssets,0),
        metAtStepSoFar(nbAssets,0),
        activeAssetsAtPeriodStart(nbAssets,0),
        activeAssetsAtPeriodStartSoFar(nbAssets,0),
        periodKOFactor(nbPeriods,0.0),
        periodIsConditionMet(nbPeriods,false),
        periodMetAtStep(nbPeriods,0),
        nbAvgOutPerPeriod(nbPeriods, 0),
        koFactorCtg(nbAssets,0.0),
        isConditionMetCtg(nbAssets, false),
        metAtStepCtg(nbAssets, 0),
        periodPerfsSoFar(nbPeriods,0.0),
        iPeriodSoFar(0),
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0),
        lastHistoricIndex(-1),
        eventRefLevels(nbPeriods, nbAssets),
        eventAssetLevels(nbPeriods, nbAssets) {
        static const string routine = "CalendarDropMC::CalendarDropMC";

        // Tie the pieces together :
        overallOption = IDoubleArrayModifierSP(inst->overallOption->getModifier(&overallPerf));
        // We need the timeAgg to give a value at the "overall strike date" ... i.e. the final sim date
        timeAgg = ITimeAggregateSP(inst->timeAggregate->getAggregate(simSeries->getAllDates().back(),
                                                                     this,
                                                                     &periodPerfs));
        timeAggComponents = IDoubleArrayModifierSP(inst->timeAggregateComponents->getModifier(&periodPerfs));
        assetFilter = IAssetFilterSP(inst->assetFilter->getFilter(simSeries->getAllDates(),
                                                                  &filterMeasures));
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetPerfs, assetFilter));
        assetBasketComponents = IDoubleArrayModifierSP(inst->assetBasketComponents->getModifier(&assetPerfs));

        bool isTrivial;
        const DateTimeArray& allSimDates = simSeries->getAllDates();
        avgMap = DateTime::createMapping(allSimDates,
                                         inst->averageOutDates,
                                         isTrivial);
        const DateTimeArray& periodDates = inst->timeAggregate->getDates();
        periodMap = DateTime::createMapping(allSimDates,
                                            periodDates,
                                            isTrivial);
        // periodEndIdx[iPeriod] is the step index for the end of the appropriate period
        periodEndIdx = IntArray(periodDates.size());
        int iStep, iPeriod;
        for(iStep=0, iStep += periodMap[iStep], iPeriod=0; 
            iStep<allSimDates.size(); 
            iStep++, iStep += periodMap[iStep], iPeriod++) {
            if (iPeriod>=periodDates.size()) {
                throw ModelException(routine,
                                     "Internal error trying to find period end indexes!");
            }
            periodEndIdx[iPeriod] = iStep;
        }

        // create an interpolated barrier to match the monitoring dates
        barrier->createInterpBarrier(inst->valueDate, monitorDates);
        // For clarity we require no monitoring dates after final period date (latter actually
        // defines limit of valuation)
        if (periodDates.back() < monitorDates.back()) {
            throw ModelException(routine,
                                 "Cannot have monitoring date (" +
                                 monitorDates.back().toString() + 
                                 ") after final period date (" +
                                 periodDates.back().toString() + ")"); 
        }
            
        // Rebate can be paid either at hit, or at the end of each period. 
        if (inst->rebatePayAtHit) {
            // Account for inst settlement, but forbid fixed settlement date
            if (CashSettleDate::TYPE->isInstance(settlement.get())) {
                throw ModelException(routine,"Rebate pay at hit does not allow instrument settlement of CashSettleDate!");
            }
            
            const DateTime& today = getToday();
            for(int i=0; i<fvFactors.size(); i++) {
                DateTime payDateNow = settlement->settles(monitorDates[i],0);
                if (payDateNow <= today) {
                    // past so consider already paid
                    fvFactors[i] = 0.0;
                } else {
                    // with timeAgg have values at each period end leaving subsequent PVing to the
                    // time agg
                    int iStepPeriodEnd = i + periodMap[i];
                    DateTime payDateMat = settlement->settles(monitorDates[iStepPeriodEnd],0);
                    fvFactors[i] /= discount->pv(payDateNow, payDateMat); 
                }
            }
        }

        iPeriod = 0; 
        for(int iAvg = 0; iAvg < inst->averageOutDates.size(); iAvg++) {    
            // number of AvgOut per period
            if(inst->averageOutDates[iAvg] > periodDates[iPeriod]) {
                iPeriod++;
            }
            nbAvgOutPerPeriod[iPeriod]++;
        }
        // cumulative number of AvgOut Dates
        if(inst->avgFromStart) {
            for(iPeriod = 1; iPeriod < periodDates.size(); iPeriod++) {
                nbAvgOutPerPeriod[iPeriod] += nbAvgOutPerPeriod[iPeriod-1];
            }
        } 
        // Check at least one average out per period 
        for(iPeriod = 0; iPeriod < periodDates.size(); iPeriod++) {
            if (nbAvgOutPerPeriod[iPeriod]<1) {
                throw ModelException(routine,
                                     "Require at least one average out date per period, check period #" +
                                     Format::toString(iPeriod+1));
            }
        }
        
        // events storage
        eventMetAtStep = IntArrayArraySP(new IntArrayArray(nbPeriods));
        
        for (int k = 0; k < nbPeriods; ++k){
            (*eventMetAtStep)[k] = IntArraySP(new IntArray(nbAssets, -1));
        }
    }

    // for ITimeAggProductView
    const DateTime& getValueDate() const {
        return getToday();
    }
    const YieldCurve* getYieldCurve() const {
        return discount;
    }
    DateTime settles(const DateTime& aDate) const {
        return settlement->settles(aDate, 0);
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset, iActive;
        bool   doingPast = pathGen->doingPast();

        int iPeriod = iPeriodSoFar;
        sum = sumSoFar;
        refLevels = refLevelsSoFar;
        periodPerfs = periodPerfsSoFar;

        // need to cope with init from past here where have awareness of 
        // periods - not in barrier
        koFactor = koFactorSoFar;
        isConditionMet = isConditionMetSoFar; // array - [nbAssets]
        metAtStep = metAtStepSoFar;// array - [nbAssets]
        
        // XXX how does the past work for this filtering??
        const IntArray& activeAssets = assetFilter->getActiveAssets(doingPast);
        if (iPeriod==0) {
            activeAssetsAtPeriodStart = activeAssets;
        } else {
            activeAssetsAtPeriodStart = activeAssetsAtPeriodStartSoFar;
        }

        for(iActive=0; iActive<activeAssets.size(); iActive++) {
            iAsset = activeAssets[iActive];
            /* Special treatment for average in for first cliq */
            if (!inst->isCliquetStyle || iPeriod == 0) {
                refLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            } else {
                refLevels[iAsset] = refLevelsSoFar[iAsset];
            }
        }

        // Form the asset basket at each Period date first
        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            // At each step, check for KO by still-active assets 
            for(iActive=0; iActive<activeAssets.size(); iActive++) {
                iAsset = activeAssets[iActive];

                barrier->hitValueAndTimeAtStep(pathGen,
                                               iStep, 
                                               iAsset,
                                               refLevels[iAsset],
                                               koFactor[iAsset],  // will be unsmoothed value (i.e. 0/1)
                                               isConditionMet[iAsset], 
                                               metAtStep[iAsset]);
                // we filter (i.e. drop) based on spot perf
                filterMeasures[iAsset] = pathGen->Path(iAsset, 0)[iStep] / refLevels[iAsset];

                // store some stuff for events reporting 
                // which will be hard to retrieve later
                if (doingPast && (isConditionMet[iAsset] == true)
                        && metAtStep[iAsset] == iStep) {
                    (*(*eventMetAtStep)[iPeriod])[iAsset] = metAtStep[iAsset];
                    eventAssetLevels[iPeriod][iAsset] = pathGen->Path(iAsset, 0)[iStep];
                    eventRefLevels[iPeriod][iAsset] = refLevels[iAsset];
                }
            }
            // Determine any dropped assets here. This happens after KO but before any baskets are made.
            // XXX hmmm - strictly need this only on drop dates ... should assetFilter allow such a flag?
            // NB relies on filterMeasures above
            assetFilter->update(iStep); 

            // Sum does not count any assets dropped this date
            if (avgMap[iStep]==0) { // true iff an averaging date
                for(iActive=0; iActive<activeAssets.size(); iActive++) {
                    iAsset = activeAssets[iActive];
                    sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];
                }
            }

            if (periodMap[iStep]==0) {  // true iff a Period date
                // Round things off for this time period

                // 1. ko info for this period
                // ----------------------------
                // The koFactor[] etc arrays may have "gaps".
                // Form contiguous arrays here where I know about "activeAssets" and
                // pass those to barrier.
                // Idea (to be confirmed) is to use all assets active when a period started
                // to determine hit condition. Assets that are dropped are not updated, so
                // will only contribute if they satisfied condition. Assets active at the end
                // of the period can contribute a smoothed hit value, which we get now :
                for(iActive=0; iActive<activeAssets.size(); iActive++) {
                    iAsset = activeAssets[iActive];
                    barrier->hitValueAndTimeAtMat(pathGen,
                                                  iStep,
                                                  iAsset,
                                                  refLevels[iAsset],
                                                  // only the current period can have a past hit
                                                  (doingPast && isConditionMet[iAsset]) ||
                                                  (iPeriod==iPeriodSoFar && isConditionMetSoFar[iAsset]), 
                                                  koFactor[iAsset], 
                                                  isConditionMet[iAsset], 
                                                  metAtStep[iAsset]);
                }
                // Now merge into single values
                koFactorCtg.resize(activeAssetsAtPeriodStart.size());
                isConditionMetCtg.resize(activeAssetsAtPeriodStart.size());
                metAtStepCtg.resize(activeAssetsAtPeriodStart.size());
                for(iActive=0; iActive<activeAssetsAtPeriodStart.size(); iActive++) {
                    iAsset = activeAssetsAtPeriodStart[iActive];
                    koFactorCtg[iActive] = koFactor[iAsset];
                    isConditionMetCtg[iActive] = isConditionMet[iAsset];
                    metAtStepCtg[iActive] = metAtStep[iAsset];
                }
                barrier->multiHelper(koFactorCtg, isConditionMetCtg, metAtStepCtg,
                                     periodKOFactor[iPeriod],
                                     periodIsConditionMet[iPeriod],
                                     periodMetAtStep[iPeriod]);

                // 2. deal with any asset basket
                // -----------------------------
                for(iActive=0; iActive<activeAssets.size(); iActive++) {
                    iAsset = activeAssets[iActive];
                
                    double avg = sum[iAsset] / nbAvgOutPerPeriod[iPeriod];
                    assetPerfs[iAsset] = avg / refLevels[iAsset];
                    assetBasketComponents->apply(iAsset); // this modifies assetPerfs[iAsset]

                    // this is a convenient place at which to reset some values for next period
                    // (we have 'avg' already)
                    if (inst->isCliquetStyle)
                    {
                        /* Reset ref level for next Period - done per asset note*/
                        if (inst->isRefPrevAvgOut) {
                            refLevels[iAsset] = avg;
                        } else {
                            refLevels[iAsset] = pathGen->Path(iAsset, 0)[iStep];
                        }
                    }
                    if (!inst->avgFromStart) {
                        sum[iAsset] = 0.0;
                    }

                    // Barrier restarts each period
                    isConditionMet[iAsset] = false; 
                }
                // Note that the assetBasket knows about 'activeAssets'
                periodPerfs[iPeriod] = assetBasket->aggregate();

                // 3. Last thing to do this step, when it's a period step, is to 
                // finish getting ready for the next Period
                activeAssetsAtPeriodStart = activeAssets;
                iPeriod++;
            }
        }

        if (doingPast) {
            koFactorSoFar = koFactor;
            isConditionMetSoFar = isConditionMet;
            metAtStepSoFar = metAtStep;
            lastHistoricIndex = endIdx - 1;

            sumSoFar = sum;
            refLevelsSoFar = refLevels;
            periodPerfsSoFar = periodPerfs;
            iPeriodSoFar = iPeriod;

            // XXX allow filter to preserve some state ... TEST ME!
            activeAssetsAtPeriodStartSoFar = activeAssetsAtPeriodStart;
            assetFilter->updateForPast();

            // In order to satisfy the BARRIER_LEVEL request we need
            // to know the ref level for each asset. This is only
            // needed for already started case. By conditioning via this
            // "if(doingPast)" we further restrict it so the BARRIER_LEVEL
            // will only be reported once a future sample is past. It is
            // far simpler to implement - otherwise getting access to the
            // ref level is non-trivial, and I'm not inclined to change
            // the design purely to satisfy this output request.
            // We can be slightly clever here and only report for still active assets
            for(iActive=0; iActive<activeAssets.size(); iActive++) {
                iAsset = activeAssets[iActive];
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }

        }
        if (!doingPast || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            // Now in the time dimension - apply perf modifiers to the perf for each period
            timeAggComponents->apply();

            // Then we have the KO feature & rebate of each period's RKO
            for (iPeriod=0; iPeriod<nbPeriods; iPeriod++) {
                // Note if condition NOT met then rebate level at period end should be used.
                double rebateValue = periodIsConditionMet[iPeriod]? 
                    (fvFactors[periodMetAtStep[iPeriod]] * rebate->getLevel(periodMetAtStep[iPeriod])) : 
                    rebate->getLevel(periodEndIdx[iPeriod]);
                // Note acting like a modifier : periodPerfs[iPeriod] is changed "in place"
                periodPerfs[iPeriod] = periodKOFactor[iPeriod] * periodPerfs[iPeriod] + 
                    (1.0 - periodKOFactor[iPeriod]) * rebateValue;
            }

            // then perform the aggregation (across periods)
            overallPerf() = timeAgg->aggregate();
            overallOption->apply();
            prices.add(inst->notional * overallPerf()); 
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
    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string routine = "CalendarDropMC::getVolInterp";
        CVolRequestLNArray reqarr(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel = inst->assetBasketComponents->getInterpLevel(iAsset);
        const DateTimeArray& periodDates = inst->timeAggregate->getDates();

        if (inst->isCliquetStyle) {
            // get hold of the future strike dates
            int numLiveCliqs = periodDates.size() - iPeriodSoFar;
            if (numLiveCliqs<=0) {
                throw ModelException(routine, "No future periods!?");
            }
            DateTimeArray liveCliqStartDates(numLiveCliqs);
            for (int iCliquet = 0; iCliquet < numLiveCliqs; iCliquet++){
                int iPeriod = iPeriodSoFar+iCliquet-1;
                liveCliqStartDates[iCliquet] = iPeriod<0?startDate:periodDates[iPeriod];
            }

            // same strike levels per cliquet (but may need to adjust first one)
            DoubleArray  strikes(numLiveCliqs, interpLevel);
            if (!fwdStarting){
                // need to set first level to absolute strike - adjusted
                // additionally for any average out samples for this cliquet
                int iStep;
                const DateTime& thisCouponStart = iPeriodSoFar==0?
                    startDate:periodDates[iPeriodSoFar-1];
                // find first avg date of this coupon
                for(iStep = 0; iStep < inst->averageOutDates.size() && 
                        inst->averageOutDates[iStep] <= thisCouponStart; iStep++) {
                    ; // empty
                }
                // then walk through counting past avg dates in this cliq
                int numRemaining = nbAvgOutPerPeriod[iPeriodSoFar];
                for(; iStep < inst->averageOutDates.size() && 
                        inst->averageOutDates[iStep] <= today; iStep++) {
                    numRemaining--;
                }
                // Can't set up refLevel earlier, 'cos need PathGen. First cliq has standard ref level
                double refLevel =  iPeriodSoFar==0 ? pathGen->refLevel(iAsset, 0) : refLevelsSoFar[iAsset];
                if (numRemaining<=0) {
                    // All averaging is done, but may still have monitoring so pick something reasonable
                    strikes[0] = interpLevel * refLevel;
                } else {
                    strikes[0] = (nbAvgOutPerPeriod[iPeriodSoFar] * refLevel * interpLevel
                                  - sumSoFar[iAsset])/ numRemaining;
                }
            }
            reqarr[0] =  CVolRequestLNSP(new CliquetVolRequest(fwdStarting, 
                                                               liveCliqStartDates, 
                                                               lastSimDate,
                                                               strikes));
        } else {
            // per asset 
            if (!fwdStarting){
                double reflvl = pathGen->refLevel(iAsset, 0);
                // some samples have fixed already (this includes averaging in)
                if (inst->avgFromStart) {
                    int numDates = inst->averageOutDates.size();
                    int numRemaining = 
                        today.numFutureDates(inst->averageOutDates);
                    
                    if (numRemaining<=0) {
                        // All averaging is done, but may still have monitoring so pick something reasonable
                        interpLevel = interpLevel * reflvl;
                    } else {
                        interpLevel = (numDates * interpLevel * reflvl
                                   - sumSoFar[iAsset])/ numRemaining;
                    }
                } else {
                    // Note we do not account for any past averaging dates within periods
                    interpLevel *= reflvl;
                }
            }
            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
        }
        return reqarr;
    }

	void overrideEventBarrier(const DateTime& eventDate) {
		barrier->overrideEventBarrier(eventDate);
	}

    // Satisfy IHandlePaymentEvents interface
    // for SPI specific record of events (pay dates, cashflows etc)
    void recordEvents(Control* control,
                      Results* results) {
        static const string method("CalendarDropMC::recordEvents");
        try {
            // PAYMENT_DATES is a list of all dates on which payments may occur
            // including past and potential future dates.
            // For CalendarDrop this means asking the TimeAggregate
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                OutputRequestUtil::recordPaymentDates(control,results,timeAgg->getPaymentDates());
            }
            
            // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
            // and be supplied for all past cash flows, and any future ones
            // that are determined.
            // For CalendarDrop this means asking the TimeAggregate
            request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished()) {
                const CashFlowArray* kfl = timeAgg->getKnownFlows();
                if (kfl && kfl->size()>0) {
                    OutputRequestUtil::recordKnownCashflows(control,
                                                            results,
                                                            discount->getCcy(),
                                                            kfl); 
                }
            }

            // BARRIER_LEVEL ...
            request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
            const DateTime& today = getToday();
            const IRefLevel* refLevel = getRefLevel();
            const DateTime& lastRefDate = refLevel->getAllDates().back();
            // Only try to satisfy this request when have some past (so easy to get ref levels).
            if (request && (today >= lastRefDate)) {
                // cliquet will require some change to this.
                if (inst->isCliquetStyle) {
                    string stuffed("Cliquet style does not support BARRIER_LEVEL request");
                    results->storeRequestResult(request,
                                                IObjectSP(new Untweakable(stuffed)),
                                                OutputNameConstSP(new OutputName("")));
                }
                else {
                    // This operates on a per-asset basis
                    for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                        // Mostly delegate to Barrier class...
                        // histRefLevels are to allow absolute barriers to be reported.
                        BarrierLevelArraySP levels =
                            barrier->reportLevels(today, histRefLevels[iAsset],
                                                  iAsset);
                        if (!levels->empty()) {
                            OutputRequestUtil::recordBarrierLevels(
                                control, results,
                                getMultiFactors()->assetGetTrueName(iAsset),
                                levels.get());
                        }
                    
                    }
                }
            }            
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // collects events which should be sitting there
    // called after past is run
    virtual void retrieveEvents(EventResults* events) const {
        StringArraySP assetNames(new StringArray(nbAssets));
        for (int i = 0; i < nbAssets; i++) {
            (*assetNames)[i] = 
                    getMultiFactors()->assetGetTrueName(i);
        }
        // if the instrument has expired then iPeriodSoFar may be too big
        barrier->getEvents(events, "Multi-asset barrier", 
                            Maths::min(iPeriodSoFar, nbPeriods-1), *eventMetAtStep,
                           *assetNames, eventAssetLevels, eventRefLevels);
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* CalendarDrop::createProduct(const MonteCarlo* model) const {

    // Validate the barrier based on the model
    Barrier* barrier = barrierUnion->getBarrier();
    if (!barrier) {
        throw ModelException("CalendarDrop::createProduct", "Failed to locate barrier instance!");
    }
    barrier->validate(model, assets->NbAssets());

    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);

    smartPtr<IMCPathConfig> pathConfig = model->getPathConfig();
    DateTimeArray monitorDates = barrier->getFutureMonitoringDates(getValueDate(),
                                                                   monitoringDates,
                                                                   pathConfig);
    simSeries->addDates(monitorDates);

    return new CalendarDropMC(this, monitorDates, simSeries);
}

// implementation of Events interface
void CalendarDrop::getEvents(const BarrierBreach* breach, IModel* model, 
                             const DateTime& eventDate,
                             EventResults* events) const {
    static const string method = "CalendarDrop::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));

			CalendarDropMC* actualProd = 
						dynamic_cast<CalendarDropMC*>(prod.get());

			// override the barrier today if necessary
			// otherwise barrier today will be the legal one
			actualProd->overrideEventBarrier(eventDate);

            // best just to run past what with the dropping
            MCPathGeneratorSP past = prod->runPast(mc);
            prod->retrieveEvents(events);
        } else {
            throw ModelException(method, 
                    "Internal error - expected Monte Carlo model for CKD pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const CalendarDrop::TYPE = CClass::registerClassLoadMethod(
    "CalendarDrop", typeid(CalendarDrop), CalendarDrop::load);

// * for class loading (avoid having header file) */
bool CalendarDropLoad() {
    return (CalendarDrop::TYPE != 0);
}

DRLIB_END_NAMESPACE









