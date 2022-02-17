//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BoostedNFB.cpp
//
//   Description : Port of RBKBonus
//
//   Date        : June 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/ITimeAggregate.hpp"
#include "edginc/HandlePaymentEvents.hpp"

DRLIB_BEGIN_NAMESPACE

class BoostedNFB: public GenericNFBase, 
                    virtual public IMCIntoProduct{
protected:
    /// fields 
    // perhaps makes sense for gen time agg to contain dates itself?
    ITimeAggregateMakerSP        timeAggregate;

    // potential bonus payments
    DateTimeArray                bonusDates;

    // determines monitoring across all periods. Subject to review after Barrier redesign
    // This barrier simply provides a hitValue
    BarrierUnionSP               barrierUnion;  // determines KO "prob" -> hitValue()
    DateTimeArray                monitoringDates; // feel this should be part of barrier class XXX
    
public:
    static CClassConstSP const TYPE;
    friend class BoostedNFBMC;

    // validation
    void validatePop2Object(){
        static const string method = "BoostedNFB::validatePop2Object";
        GenericNFBase::validatePop2Object();
  
        try
        {
            (barrierUnion->getBarrier())->validate(assets->NbAssets());

        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return monitoringDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    BoostedNFB(): GenericNFBase(TYPE) {} // for reflection
    BoostedNFB(const BoostedNFB& rhs);     // not implemented
    BoostedNFB& operator=(const BoostedNFB& rhs); // not implemented

    static IObject* defaultBoostedNFB(){
        return new BoostedNFB();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BoostedNFB, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultBoostedNFB);
        FIELD(timeAggregate, "How to combine values per time period");
        FIELD(bonusDates, "subset of period dates");
        FIELD(barrierUnion, "Barrier Union");
        FIELD(monitoringDates, "Monitoring Dates for barrier");
    }
};

//////////////////////////////////////////////////////////////////
// Tailored version for RBKBonus port from EDG
// Performance of the multi-period case was poor for the general
// BoostedNFB so form a specialised version.
// I am sure there are clean ways to share routines (such as vol int
//////////////////////////////////////////////////////////////////
class BoostedNFBMC: public IMCProduct,
                    virtual public IMCProductLN,
                    virtual public IMCProductImplied,
                    virtual public IHandlePaymentEvents,
                    virtual public ITimeAggProductView {
private:
    const BoostedNFB*         inst;
    BarrierSP                 barrier;   // non-const note
    int                       nbAssets;  // convenient
    int                       nbPeriods;  // convenient

    SimpleDoubleArray         periodPerfs; // time dim across periods

    ITimeAggregateSP          timeAgg; // how to turn per-period gen perfs into a single value

    IntArray                  periodMap; // to track periods
    IntArray                  periodEndIdx; // if rebate not paid at hit then paid at end of each period
    DoubleArray               fvFactors; // if rebate at hit the this FVs to mat. 
                                         // Past rebates drop out via fvFactor=0.0
    // per-asset values (and pairs for past)
    DoubleArray               koFactor;
    DoubleArray               koFactorSoFar;
    DoubleArray               refLevels;
    BoolArray                 isConditionMet;
    BoolArray                 isConditionMetSoFar;

    // per-period values
    DoubleArray               periodKOFactor;
    BoolArray                 periodIsConditionMet;
    BoolArray                 isBonusPeriod;

    // for past
    SimpleDoubleArray         periodPerfsSoFar;
    int                       iPeriodSoFar;
    double                    bonusSoFar;
public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to BoostedNFB) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    BoostedNFBMC(const BoostedNFB*         inst,
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
        barrier(copy(inst->barrierUnion->getBarrier())),
        nbAssets(getNumAssets()),
        nbPeriods(inst->timeAggregate->getDates().size()),
        periodPerfs(nbPeriods, 0.0),
        fvFactors(monitorDates.size(), 1.0), // default is to leave to someone else
        koFactor(nbAssets, 0.0),
        koFactorSoFar(nbAssets, 0.0),
        refLevels(nbAssets,0.0),
        isConditionMet(nbAssets,false),
        isConditionMetSoFar(nbAssets,false),
        periodKOFactor(nbPeriods,0.0),
        periodIsConditionMet(nbPeriods,false),
        isBonusPeriod(nbPeriods, false),
        periodPerfsSoFar(nbPeriods,0.0),
        iPeriodSoFar(0),
        bonusSoFar(0.0) {
        static const string routine = "BoostedNFBMC::BoostedNFBMC";

        // Tie the pieces together :
        // We need the timeAgg to give a final value
        timeAgg = ITimeAggregateSP(inst->timeAggregate->getAggregate(simSeries->getAllDates().back(),
                                                                     this,
                                                                     &periodPerfs));

        bool isTrivial;
        const DateTimeArray& allSimDates = simSeries->getAllDates();
        const DateTimeArray& periodDates = inst->timeAggregate->getDates();
        periodMap = DateTime::createMapping(allSimDates,
                                            periodDates,
                                            isTrivial);
        // periodEndIdx[iPeriod] is the step index for the end of the appropriate period
        periodEndIdx = IntArray(periodDates.size());
        int iStep, iPeriod, iBonus=0;
        for(iStep=0, iStep += periodMap[iStep], iPeriod=0; 
            iStep<allSimDates.size(); 
            iStep++, iStep += periodMap[iStep], iPeriod++) {
            if (iPeriod>=periodDates.size()) {
                throw ModelException(routine,
                                     "Internal error trying to find period end indexes!");
            }
            periodEndIdx[iPeriod] = iStep;

            // subset of period dates are the bonus dates :
            if (iBonus<inst->bonusDates.size() &&
                inst->bonusDates[iBonus] == periodDates[iPeriod]) {
                isBonusPeriod[iPeriod] = true;
                iBonus++;
            }
        }
        if (iBonus<inst->bonusDates.size()) {
            throw ModelException(routine,
                                 "Bonus dates must be a subset of period dates but problem after #" +
                                 Format::toString(iBonus+1) + " = " +
                                 inst->bonusDates[iBonus].toString());
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
        int    iAsset;
        bool   doingPast = pathGen->doingPast();

        double bonus = bonusSoFar;
        int iPeriod = iPeriodSoFar;
        periodPerfs = periodPerfsSoFar;

        // need to cope with init from past here where have awareness of 
        // periods - not in barrier
        koFactor = koFactorSoFar;
        isConditionMet = isConditionMetSoFar; // array - [nbAssets]
        int metAtStep; // no actually needed
        
        // Form the asset basket at each Period date first
        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            // At each step, check for KO by still-active assets 
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                barrier->hitValueAndTimeAtStep(pathGen,
                                               iStep, 
                                               iAsset,
                                               pathGen->refLevel(iAsset, 0),
                                               koFactor[iAsset],  // will be unsmoothed value (i.e. 0/1)
                                               isConditionMet[iAsset], 
                                               metAtStep);
            }

            if (periodMap[iStep]==0) {  // true iff a Period date
                // Round things off for this time period

                // 1. ko info for this period
                // ----------------------------
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    barrier->hitValueAndTimeAtMat(pathGen,
                                                  iStep,
                                                  iAsset,
                                                  pathGen->refLevel(iAsset, 0),
                                                  iPeriod==iPeriodSoFar && isConditionMetSoFar[iAsset], // only the current period can have a past hit
                                                  koFactor[iAsset], 
                                                  isConditionMet[iAsset], 
                                                  metAtStep);
                }
                double periodBinaryValue;
                barrier->multiHelper2(koFactor, isConditionMet, 
                                      periodKOFactor[iPeriod],
                                      periodBinaryValue);  // unsmoothed value (i.e. 0 or 1) so indicates whether payment made or not

                // 2. Figure out any payment and bonus
                double toPay = 1.0;
                if (Maths::isPositive(periodBinaryValue)) {
                    if (isBonusPeriod[iPeriod]) {
                        // supplement with any pending bonus
                        toPay += bonus;
                        bonus = 0.0;
                    }
                } else {
                    // accumulate in bonus
                    bonus += 1.0;
                }
                // note this smoothes any bonus payment too
                periodPerfs[iPeriod] = periodKOFactor[iPeriod] * toPay;

                // 3. Next Period
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    // Barrier restarts each period
                    isConditionMet[iAsset] = false; 
                }
                iPeriod++;
            }
        }

        if (doingPast) {
            koFactorSoFar = koFactor;
            isConditionMetSoFar = isConditionMet;

            bonusSoFar = bonus;
            periodPerfsSoFar = periodPerfs;
            iPeriodSoFar = iPeriod;

            // In order to satisfy the BARRIER_LEVEL request we need
            // to know the ref level for each asset. This is only
            // needed for already started case. By conditioning via this
            // "if(doingPast)" we further restrict it so the BARRIER_LEVEL
            // will only be reported once a future sample is past. It is
            // far simpler to implement - otherwise getting access to the
            // ref level is non-trivial, and I'm not inclined to change
            // the design purely to satisfy this output request.
            // We can be slightly clever here and only report for still active assets
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                refLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }

        }
        if (!doingPast || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            // Perform the aggregation (across periods)
            double payoff = timeAgg->aggregate();
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
    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string routine = "BoostedNFBMC::getVolInterp";
        CVolRequestLNArray reqarr(1);

        // Interp at barrier ...
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();
        bool               fwdStarting = startDate.isGreater(today);
        double interpLevel = barrier->getInterpLevel(pathGen,
                                                     this,
                                                     iAsset);
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                 startDate,
                                                                 lastSimDate,
                                                                 fwdStarting));
        return reqarr;
    }

    // Satisfy IHandlePaymentEvents interface
    void recordEvents(Control* control,
                      Results* results) {
        static const string method("BoostedNFBMC::recordEvents");
        try {
            // PAYMENT_DATES is a list of all dates on which payments may occur
            // including past and potential future dates.
            // For BoostedNFB this means asking the TimeAggregate
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                OutputRequestUtil::recordPaymentDates(control,results,timeAgg->getPaymentDates());
            }
            
            // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
            // and be supplied for all past cash flows, and any future ones
            // that are determined.
            // For BoostedNFB this means asking the TimeAggregate
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
            const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
            // Only try to satisfy this request when have some past (so easy to get ref levels).
            if (request && (today >= lastRefDate)) {
                // This operates on a per-asset basis
                for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                    // Mostly delegate to Barrier class...
                    // histRefLevels are to allow absolute barriers to be reported.
                    BarrierLevelArraySP levels =
                        barrier->reportLevels(today, refLevels[iAsset],
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
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

};



/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* BoostedNFB::createProduct(const MonteCarlo* model) const {

    // Validate the barrier based on the model
    Barrier* barrier = barrierUnion->getBarrier();
    if (!barrier) {
        throw ModelException("BoostedNFB::createProduct", "Failed to locate barrier instance!");
    }
    barrier->validate(model, assets->NbAssets());

    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    smartPtr<IMCPathConfig> pathConfig = model->getPathConfig();
    DateTimeArray monitorDates = barrier->getFutureMonitoringDates(getValueDate(),
                                                                   monitoringDates,
                                                                   pathConfig);
    simSeries->addDates(monitorDates);

    return new BoostedNFBMC(this, monitorDates, simSeries);
}


CClassConstSP const BoostedNFB::TYPE = CClass::registerClassLoadMethod(
    "BoostedNFB", typeid(BoostedNFB), BoostedNFB::load);

// * for class loading (avoid having header file) */
bool BoostedNFBLoad() {
    return (BoostedNFB::TYPE != 0);
}

DRLIB_END_NAMESPACE









