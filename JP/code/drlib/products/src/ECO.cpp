//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ECO.cpp
//
//   Description : Payoff = Notional * ( OverallPerf - Cash)
//                 OverallPerf = GenPerf(Sum(GenPerf(Asset[1,...,n]))) 
//
//   Date        : July 2003
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
class ECO: public GenericNFBase, 
           virtual public IMCIntoProduct{
protected:
    /// fields 
    IDoubleArrayModifierMakerSP  overallOption; // call/put/digital etc on sum of ...
    IDoubleArrayModifierMakerSP  assetPerfs;    // call/put/digital etc on each asset
    DateTimeArray                averageOutDates; // purely for getting average perf of each asset, XXX what about extremes?
    IRebateMakerSP               rebate;        // the "opposite" side of the payoff from 'assetPerfs'

    BarrierUnionSP               barrierUnion;  // determines KO "prob" per asset SOMEHOW!
    DateTimeArray                monitoringDates; // feel this should be part of barrier class XXX

public:
    static CClassConstSP const TYPE;
    friend class ECOMC;

    // validation
    void validatePop2Object(){
        static const string method = "ECO::validatePop2Object";
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
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
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

private:
    ECO(): GenericNFBase(TYPE) {} // for reflection
    ECO(const ECO& rhs);     // not implemented
    ECO& operator=(const ECO& rhs); // not implemented

    static IObject* defaultECO(){
        return new ECO();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ECO, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultECO);
        FIELD(overallOption, "Overall Perf Modifier");
        FIELD(assetPerfs, "Asset Perf Modifier");
        FIELD(averageOutDates, "Averaging out dates for basket");
        FIELD(rebate, "Opposite side of payoff to assetPerfs");
        FIELD_MAKE_OPTIONAL(rebate);
        FIELD(barrierUnion, "Barrier Union");
        FIELD(monitoringDates, "Monitoring Dates for barrier");
    }
};

/* MC product class for ECO */
class ECOMC: public IMCProduct,
             virtual public IMCProductLN,
             virtual public IMCProductImplied{
private:
    const ECO*                inst;
    BarrierSP                 barrier;   // non-const note
    int                       nbAssets;  // convenient

    DoubleArray               koFactors; // [nbAssets]
    SimpleDoubleArray         assetPerfs;
    IRebateSP                 rebate;
    IDoubleArrayModifierSP    basketComponents;
    TrivialDoubleArray        basket;  // just a double ... being the sum of per-asset "perfs"
    IDoubleArrayModifierSP    overallOption;

    IntArray                  avgMap;    // to track averaging dates
    DoubleArray               sum;       // [nbAssets], saves alloc later

    DoubleArray               basketParts; // [nbAssets]

    // for past
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 

    // To allow caching of product-level info across greeks
    class PricesECO: public MCPricesSimple{
        // cached params
        DoubleArrayArraySP cachedBasketParts; // [path idx][asset idx]

        // working params
        IntArray           changedAssets;
        bool               doingGreek;

    protected:
        IMCPrices* emptyConstructor() const{
            return new PricesECO(false, 1, 1, 1);
        }

    public:
        //// the assets whose path has changed
        const IntArray& getChangedAssets(){
            return changedAssets;
        }
        //// returns the cached data
        const DoubleArray& cacheRead(int pathIdx) {
            return (*cachedBasketParts)[pathIdx];
        }
        //// can we read from the cache?
        bool cacheValid(){
            return (doingGreek && cachedBasketParts.get());
        }
        //// can we write to the cache?
        bool cacheUpdateAllowed(){
            return (!doingGreek && cachedBasketParts.get());
        }
        //// save the supplied basketParts
        void cacheWrite(int                pathIdx,
                        const DoubleArray& basketParts){
            (*cachedBasketParts)[pathIdx] = basketParts;
        }
        //// number of bytes used per path in current configuration
        virtual int storagePerPath(IMCProduct* product) const {
            int n = MCPricesSimple::storagePerPath(product);
            if (cachedBasketParts.get()){
                n += sizeof(double) * product->getNumAssets();
            }
            return n;
        }

        /** Returns a deep copy of this object */
        IMCPrices* clone() const{
            PricesECO& copy = dynamic_cast<PricesECO&>(*MCPricesSimple::clone());
            copy.cachedBasketParts = cachedBasketParts; // shallow copy
            copy.doingGreek = doingGreek;
            return &copy;
        }
        PricesECO(bool              useCache,
                  int               numAssets,
                  int               nbIter,
                  int               nbSubSamples): 
            MCPricesSimple(nbIter, nbSubSamples), 
            changedAssets(numAssets), doingGreek(false){
            if (useCache){
                cachedBasketParts = DoubleArrayArraySP(new DoubleArrayArray(nbIter));
                for (int i = 0; i < nbIter; i++){
                    (*cachedBasketParts)[i].resize(numAssets);
                }
            }
            // fill changedAssets with 0, 1, 2, ...
            for (int i = 0; i < numAssets; i++){
                changedAssets[i] = i;
            }
        }
        /** invoked before each greek calculation (but not before initial
            pricing run) */
        virtual void configureCache(const IntArray& changedAssets){
            if (cachedBasketParts.get()){
                this->changedAssets = changedAssets;
            }
            doingGreek = true;
        }            

        virtual ~PricesECO() {}
        
    };

public:

    /**  */
    ECOMC(const ECO*            inst,
          const DateTimeArray&  monitorDates,
          const SimSeriesSP&    simSeries):
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
        koFactors(nbAssets, 0.0),
        assetPerfs(nbAssets, 0.0),
        basket(0.0),
        sum(nbAssets, 0.0),
        basketParts(nbAssets, 0.0),
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0) {

        // allow no rebate to be backwards compatible for test files
        if (inst->rebate.get()) {
            rebate = inst->rebate->getRebate(simSeries->getAllDates(),discount);
        } else {
            IRebateMakerSP rm = IRebateMakerSP(new FlatRebateMaker(0.0));
            rebate = rm->getRebate(simSeries->getAllDates(),discount);
        }

        basketComponents = IDoubleArrayModifierSP(inst->assetPerfs->getModifier(&assetPerfs));
        overallOption = IDoubleArrayModifierSP(inst->overallOption->getModifier(&basket));

        bool isTrivial;
        avgMap = DateTime::createMapping(simSeries->getAllDates(),
                                         inst->averageOutDates,
                                         isTrivial);

        // create an interpolated barrier to match the monitoring dates
        barrier->createInterpBarrier(inst->valueDate, monitorDates);
    }

    virtual ~ECOMC() {}

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        static const string method = "ECO::payoff";
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset, i;
        int    pathIdx = pathGen->getPathIndex();

        // 1. determine raw performances for each asset
        // 2. establish KO multipliers ("probs") for each asset
        // 3. apply any performance modifiers (call/put/etc)
        // 4. form "expected payoff per asset" & aggregate into "basket" - XXX what about "weights"?
        // 5. apply any modifier (call/put/etc) to basket
        // Not yet considered:
        // CashFlow

        // Initialise for potentially cached ones
        PricesECO& myPrices = static_cast<PricesECO&>(prices);
        if (myPrices.cacheValid()){
            basketParts = myPrices.cacheRead(pathIdx);
        }
        const IntArray& changedAssets = myPrices.getChangedAssets();

        for(i=0; i<changedAssets.size(); i++) {
            iAsset = changedAssets[i];

            // Initialise for non-cached values, 
            sum[iAsset] = sumSoFar[iAsset];

            // 1a. determine raw performances for each asset
            const double* path = pathGen->Path(iAsset, 0);
            for (int iStep=beginIdx+avgMap[beginIdx]; iStep<endIdx; 
                 iStep++, iStep += avgMap[iStep]) {
                sum[iAsset] += path[iStep];
            }

            // 2. establish KO multipliers ("probs") for each asset
            // Note this is outside the "if past" clause, which allows
            // the barrier to handle its own past.
            koFactors[iAsset] = barrier->hitValue(pathGen, iAsset);
        }

        if (pathGen->doingPast()){ 
            sumSoFar = sum;
            // This for BARRIER_LEVEL request. See RainbowKO.cpp for more complete comment
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }
        }
        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().

            double rebateValue = rebate->getLevel(endIdx-1);

            // 1b. finish this job
            for(i=0; i<changedAssets.size(); i++) {
                iAsset = changedAssets[i];
                double avg = sum[iAsset]/inst->averageOutDates.size();
                assetPerfs[iAsset] = avg/pathGen->refLevel(iAsset, 0);

                // 3. apply any modifiers (call/put/etc)
                // NB This CHANGES 'assetPerfs'
                // A little dodgy : passing index is potentially risky...
                basketComponents->apply(iAsset);

                // 3b. form "expected payoff per asset"
                basketParts[iAsset] = koFactors[iAsset] * assetPerfs[iAsset] +
                    (1.0 - koFactors[iAsset]) * rebateValue;
            }

            if (myPrices.cacheUpdateAllowed()){
                myPrices.cacheWrite(pathIdx, basketParts);
            }
            
            // 4. aggregate into "basket". No more caching from here
            double& bask = basket();
            bask = 0.0;
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                bask += basketParts[iAsset];
            }
            // KKK could move this into basketParts[] above
            bask /= nbAssets; // equal weights for now...

            // 5. apply any modifier (call/put/etc) to basket
            overallOption->apply();
            
            // Done...
            prices.add(inst->notional * basket()); 
        }
    }

    IMCPrices* createOrigPrices(int  nbIter,
                                        int  nbSubSamples,
                                        int  mode) {
        return new PricesECO((mode & CACHE_PRODUCT_BIT)? true: false, // use cache?
                             nbAssets,
                             nbIter, 
                             nbSubSamples);
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
        static const string method = "ECO::getVolInterp";

        try
        {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime&    today = getToday();
            const DateTime&    startDate = getRefLevel()->getAllDates().front();
            bool               fwdStarting = startDate.isGreater(today);
            const DateTime&    lastSimDate = getSimSeries()->getLastDate();
            double interpLevel = inst->assetPerfs->getInterpLevel(iAsset) * 
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
        the IMCProduct, but need to add BARRIER_LEVEL support. Should review. XXX 
        See also RainbowKO.cpp */
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

};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* ECO::createProduct(const MonteCarlo* model) const {

    // Validate the barrier based on the model
    Barrier* barrier = barrierUnion->getBarrier();
    if (!barrier) {
        throw ModelException("ECO::createProduct", "Failed to locate barrier instance!");
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

    return new ECOMC(this, monitorDates, simSeries);
}


CClassConstSP const ECO::TYPE = CClass::registerClassLoadMethod(
    "ECO", typeid(ECO), ECO::load);

// * for class loading (avoid having header file) */
bool ECOLoad() {
    return (ECO::TYPE != 0);
}

DRLIB_END_NAMESPACE









