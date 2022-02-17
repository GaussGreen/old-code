//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TriggerECO.cpp
//
//   Description : Payoff = Notional * ( OverallPerf - Cash)
//                 OverallPerf = GenPerf(Sum(GenPerf(Asset[1,...,n]))) 
//                 Trigger Condition on Basket Level.
//                 Based on ECO with Revision 1.14
//
//   Date        : August 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/IRebate.hpp"
#include "edginc/KOFixedLeg.hpp"
#include "edginc/KOLiborLeg.hpp"

DRLIB_BEGIN_NAMESPACE
class TriggerECO: public GenericNFBase, 
           virtual public IMCIntoProduct{
protected:
    /// fields 
    IDoubleArrayModifierMakerSP  overallOption; // call/put/digital etc on sum of ...
    IDoubleArrayModifierMakerSP  assetPerfs;    // call/put/digital etc on each asset
    DateTimeArray                averageOutDates; // purely for getting average perf of each asset, XXX what about extremes?
    IRebateMakerSP               rebate;        // the "opposite" side of the payoff from 'assetPerfs'
    bool                         rebatePayAtHit;  // ?put this inside rebate?

    BarrierUnionSP               barrierUnion;  // determines KO "prob" per asset SOMEHOW!
    DateTimeArray                monitoringDates; // feel this should be part of barrier class XXX

    BarrierUnionSP               triggerUnion;  // determines Trigger "prob" per asset SOMEHOW!
    IAggregateMakerSP            triggerBasket; // generic performance for Trigger Condition.

    bool                         hasFloater;
    bool                         hasFixedLeg;
    bool                         withPlainSwap; //True = add plain swap value, so as to return KO Swap.  
    LiborLegSP                   floater;
    FixedLegSP                   fixedLeg;
    string                       koStubRule;
    double                       redemption;    

public:
    static CClassConstSP const TYPE;
    friend class TriggerECOMC;

    // validation
    void validatePop2Object(){
        static const string method = "TriggerECO::validatePop2Object";
        GenericNFBase::validatePop2Object();

        try
        {
            (barrierUnion->getBarrier())->validate(assets->NbAssets());

            (triggerUnion->getBarrier())->validate(assets->NbAssets());

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

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
        // parent
        GenericNFBase::GetMarket(model, market);
        // relevant elements of self
        rebate->getMarket(model, market.get());
        if (hasFloater)
            floater->getMarket(model, market.get());
    }

    //  For SwapLeg, to feed fixing level when ValueDate <= fixing <= Shifted ValueDate 
    //  This would be moved to SwapLegIntFace, soon.
    bool sensShift(Theta* shift) {
        static const string method = "TriggetECO::sensShift";
        try  {
            const DateTime& newDate = shift->rollDate(valueDate);
            // and fixings
            if (hasFloater) {
                floater->setFixingforThetaShift(valueDate,
                    discount.get(),
                    newDate);
            }
        
            // roll the parent (updates value date etc)
            GenericNFBase::sensShift(shift);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
        return true; // our components have theta type sensitivity
    }

private:
    TriggerECO(): GenericNFBase(TYPE) {} // for reflection
    TriggerECO(const TriggerECO& rhs);     // not implemented
    TriggerECO& operator=(const TriggerECO& rhs); // not implemented

    static IObject* defaultTriggerECO(){
        return new TriggerECO();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TriggerECO, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultTriggerECO);
        FIELD(overallOption, "Overall Perf Modifier");
        FIELD(assetPerfs, "Asset Perf Modifier");
        FIELD(averageOutDates, "Averaging out dates for basket");
        FIELD(rebate, "Opposite side of payoff to assetPerfs");
        FIELD(barrierUnion, "Barrier Union");
        FIELD(monitoringDates, "Monitoring Dates for barrier");
        FIELD(triggerUnion, "Trigger, Barrier Union For Trigger Condition")
        FIELD(triggerBasket, "triggerBasket.  The performance to be used to judge Trigger or not.");
        FIELD(hasFloater, "has Libor Leg?")
        FIELD(hasFixedLeg, "has Fixed Leg?")
        FIELD(koStubRule, "KO Stub rule for swap leg stream");
        FIELD(withPlainSwap, "True = CashFlow is KO Swap.  False (Default)= CashFlow is -KI Swap");
        FIELD(fixedLeg, "Fixed Coupon Leg");
        FIELD_MAKE_OPTIONAL(fixedLeg);
        FIELD(floater, "libor leg");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(redemption, "redemption, paid at maturity.");
        FIELD(rebatePayAtHit, "False=> at maturity");
    }
};

/* MC product class for TriggerECO */
class TriggerECOMC: public IMCProduct,
             virtual public IMCProductLN,
             virtual public IMCProductImplied{
private:
    const TriggerECO*                inst;
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

    BarrierSP                 trigger;   
    IAggregateSP              triggerBasket;
    double                    redemption;
    DoubleArray               cashFlowAfterKO;  // the value after KO at each step.
    DoubleArray               fvFactors; // if rebate at hit the this FVs to mat. 
                                         // Past rebates drop out via fvFactor=0.0
    KOFixedLegSP              koFix;
    KOLiborLegSP              koFlt;
    KOStubRuleSP              koRule;    //need to declare here, otherwise it's only "alive" in specific scope.

    // for past
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 

    // flag to know already triggered or not.
    bool                      isTriggeredInPast;
    int                       iStepTrgInPast;

    // To allow caching of product-level info across greeks
    class PricesTriggerECO: public MCPricesSimple{
        // cached params
        DoubleArrayArraySP cachedBasketParts; // [path idx][asset idx]

        // working params
        IntArray           changedAssets;
        bool               doingGreek;

    protected:
        IMCPrices* emptyConstructor() const{
            return new PricesTriggerECO(false, 1, 1, 1);
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
            PricesTriggerECO& copy = dynamic_cast<PricesTriggerECO&>(*MCPricesSimple::clone());
            copy.cachedBasketParts = cachedBasketParts; // shallow copy
            copy.doingGreek = doingGreek;
            return &copy;
        }
        PricesTriggerECO(bool              useCache,
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

        virtual ~PricesTriggerECO() {}
        
    };

public:

    /**  */
    TriggerECOMC(const TriggerECO*            inst,
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
        trigger(copy(inst->triggerUnion->getBarrier())),
        fvFactors(inst->monitoringDates.size(), 1.0), // default is pay at mat
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0),
        isTriggeredInPast(false),
        iStepTrgInPast(0) {        

        rebate = inst->rebate->getRebate(simSeries->getAllDates(),discount);
        basketComponents = IDoubleArrayModifierSP(inst->assetPerfs->getModifier(&assetPerfs));
        overallOption = IDoubleArrayModifierSP(inst->overallOption->getModifier(&basket));

        bool isTrivial;
        avgMap = DateTime::createMapping(simSeries->getAllDates(),
                                         inst->averageOutDates,
                                         isTrivial);

        // create an interpolated barrier to match the monitoring dates
        barrier->createInterpBarrier(inst->valueDate, inst->monitoringDates);


        // create an interpolated trigger to match the monitoring dates
        // To Do : We'd like to allow different Trigger Schedule (but still subset) from monitoring.
        trigger->createInterpBarrier(inst->valueDate, inst->monitoringDates);
        
        // make the trigger Basket
        triggerBasket = IAggregateSP(inst->triggerBasket->getAggregate(&assetPerfs));

        // set up KO legs
        InstrumentSettlementSP mySettle(inst->instSettle.get());
        string payTimingRule = "AsSchedule";
        if (inst->koStubRule == "S")
            payTimingRule = "AsInstSettle";
        koRule = KOStubRuleSP(new KOStubRule(inst->koStubRule,
                                            false,              //isAccrueUpToSettle,
                                            payTimingRule));      //payTimingRule
        if (inst->hasFixedLeg){
            koFix = KOFixedLegSP(new KOFixedLeg(inst->fixedLeg,koRule));
        }
        if (inst->hasFloater){
            koFlt = KOLiborLegSP(new KOLiborLeg(inst->floater, koRule, mySettle, inst->valueDate, inst->discount));
        }

    }

    virtual ~TriggerECOMC() {}

    // Prepare the value after KO on each monitoring date, at here. 
    // I guess it's less comupte time doing at payoff.
    void finaliseSimulationDates(){
        static const string method = "TriggerECO::finaliseSimulationDates";
        int i;
        // rebate
        if (inst->rebatePayAtHit) {
            // Account for inst settlement, but forbid fixed settlement date
            if (CashSettleDate::TYPE->isInstance(settlement.get())) {
                throw ModelException(method,"Rebate pay at hit does not allow instrument settlement of CashSettleDate!");
            }
        }
        redemption = inst->redemption;
        
        // prepare pvFactors and cash flow
        try
        {
            const DateTime& today = getToday();
            cashFlowAfterKO = DoubleArray(inst->monitoringDates.size());
            for (i = 0; i < inst->monitoringDates.size(); i++){

                DateTime payDate = settlement->settles(inst->monitoringDates[i],0);
                if (payDate <= today) {
                    // past so already paid
                    fvFactors[i] = 0.0;
                } else {
                    // infrastructure will pv from payment date, so override for this case
                    // pvFromPaymentDate() can be overwriten in this class, but I didn't.
                    fvFactors[i] = discount->pv(payDate) / pvFromPaymentDate(); 
                }
                
                cashFlowAfterKO[i] = 0.0;
                if(inst->hasFixedLeg){
                    cashFlowAfterKO[i] -= koFix->getKOPV(inst->monitoringDates[i],discount);
//                    cashFlowAfterKO[i] -= inst->fixedLeg->getKOPV(inst->monitoringDates[i], 
//                                                                  discount, 
//                                                                  inst->koStubRule);
                }
                if(inst->hasFloater){
                    cashFlowAfterKO[i] -=  koFlt->getKOPV(inst->monitoringDates[i]);
//                    cashFlowAfterKO[i] -= inst->floater->getKOPV(inst->valueDate, 
//                                                                 inst->monitoringDates[i], 
//                                                                 discount, 
//                                                                 inst->koStubRule);
                }
                // getKOPV returns the value at monitoringDate.  
                // need to shift at last payment Date (Maturity + X), as well.
                // fvFacotrs[] are not possible to use here, because
                // cashFlow do not depend on the inst settlement but own paydates.
                cashFlowAfterKO[i] *= discount->pv(inst->monitoringDates[i]) / pvFromPaymentDate();
            }
        }
        catch (exception& e)
        {            
            throw ModelException(e, method, "Swap Leg has error in Schedule.");
        }
/*
/// Degbug ///        
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file (keiji)

    sprintf(fname, "c:\\temp\\debug.dat");
    DumpFile = fopen(fname, "w");
    if (DumpFile)
    {
        sprintf(buffer, "iMonitorDate, fvFactors, cashFlowAfterKO \n");
        fprintf(DumpFile, buffer);
        for (i=0; i<inst->monitoringDates.size(); i++)
        {
            sprintf (buffer, "%d, %f, %f \n", i, fvFactors[i], cashFlowAfterKO[i]);
            fprintf(DumpFile, buffer);
        }
        sprintf (buffer, "\n");
        fprintf (DumpFile, buffer);
        fclose(DumpFile);
        DumpFile = 0;
    }

///////////////
*/
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        static const string method = "TriggerECO::payoff";
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset, i;
        int    pathIdx = pathGen->getPathIndex();

        // 1. determine raw performances for each asset
        // 2. establish KO multipliers ("probs") for each asset
        // new.  Check Trigger on the node or not.
        // 3. apply any performance modifiers (call/put/etc)
        // 4. form "expected payoff per asset" & aggregate into "basket" - XXX what about "weights"?
        // 5. apply any modifier (call/put/etc) to basket
        // Not yet considered:
        // CashFlow

        // Initialise for potentially cached ones
        PricesTriggerECO& myPrices = static_cast<PricesTriggerECO&>(prices);
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

        double koTrigger;
        bool isConditionMet;
        int metAtStep;
        double conditionalCashFlow = 0.0;

        trigger->hitValueAndTimeAsBasket(pathGen,
                                         inst->triggerBasket,
                                         koTrigger,
                                         isConditionMet,
                                         metAtStep);

        if (pathGen->doingPast()){ 
            sumSoFar = sum;
            // This for BARRIER_LEVEL request. See RainbowKO.cpp for more complete comment
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }
            // calculate cashflow
            finaliseSimulationDates();
            // if past path is already touched, it's already triggered in past values.
            // those flag are used in only recordExtraOutput, to add or not PlainSwap.
            isTriggeredInPast = isConditionMet;
            iStepTrgInPast = metAtStep;
        }

        if (isConditionMet){
            //the swap leg value
            conditionalCashFlow = cashFlowAfterKO[metAtStep];
            //add the rebate value (bonus rebate)
            conditionalCashFlow += (inst->rebatePayAtHit ? fvFactors[metAtStep] : 1.0) 
                                    * rebate->getLevel(metAtStep);
        }

        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().

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
                basketParts[iAsset] = koFactors[iAsset] * assetPerfs[iAsset];
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
            bask /= nbAssets; // equal weights for now...

            // 5. apply any modifier (call/put/etc) to basket
            overallOption->apply();
            
            // Done...
            //prices.add(inst->notional * basket()); 
            prices.add(inst->notional * ( koTrigger*(basket()+redemption) +  (1.0-koTrigger)*conditionalCashFlow) ); 
        }
    }

    IMCPrices* createOrigPrices(int  nbIter,
                                        int  nbSubSamples,
                                        int  mode) {
        return new PricesTriggerECO((mode & CACHE_PRODUCT_BIT)? true: false, // use cache?
                             nbAssets,
                             nbIter, 
                             nbSubSamples);
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment     */
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
        static const string method = "TriggerECO::getVolInterp";

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
        // add plain cash flow.  
        bool addFloat = inst->hasFloater;
        bool addFixed = inst->hasFixedLeg;
        DateTime hitDate, payDate; //= inst->monitoringDates[inst->monitoringDates.size()-1];
        if (isTriggeredInPast){
            hitDate = inst->monitoringDates[iStepTrgInPast];
            payDate = settlement->settles(hitDate,0);
        }
        
        if (inst->withPlainSwap)
        {
            if (isTriggeredInPast){
                DateTime nextPayDate;
                if (inst->hasFloater){
                    if(inst->floater->getNextPayDate(hitDate, nextPayDate)){
                        if (nextPayDate < inst->valueDate){
                            addFloat = false;
                        }
                    }
                }
                if (inst->hasFixedLeg){
                    if(inst->fixedLeg->getNextPayDate(hitDate, nextPayDate)){
                        if (nextPayDate < inst->valueDate){
                            addFixed = false;
                        }
                    }
                }
            }

            double plainCashFlow = 0.0;
            plainCashFlow +=  addFloat? inst->floater->getPV(inst->valueDate, discount) : 0.0 ;
            plainCashFlow +=  addFixed? inst->fixedLeg->getPV(inst->valueDate, discount) : 0.0 ;
            results->storePrice(plainCashFlow * inst->notional + results->retrievePrice(),
                results->getCcyName());           
        }
        else{
            addFloat = addFixed = false;
        }

        // prepare for KNOWN CASH FLOW and Payment Dates
        CashFlowArraySP knownCfl (new CashFlowArray(0));
        CashFlowArraySP knownCflL (new CashFlowArray(0));
        if(addFloat){
            knownCflL = inst->floater->getCashFlowArray(inst->valueDate, inst->valueDate, discount);                
        }
        CashFlowArraySP knownCflF (new CashFlowArray(0)); ;
        if(addFixed){
            knownCflF = inst->fixedLeg->getCashFlowArray(DateTime(0,0));                
        }
        knownCfl = CashFlow::merge(knownCflF, knownCflL);
        // add KO cancelled leg.
        if(isTriggeredInPast){
            CashFlowArraySP koCflL (new CashFlowArray(0));
            CashFlowArraySP koCflF (new CashFlowArray(0));
            CashFlowArraySP rebCfl (new CashFlowArray(0));
            if (isTriggeredInPast){
                if(inst->hasFloater)
                    koCflL = koFlt->getKOCashFlowArray(hitDate, discount, inst->koStubRule, true);
                    //koCflL = inst->floater->getKOCashFlowArray(hitDate, discount, inst->koStubRule, true);
                if(inst->hasFixedLeg)
                    koCflF = koFix->getKOCashFlowArray(hitDate, discount, inst->koStubRule, true);
                    //koCflF = inst->fixedLeg->getKOCashFlowArray(hitDate, discount, inst->koStubRule, true);
                rebCfl->push_back(CashFlow(payDate, rebate->getLevel(iStepTrgInPast)));
            }                                                
            CashFlowArraySP koCfl(CashFlow::merge(koCflF, koCflL));
            knownCfl = CashFlow::merge(CashFlow::merge(koCfl, rebCfl),knownCfl);
        }
        
        // KNOWN CASH FLOW
        OutputRequest*  request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished())
        {   // add plain swap as known cash flow.  (Although, could be cancelled....) is it good?
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    knownCfl.get());   
        }

        // PAYMENT_DATES  
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            DateTimeArray payDateArray;
            for (int i=0;i<(*knownCfl).size();i++){
                payDateArray.push_back((*knownCfl)[i].date);
            }
            if(isTriggeredInPast)
                payDateArray.push_back(payDate);
            OutputRequestUtil::recordPaymentDates(control,results,&payDateArray);
        }
                   
        // BARRIER_LEVEL ...
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        const DateTime& today = getToday();
        const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
        // Only try to satisfy this request when have some past (so easy to get ref levels).
        if (request && (today >= lastRefDate)) {
            int iAsset;
//            DoubleArray assetSpots(nbAssets,0.0);
//            for(iAsset=0; iAsset<nbAssets; iAsset++) 
//                assetSpots[iAsset] = getMultiFactors()->assetGetSpot(iAsset)/histRefLevels[iAsset];
            // This operates on a per-asset basis
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                // Mostly delegate to Barrier class...
                // histRefLevels are to allow absolute barriers to be reported.
                BarrierLevelArraySP levels = 
                    barrier->reportLevels(today, histRefLevels[iAsset], iAsset);

                // calculate for up out triger.
                // currently, calc the barrier level by assuming one asset.  Not good as
                // no necessary Barrier Breach report would be generated.
                BarrierLevelArraySP trgLevels = 
                    trigger->reportLevels(today, histRefLevels[iAsset], iAsset);
//                BarrierLevelArraySP trgLevels = trigger->reportLevelsAsBasket(inst->triggerBasket,
//                                                                        today,
//                                                                        assetSpots,
//                                                                        histRefLevels[iAsset],
//                                                                        iAsset);
                // migrate two barrier levels
                for (int j=0;j<trgLevels->size();j++)
                    levels->push_back((*trgLevels)[j]);

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
IMCProduct* TriggerECO::createProduct(const MonteCarlo* model) const {
    // Validate the barriers based on the model
    barrierUnion->getBarrier()->validate(model, assets->NbAssets());
    triggerUnion->getBarrier()->validate(model, assets->NbAssets());

    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    simSeries->addDates(monitoringDates);
    return new TriggerECOMC(this, simSeries);
}

CClassConstSP const TriggerECO::TYPE = CClass::registerClassLoadMethod(
    "TriggerECO", typeid(TriggerECO), TriggerECO::load);

// * for class loading (avoid having header file) */
bool TriggerECOLoad() {
    return (TriggerECO::TYPE != 0);
}

DRLIB_END_NAMESPACE









