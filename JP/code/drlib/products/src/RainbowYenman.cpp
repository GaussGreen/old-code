//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RainbowYenman.cpp
//
//   Description : Swap with KO (NFBinary type KO).  Fixed Coupon Stream + Eq Link Cpn 
//                 based on Rainbow Basket.
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Format.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/IRebate.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/KOLiborLeg.hpp"
#include "edginc/EqLinkCashFlow.hpp"

// for Tree
#include "edginc/CallableEquityKOSwap.hpp"

DRLIB_BEGIN_NAMESPACE
class RainbowYenman: public GenericNFBase, 
                     virtual public IMCIntoProduct,
                     // for tree
                     virtual public FDModel::IIntoProduct,
                     virtual public LastSensDate 
{
public:
    static CClassConstSP const TYPE;
    friend class RainbowYenmanMC;
    friend class RainbowYenmanSVMC;
    
public:
    //---------------------------------------------//
    //  Class for Simple Barrier  (SimHit)         //
    //---------------------------------------------//
    class     BarSimHit: public CObject {
    public:
        static CClassConstSP const TYPE;  

        void validatePop2Object(){
            static const string method("RainbowYenman::BarSimHit::validatePop2Object");
           
            try {

            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        };

        BarrierScheduleSP makeBarrierSchedule() const{
            BarrierScheduleSP      barSch = BarrierScheduleSP(new BarrierSchedule(levels,
                                         true,                  // always isOut
                                         isUp,economicLevels,
                                         numHits,isHit,hitDate,
                                         "E",   // always European
                                         false,0.0,0.0,         // noSmoothing, no spread
                                         true,                          // simHit!!
                                         true));                        // isInternalMonitoringDates

            //barSch->validatePop2Object();

            return barSch;
        }
        
        BarSimHit():CObject(TYPE){}; 

        ScheduleSP getBarSched() const{
            return levels;
        };

        ScheduleSP getEcoBarSched() const{
            return economicLevels;
        };

    private:
        // data
//        BarrierScheduleSP      barSch;           // use barrier schedule class for real operation. $unregistered

        int                    numHits;           // How many hits trigger a breach
        mutable IntArraySP     isHit;             // has the barrier already been breached for each asset
        DateTime               hitDate;           // when it was breached

        ScheduleSP                  levels;            // the schedule of dates and levels
        bool                        isUp;
    
        ScheduleSP                  economicLevels;

   
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("RainbowYenman");
            REGISTER(RainbowYenman::BarSimHit, clazz);
            SUPERCLASS(CObject);
            FIELD(numHits, "number of hits");
            FIELD(isHit, "is hit ? per asset");
            FIELD(hitDate, "breach date");
            FIELD(levels, "barrier schedule");
            FIELD(economicLevels, "Economic barrier schedule");
            FIELD(isUp, "breach above ?");
            EMPTY_SHELL_METHOD(defaultBarSimHit);
        };
        
        static IObject* defaultBarSimHit(){
            return new RainbowYenman::BarSimHit();
        };
    };


    // validation
    void validatePop2Object(){
        static const string method = "RainbowYenman::validatePop2Object";
        GenericNFBase::validatePop2Object();
        try
        {
            if (hasFloater && !floater.get()){
                throw ModelException(method,
                                     "hasFloater = true, but floater is not found");
            }

            // to do :  what kind of validation do I need???
            // Make monitoringDate.
            // get the barrier sampling date from barrier union
            BarrierSP realBar = BarrierSP(barSimHit.makeBarrierSchedule());
            realBar->validatePop2Object();  // set up barrier Schedule for BarrierSimHit;
            ScheduleSP barSch = realBar->getBarrierSchedule();
            DateTimeArray barDates = barSch->getDateArray();

            DateTimeArray eqDates = eqCpn->getObservDates();

            if (matDates.size()>0){
                DateTime::ensureStrictlyIncreasing(matDates,"matDates must be strickly increasing",true);
            }
            // check input for inAdvance chase.
            if (inAdvance){
                if (matDates.size()<=0)
                    throw ModelException(method, "matDates is necessary when isAdvance = true");
                if (matDates[0] < eqDates[eqDates.size()-1])
                    throw ModelException(method, "matDates should be later than the last of equity observation dates.");
            }
            // painful but too much code implicitly depending on this in Barrier class
            // require all simulation dates are monitoring dates - so ...
            if (!DateTime::isSubset(eqDates, barDates)) {
                throw ModelException(method,
                                     "Require barrier schedule to be a subset of equity link coupon's observation dates");
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
        DateTimeArray eqObservDates = eqCpn->getObservDates();
        if (inAdvance)
            eqObservDates = DateTime::merge(eqObservDates, matDates);

        // get the barrier sampling date from barrier union
        BarrierSP realBar = BarrierSP(barSimHit.makeBarrierSchedule());
        realBar->validatePop2Object();
        ScheduleSP barSch = realBar->getBarrierSchedule();
        DateTimeArray barDates = barSch->getDateArray();

        // migrate bar dates & eq link coupon observation date
        return DateTime::merge(eqObservDates, barDates);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
    static const string method = "RainbowYenman::GetMarket";
    try{
        // parent
        GenericNFBase::GetMarket(model, market);
        // relevant elements of self
        rebate->getMarket(model, market.get());
        if(hasFloater)
        {
            floater->setCouponCurve(discount); // in case floater gets coupon curve from inst
            floater->getMarket(model, market.get());
        }

        /** Set Up CallableEquitKOSwap products */
        if (CTree1f::TYPE->isInstance(model)){
            setOriginalCKS();
            originalCKS->GetMarket(model, market);
        }

    } catch (exception& e){
        throw ModelException(e, method);
        }
    }

protected:
    /// fields 
    IDoubleArrayModifierMakerSP  matPerformance;    // option at maturity
    IAggregateMakerSP            assetBasket;       // pct/rainbow/product basket etc of ...
    IRebateMakerSP               rebate;            // the "opposite" side of the payoff from 'assetBasket'
    bool                         hasFloater;    
    
    EqLinkCashFlowSP             eqCpn;             //interface
    bool                         inAdvance;         // [false (default)]:Eq Observ and KO date is same.  [true]:Eq observ Dates is earlier, and KO sample is later.
    DateTimeArray                matDates;          // inAdvance = "Y", then need to specify the matDates.

    RainbowYenman::BarSimHit     barSimHit;      // Barrier Union Class.
    
    double                       redemption;        // redemption amount paid at exit (maturity or early terminate).
    LiborLegSP                   floater;
    
private:
    RainbowYenman(): GenericNFBase(TYPE) {} // for reflection
    RainbowYenman(const RainbowYenman& rhs);     // not implemented
    RainbowYenman& operator=(const RainbowYenman& rhs); // not implemented

    //-----------------------------for Tree--------------------------------------------------//
    void setOriginalCKS(){
    static const string method = "RainbowYenman::GetMarket";
    try{
        // validation
        int nbAssets = assets->NbAssets();
        if (nbAssets != 1){
            throw ModelException(method, "Tree can have only 1 asset.  This inst has " 
                                        + Format::toString(nbAssets) + " assets.");
        }

        // prepare two importand date array, pay dates & observation dates.
        DateTimeArray eqPayDates = eqCpn->getPaymentDates();
        if (eqPayDates.size() <=0 )
            throw ModelException(method, "Equity Linked coupon size should not be empty!!");
        DateTimeArray   eqObsDts = eqCpn->getObservDates();
        
        //--------------- Setting up CKS ---------------//

        ///////////////////////////////////////////
        //  convert from refLevel to CashFlowArray (average in).
        DateTimeArray   startDates = refLevel->getAllDates();
        DoubleArray     smplLvls(0);
        if (startDates.size()>1){
            throw ModelException(method, "reference dates have " 
                                        + Format::toString(startDates.size())
                                        + " dates.  Should be one single date."
                                        + "  Average-in is not available with tree.");
        }
        double  initialSpot = 0.0;  
        if (valueDate>=startDates[0]){  // started case
            DoubleArray refLevels = pastValues->getPastValues(startDates, 0, valueDate);
            initialSpot = refLevels[0];
            smplLvls = pastValues->getPastValues(eqObsDts, 0, valueDate);
        }else{                          // forward starting => cks will take care about this.            
        }
        CashFlowArray samples = CashFlowArray(eqObsDts.size());
        int smplSize = smplLvls.size();
        for (int i=0;i<eqObsDts.size();i++){
            if (i<smplSize)
                samples[i] = CashFlow(eqObsDts[i], smplLvls[i]);
            else
                samples[i] = CashFlow(eqObsDts[i], 0.0);
        }
        CashFlowArray avgIn = CashFlowArray(1, CashFlow(startDates.back(), initialSpot));

        ///////////////////////////////////////////
        // dummy null objects
        KOFixedLegMaker nullKOFixed;
        CallableEquityKOSwap::Trigger nullTrigger;
        CallableEquityKOSwap::CallSchedule nullCallSched;

        // koStubRule = "Bond" type
        KOStubRuleSP koRule = KOStubRuleSP(new KOStubRule("B",false,"AsSchedule"));     // this settlemetn timing rule is not used, yet.
        
        ///////////////////////////////////////////
        //making EqCpn by adding koRule.
        DateTime initialAccrueDate = eqObsDts[0].rollDate(-1);   // accrue start just one day before.
        EqCpnKOMakerSP koELC = EqCpnKOMakerSP(new EqCpnKOMaker(eqCpn,initialAccrueDate,koRule, inAdvance));

        ///////////////////////////////////////////
        //making koLiborLeg using above (same as eqCpn) koRule.
        bool hasFloater = (!!floater && floater->getSize() > 0);
        KOLiborLegMaker koLiborLeg;
        if (hasFloater)
            koLiborLeg = KOLiborLegMaker(floater, koRule);

        // prepare barrier as ScheduleSP
        ScheduleSP barSch = barSimHit.getBarSched();
        ScheduleSP ecoBarSch = barSimHit.getEcoBarSched();

        ///////////////////////////////////////////
        //making rebate.  
        // When KO, return redemption.  
        // So barrier schedule should not have any KO date 
        // after the last coupon accrue started!! (in LiborLeg and EqCpn)
        DateTime lastAccrueStart = eqObsDts[eqObsDts.size()-2];
        if (hasFloater){
            if (lastAccrueStart > floater->AccrualDates[floater->getSize()-1])
            lastAccrueStart = floater->AccrualDates[floater->getSize()-1];
        }
        // convert from IRebate to ScheduleSP        
        ScheduleSP rebSch;
        bool hasRebate = false;
        if (!rebate.get()){
            rebSch = ScheduleSP(0); // null
        }else{
            DateTimeArray rebDates = barSch->getDates();
            DoubleArray   rebValues(rebDates.size(),redemption);
            string        rebInterpType = "N";
            if (FlatRebateMaker::TYPE->isInstance(rebate)){
                FlatRebateMaker* flatRebate = dynamic_cast<FlatRebateMaker*>(rebate.get());
                double rebVal = flatRebate->getRebateAmount();
                for (int i=0; i<rebDates.size(); i++){
                    rebValues[i] += rebVal;
                }
            }else if (ScheduleRebateMaker::TYPE->isInstance(rebate)){
                ScheduleRebateMaker* rebMaker = dynamic_cast<ScheduleRebateMaker*>(rebate.get());
                ScheduleSP scheRebate = rebMaker->makeSchedule();
                rebInterpType = scheRebate->getInterp();
                rebDates = scheRebate->getDates();                        
                DoubleArray tmp = scheRebate->getValues();
                for (int i=0; i<tmp.size(); i++){
                    rebValues[i] += tmp[i];
                }
            }else{
                throw ModelException(method, "Invalid type of rebate.  Only Flat or Schedule Rebate "
                                                "is acceptable for Tree." );
            }
            rebSch = ScheduleSP(new Schedule(rebDates, rebValues, rebInterpType));
        }       
                                                    
        ///////////////////////////////////////////
        //making trigger
        CallableEquityKOSwap::Trigger upTrigger(true,    // isUp
                                                barSch,
                                                ecoBarSch, 
                                                false,          //IntraDayMonitor
                                                true,           // hasRebate
                                                rebSch,
                                                "NO_SMOOTHING",     //trigger.smoothingType
                                                0.0,                //trigger.peakDelta
                                                eqPayDates);
        
        ///////////////////////////////////////////
        // making finalPerf
        //Generic1Factor::GetMarket(model, market);   // Need to get market for finalPerf.
        DateTime fpSettleDate = eqPayDates.back();            

        if (matDates.size()>1){
            // lazy to allow average-out, because AvgOutPerf is not using generalized perf type....
            //DateTime::ensureStrictlyIncreasing(matDates,"matDates must be strickly increasing",true);
            throw ModelException(method, "matDates should be single date when you use Tree.");
        }

        // check input for inAdvance chase.
        if (inAdvance){
            if (matDates.size()<=0)
                throw ModelException(method, "matDates is necessary when isAdvance = true");
            if (matDates[0] < eqObsDts.back())
                throw ModelException(method, "matDates should be later than the last of equity observation dates.");
        }

        DateTime matDate = inAdvance ? matDates.back() : eqObsDts.back();
        if (fpSettleDate < matDate){
            throw ModelException(method, "settle date for FinalPerf[" 
                                        + fpSettleDate.toString() +
                                            "] should be later or equal than last of the matDates["
                                        + matDate.toString() + "] .");
        }
        // don't use instrument's settlement rule.
        InstrumentSettlementSP myFinalSettle = smartPtr<CashSettleDate> (new CashSettleDate(fpSettleDate));        
        VanillaPerfMakerSP finalPerf = VanillaPerfMakerSP(new VanillaPerfMaker(matPerformance, matDate));

        // Create a TargetRedemptionNote to make use of common validation and for pricing
        bool oneContract = false;

        InstrumentSettlementSP  premiumSettle;    /* When premium is paid */

        // convert from assets (multi asset) to CAssetWrapper (single asset) 
        CAssetWrapper myAsset = CAssetWrapper(assets->getName(0));

        // get ccyTreatment from the first asset.
        const IAsMultiFactors* asMulti = dynamic_cast<const IAsMultiFactors*>(assets.get());
        if (!asMulti){
            throw ModelException(method, "assets as object of type "
                                        + assets.get()->getClass()->getName()
                                        +" is not supported for Tree");
        }
        IMultiFactorsSP mAsset(asMulti->asMultiFactors());
        string ccyTreatment = mAsset->assetGetCcyTreatment(0);

//            const IMultiFactors* ma = dynamic_cast<const IMultiFactors*>(assets.get());
//            ccyTreatment = ma->assetGetCcyTreatment(0);

        originalCKS = CallableEquityKOSwapSP(
                    new     CallableEquityKOSwap(
                                    // Generic1Factor
                                    //valueDate, oneContract, notional, initialSpot, 
                                    //ccyTreatment,instSettle, premiumSettle, asset, discount,                                    
                                    valueDate, oneContract, notional, initialSpot, 
                                    ccyTreatment,instSettle, premiumSettle, myAsset, discount,                                    
                                    // for CKS
                                    true,                           //hasEqLeg,
                                    koELC,
                                    finalPerf,                 // zero performance
                                    redemption,    
                                    false,                           // hasFixedLeg,
                                    hasFloater,                          // hasFloater,
                                    nullKOFixed,                     // 
                                    koLiborLeg,                           // koFloater
                                    true,                           // hasUpTrigger,
                                    false,                          // hasLowTrigger,
                                    upTrigger,
                                    nullTrigger,                           // lowTrigger
                                    false,                          // hasCallSchedule,
                                    nullCallSched,                        // callSchedule,    
                                    true,                           // is average in? => always use avgIn.
                                    avgIn,                                // averageIn.avgIn,
                                    false,                           //averageIn.useIniSpotAsRefLvl,
                                    samples,
                                    1,                              //numOfAvgInSV,
                                    myFinalSettle
#ifdef  TREE_THETA_CAP
                                    ,
                                    0.0,                            //thetaSmoothThrehold,
                                    false                          //isPositiveThetaSmooth
#endif
                                    ));

        originalCKS->validatePop2Object();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
    }
    //-----------------------------End of for Tree----------------------------------------------//

    static IObject* defaultRainbowYenman(){
        return new RainbowYenman();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RainbowYenman, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultRainbowYenman);
        FIELD(matPerformance, "Overall Perf Modifier");
        FIELD(assetBasket, "How to aggregate asset perfs");
        FIELD(eqCpn, "Equity Linked Cash Flow");
        FIELD(inAdvance, "[false(default)]:When KO occurs, the sampled coupon on the same date will be paid.  [true]When KO occurs, the sampled coupon is cancelled.");
        FIELD(barSimHit, "KO barrier.  Sampling on M (hit asset) of N-asset");
        FIELD(redemption, "redemption amount at terminate. regardless to be paid early term or not.");
        FIELD(rebate, "Opposite side of payoff to assetBasket");
        FIELD_MAKE_OPTIONAL(rebate);
        FIELD(hasFloater, "has Libor Leg?")
        FIELD(floater, "libor leg");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(matDates, "maturity performance observation date.  Necessary only for isAdvance = true");
        FIELD_MAKE_OPTIONAL(matDates);
        FIELD(originalCKS, "Internal CallableEquityKOSwap");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(originalCKS);  
    }

    //-----------------------------for Tree--------------------------------------------------//
private: 
    // for internal use  
    CallableEquityKOSwapSP       originalCKS;    

public:

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const{
        originalCKS->Validate();
        return originalCKS->createProduct(model);
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const{
        DateTime end = eqCpn->getPaymentDates().back();
        return end;    
    }

    // make originalCKS, here.  createProd is too late.
    bool priceDeadInstrument(CControl* control, CResults* results) const{
        return originalCKS->priceDeadInstrument(control, results);
    }

    bool sensShift(Theta* theta){
        return GenericNFBase::sensShift(theta);
// no need to call CallableEquityKOSwap's because it's called automatically... 
//        if (!!originalCKS.get())
//            return originalCKS->sensShift(theta);            
    }


    //-----------------------------End of for Tree----------------------------------------------//
};

//////////////////////////////////////////////////////////////////////////////////

/* MC product class for RainbowYenman */
class RainbowYenmanMC: public IMCProduct,
                   virtual public IMCProductLN,
                   virtual public IMCProductImplied{
private:
    const RainbowYenman*          inst;
    BarrierSP                   barSP;       // non-const note
    DateTimeArray               monitoringDates;    
    DateTimeArray               payDates;       // eqLink pay dates

    int                       nbAssets;         // convenient

    SimpleDoubleArray         assetComps;
    SimpleDoubleArray         perfForELC;
    IntArray                  eqObsMap;
    IntArray                  obsIndex;
    IAggregateSP              assetBasket;
    TrivialDoubleArray        basket;  // just a double ... being the aggregation of per-asset "perfs"
    IRebateSP                 rebate;
    IDoubleArrayModifierSP    overallOption;

    IntArray                  avgMap;    // to track averaging dates
    DoubleArray               sum;       // [nbAssets], saves alloc later
    DoubleArray               fvFactors; // if rebate at hit the this FVs to mat. 
                                         // Past rebates drop out via fvFactor=0.0

    // for past
    DoubleArray               sumSoFar;  // [nbAssets]
    // From a past pricing; used for BARRIER_LEVEL request
    DoubleArray               histRefLevels; // [nbAssets] 

    EqLinkCashFlowSP         eqCpn;
    EqCpnKOLegSP             eqKOLeg;
    
    double                   redemption;    
    LiborLegSP               floater;

    // need to daclare here to avoid loosing the object entire life.
    KOLiborLegSP             koFlt;         // for intenal usage
            
    DoubleArray              liborValues;   // libor value when it's knocked out at each time step.
    double                   plainLibor;    // plain libor value
    bool                     hasFloater;    
    
    CashFlowSP               knownRebate;   // if already touched, for known cash flow

public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to RainbowYenman) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    RainbowYenmanMC(const RainbowYenman*          inst,
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
        perfForELC(inst->eqCpn->getPaymentDates().size(),0.0),
        basket(0.0),
        barSP(inst->barSimHit.makeBarrierSchedule()),
        rebate(inst->rebate->getRebate(simSeries->getAllDates(), discount)),
        sum(nbAssets, 0.0),
        sumSoFar(nbAssets, 0.0),
        histRefLevels(nbAssets, 0.0){

        int i=0,j=0;
        hasFloater = inst->hasFloater;
        
        // overwrite paymentDate, because FixedLeg or EqLink could have paydates which
        // is longer then instSettle.
        DateTimeArray eqPayDates = inst->eqCpn->getPaymentDates();
        DateTime lastPayDate = eqPayDates.back();
        if (lastPayDate > paymentDate)
            paymentDate = lastPayDate;
        if (hasFloater){
            lastPayDate = inst->floater->PayDates.back();
            if (lastPayDate > paymentDate)
                paymentDate = lastPayDate;
        }
        
        // -----------  Create Monitoring Array ----------//
        // make the monitoringDates by merging barrier and eqCpn's monitoring daes
        // ?????????????????????
        // However, souceBarDates are assumed (validate) to be subset of eqObservDates.
        // Thus, for time being, there is not much meaning to make monitoringDate by merge.        
        barSP->validatePop2Object();    //setup internal barrierSchedule.
        ScheduleSP barSch = barSP->getBarrierSchedule();            
        DateTimeArray barDates = barSch->getDateArray();
        DateTimeArray eqObservDates = inst->eqCpn->getObservDates();
        if (inst->inAdvance)
            eqObservDates = DateTime::merge(eqObservDates, inst->matDates);
        monitoringDates = DateTime::merge(eqObservDates, barDates);
        int nbSteps = monitoringDates.size();

        // -----------  Set up Barrier Class ----------//
        barSP->validate(nbAssets);
        // create an interpolated barrier to match the monitoring dates
        barSP->createInterpBarrier(inst->valueDate, monitoringDates);


        // -----------  Set Up Basket for Final Option / Range Sampling ----------//
        // Tie the pieces together :
        // First creating the per-asset performances, then aggregating them into a single number and
        // finally turning that into an overall perf
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
        overallOption = IDoubleArrayModifierSP(inst->matPerformance->getModifier(&basket));


        //default koRule for all legs.
        KOStubRuleSP koRule(new KOStubRule("B",false,"AsSchedule"));    //koStubRule,isAccrueUpToSettle,payTimingRule
        
        // -----------  Make EqCpnKOLeg class ----------//
        eqKOLeg = EqCpnKOLegSP(new EqCpnKOLeg(inst->eqCpn,perfForELC,inst->valueDate,koRule,inst->inAdvance));
        eqObsMap = eqKOLeg->setMonStep(simSeries->getAllDates());
        eqKOLeg->setPayDatesAndDiscFacts(monitoringDates, discount, paymentDate);
        // make a IndexArray to know the position of iObs against iStep.
        int iObs = 0;
        obsIndex.resize(eqObsMap.size()-1);
        for (int iStep = 0; iStep < eqObsMap.size()-1;  iStep++){
            obsIndex[iStep] = iObs;
            if (eqObsMap[iStep] == 0)
                iObs++;
        }
        
        
        // -----------  Set up Product.... ----------//        
        double df =  pvFromPaymentDate();                    
        redemption = inst->redemption;

        // -----------  Set up pvFactors for Rebate ----------//
        const DateTime& today = getToday();
        payDates = eqKOLeg->getPayDatesOnTimeLine();
        int nbCpns = payDates.size();
        // currently, I assume # of monitoring date is same to # of coupon pay dates of EqCpn.
        // If we want to allow the Average sample as each coupon performance, we need think about it.
        if (nbSteps != nbCpns)
            throw ModelException("RainbowYenmanMC", "number of monitoring dates should be same to coupon payment dates in eq link coupon.");
        fvFactors = DoubleArray(nbCpns, 1.0); // default is pay at mat
        for(i=0; i<fvFactors.size(); i++) {
            if (payDates[i] <= today) {
                // past so already paid
                fvFactors[i] = 0.0;
            } else {
                // infrastructure will pv from payment date, so override for this case
                // pv(payDate[i]) or pv(today, payDates[i]) are no difference.
                // for easy to read the code, I used pv(today, payDate).
                fvFactors[i] = discount->pv(today, payDates[i]) / df; 
            }
        }

        // -----------  Set up Libor Leg ----------//
        liborValues = DoubleArray(nbSteps,0.0);     //initialize by zero value.
        if (hasFloater){
            floater = inst->floater;
            KOSettleSP koSettle(new KOSettle(koRule,inst->instSettle,
                                              floater->PayDates,
                                              floater->AccrualDates, 
                                              false));       // accrue up to Settle
            
            koFlt = KOLiborLegSP(new KOLiborLeg(floater,koSettle,inst->valueDate,discount));

            // when all samples are in past, df = 0.0.
            double factor = df >0.0 ? 1/df : 0.0;   
            
            // get plain libor leg value and value at each time line
            // values are converted as of last settlement date.
            plainLibor = koFlt->getPV(inst->valueDate, discount)*factor;
            liborValues = koFlt->getPVsAlongKOTimeStep(monitoringDates,factor);      

            // prepare KNOWN_CASH_FLOW.  plain libor leg at first level.
            koFlt->makeKnownCashFlows(inst->valueDate, false);
        }


        // -----------  prepare for Average-In ----------//
        bool isTrivial;
        avgMap = DateTime::createMapping(simSeries->getAllDates(),
                                         eqObservDates,
                                         isTrivial);

    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        static const string method = "RainbowYenman::payoff";
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
        iStep = beginIdx + eqObsMap[beginIdx];
        if (iStep<endIdx){
            int iObs = obsIndex[iStep];
            for (; iStep<endIdx; iStep += eqObsMap[iStep+1]+1, iObs ++) {
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetComps[iAsset] = pathGen->Path(iAsset, 0)[iStep]/pathGen->refLevel(iAsset, 0);
                }
                perfForELC[iObs] = assetBasket->aggregate();            
            }
        }
        if (inst->inAdvance && endIdx > 0){
            // need to calc the final performance.
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                assetComps[iAsset] = pathGen->Path(iAsset, 0)[endIdx -1]/pathGen->refLevel(iAsset, 0);
            }
        }

        // 2. establish KO multiplier
        // Note this is outside the "if past" clause, which allows
        // the barrier to handle its own past.
        double koFactor;
        bool isConditionMet;
        int metAtStep;
        barSP->hitValueAndTime(pathGen,
                                 koFactor,
                                 isConditionMet,
                                 metAtStep);

        int atStep = isConditionMet ? metAtStep : endIdx;
        atStep += inst->inAdvance ? 0 : 1;      // inAdvance, the coupon, whose sample is KO date, is cancelled.
        //capped so as to avoid ABR error.
        //"atStep" is used for only eqCashFlow
        //isAdvance = true casee, it won't be > endIdx.
        //isAdvance = false case, the last KO event is not sensitive to eqLeg results
        // as stub Rule = "B". 
        atStep = atStep > endIdx ? endIdx : atStep; 

        if (pathGen->doingPast()){ 
            sumSoFar = sum;
            // This for BARRIER_LEVEL request. See RainbowKO.cpp for more complete comment
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }
            // prepare the KNOWN_CASH_FLOW by looking at the past values.
            eqKOLeg->makeKnownCashFlow(beginIdx, atStep, nbAssets, hasFuture());
            // re-set libor KNOWN_CASH_FLOW when it's already breached.  
            if (hasFloater && isConditionMet)
                koFlt->makeKnownCashFlows(monitoringDates[metAtStep], isConditionMet);
            if (isConditionMet) //add rebate + redemption
                knownRebate = CashFlowSP(new CashFlow(payDates[metAtStep],rebate->getLevel(metAtStep)+redemption));
        }


        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().

            // 3(new). calculate eqLink coupon.  
            // eqKOLeg class has value as of today.
            //double eqlinkPay = eqKOLeg->getValue(pathGen, atStep, assetBasket, assetComps)/df;            
            double eqlinkPay = eqKOLeg->getValue(beginIdx, atStep, nbAssets);
            
            // 4. aggregate into a basket
            basket() = assetBasket->aggregate();

            // 5. apply any modifier (call/put/etc) to basket
            overallOption->apply();

            // 6. form "expected payoff" combining overall option and rebate
            double rebateValue = isConditionMet? 
                (fvFactors[metAtStep] * (rebate->getLevel(metAtStep) + redemption)) : 0.0;
            double payoff = koFactor * (basket() + redemption) + (1.0 - koFactor) * rebateValue + eqlinkPay;
            if (hasFloater)
                payoff +=  isConditionMet ? liborValues[metAtStep] : plainLibor ;
            
            // Done...
            prices.add(inst->notional * payoff); 
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // only for European, so nothing to do.
    }

    /** Use this opportunity to do any Implied driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustments */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{
        // only for European, so nothing to do.
    }

    // for the LogNormal path generator
    // XXX since there could be several interp levels needed we should
    // XXX defer construction of reqarr to a method on the barrier. TBD
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "RainbowYenman::getVolInterp";

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
        static const string method = "RainbowYenman::recordExtraOutput";
        try{
            CashFlowArraySP knownCfl = CashFlowArraySP(new CashFlowArray(0));

            // KNOWN CASH FLOW
            // not yed implemented after Maturity...
            OutputRequest* request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished())
            {   // add plain swap as known cash flow.  (Although, could be cancelled....) 
                //CashFlowArraySP rng = CashFlowArraySP(new CashFlowArray(eqKOLeg.getKnownCashFlows()));
                CashFlowArraySP rng = eqKOLeg->getKnownCashFlows();
                if (hasFloater)
                    knownCfl = CashFlow::merge(rng,koFlt->getKnownCashFlows());
                else 
                    knownCfl = rng;
                // when KO occurs, knownRebate is exiting and it's already include the redemption.
                // If all simulation date are in past, then need to add redemption at final.
                if (!!knownRebate)  
                    knownCfl->append(knownRebate);
                else if (!hasFuture()){  
                    rng->back().amount += redemption;
                }
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
                DateTimeArray listPayDates = eqKOLeg->getPayDates();
                DateTimeArray payDateArray = DateTime::merge(listPayDates,cashPayDates);
                OutputRequestUtil::recordPaymentDates(control,results,&payDateArray);
            }

            // BARRIER_LEVEL ...
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
                        barSP->reportLevels(today, histRefLevels[iAsset], iAsset);
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


};


//////////////////////////////////////////////////////////////////////////////////

/* MC product class for RainbowYenman state vars - with isRefAvgOut false */
class RainbowYenmanSVMC : public MCProductClient,
                  virtual public IMCProductLN{
private:
    const RainbowYenman*          inst;
    BarrierSP                   barSP;       // non-const note
    DateTimeArray               monitoringDates;    
    DateTimeArray               payDates;       // eqLink pay dates

    int                       nbAssets;         // convenient

    SimpleDoubleArray         assetComps;
    SimpleDoubleArray         perfForELC;
    IntArray                  eqObsMap;
    IntArray                  obsIndex;
    IAggregateSP              assetBasket;
    TrivialDoubleArray        basket;  // just a double ... being the aggregation of per-asset "perfs"
    IRebateSP                 rebate;
    IDoubleArrayModifierSP    overallOption;

    IntArray                  avgMap;    // to track averaging dates
    DoubleArray               sum;       // [nbAssets], saves alloc later
    DoubleArray               fvFactors; // if rebate at hit the this FVs to mat. 
                                         // Past rebates drop out via fvFactor=0.0
    
    // for past
    DoubleArray               sumSoFar;  // [nbAssets]

    EqLinkCashFlowSP         eqCpn;
    EqCpnKOLegSP             eqKOLeg;
    
    double                   redemption;    

    // need to daclare here to avoid loosing the object entire life.
    KOLiborLegSVSP           koFltSV;         // for intenal usage
    KOStubRuleSP             koRule;
            
    double                   plainLibor;    // plain libor value
    
    CashFlowSP               knownRebate;   // if already touched, for known cash flow

    // specific for SVMC
    DoubleArray             refLevel;
    
    SVGenSpotSP                  spotGen;      //!< Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   //!< Ref level state variable
    SVGenDiscFactorSP            matDfGen;        //!< Generator for discount factors
    SVDiscFactorSP         matDfSV;         //!< Df state variable
    LiborLeg::LiborLegSVGenSP liborLegGen;  //!< Generator for Libor flows
    LiborLeg::LiborLegSVSP    liborLegSV;   //!< Libor flows state variable
    SVGenDiscFactorSP            couponDfGen;  //!< Generator for discount factors from coupon dates
    SVDiscFactorSP         couponDfSV;   //!< Df state variable

    SVGenBarrierHVTSP            barrierHVTGen; //!< Barrier generator
    SVGenBarrierHVT::IStateVarSP  barrierHVTSV;  //!< Barrier state variable

    // for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS
    DateTimeArray            instPaymentDates;          
    OutputRequestUtil::KnownCashFlows knownCFs;
    int                      iFirstUnKnownFloater;  //first unknown floater cash flow.

protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "RainbowYenmanSVMC::pathGenUpdated";

        try {
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            matDfSV = matDfGen->getSVDiscFactor(newPathGen);
            couponDfSV = couponDfGen->getSVDiscFactor(newPathGen);
            if (inst->hasFloater) {
                liborLegSV = liborLegGen->getLiborLegSV(liborLegSV, newPathGen);
            }
            barrierHVTSV = barrierHVTGen->getHitValueTimeSV(barrierHVTSV, newPathGen);

            // -----------  Set up Libor Leg ----------//
            if (inst->hasFloater){
                // Build KOLiborLegSV, using the latest SV Generator
                KOSettleSP koSettle(new KOSettle(koRule,inst->instSettle,
                                                inst->floater->PayDates,
                                                inst->floater->AccrualDates, 
                                                false));       // accrue up to Settle                
                
                koFltSV = KOLiborLegSVSP(new KOLiborLegSV(liborLegSV,koSettle,inst->valueDate,discount));

                // prepare KNOWN_CASH_FLOW.  plain libor leg at first level.
                // koFltSV->makeKnownCashFlows(inst->valueDate, false);
            }

        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };
   
public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(matDfGen.get());
        svCollector->append(couponDfGen.get());
        svCollector->append(barrierHVTGen.get());    // barrier structure
        if (inst->hasFloater) {
            svCollector->append(liborLegGen.get());
        }
    }
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    RainbowYenmanSVMC(const RainbowYenman*            inst,
               const SimSeriesSP&       simSeries):
        MCProductClient(inst->assets.get(),
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
        perfForELC(inst->eqCpn->getPaymentDates().size(),0.0),
        basket(0.0),
        barSP(inst->barSimHit.makeBarrierSchedule()),
        rebate(inst->rebate->getRebate(simSeries->getAllDates(), discount)),
        sum(nbAssets, 0.0),
        sumSoFar(nbAssets, 0.0),
        spotGen(new SVGenSpot(simSeries)),
        refLevel(nbAssets, 0.0),
        redemption(inst->redemption),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())),
        matDfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())){

        int i=0,j=0;
        
        // overwrite paymentDate, because FixedLeg or EqLink could have paydates which
        // is longer then instSettle.
        DateTimeArray eqPayDates = inst->eqCpn->getPaymentDates();
        DateTime lastPayDate = eqPayDates.back();
        if (lastPayDate > paymentDate)
            paymentDate = lastPayDate;
        if (inst->hasFloater){
            lastPayDate = inst->floater->PayDates.back();
            if (lastPayDate > paymentDate)
                paymentDate = lastPayDate;
            // Build the SV Generator
            liborLegGen = LiborLeg::LiborLegSVGenSP(inst->floater->createLiborLegSVGen(inst->discount.getSP()));
        }
        
        // -----------  Create Monitoring Array ----------//
        // make the monitoringDates by merging barrier and eqCpn's monitoring daes
        // ?????????????????????
        // However, souceBarDates are assumed (validate) to be subset of eqObservDates.
        // Thus, for time being, there is not much meaning to make monitoringDate by merge.        
        barSP->validatePop2Object();    //setup internal barrierSchedule.
        ScheduleSP barSch = barSP->getBarrierSchedule();            
        DateTimeArray barDates = barSch->getDateArray();
        DateTimeArray eqObservDates = inst->eqCpn->getObservDates();
        if (inst->inAdvance)
            eqObservDates = DateTime::merge(eqObservDates, inst->matDates);
        monitoringDates = DateTime::merge(eqObservDates, barDates);
        int nbSteps = monitoringDates.size();

        // -----------  Set up Barrier Class ----------//
        barSP->validate(nbAssets);
        // create an interpolated barrier to match the monitoring dates
        barSP->createInterpBarrier(inst->valueDate, monitoringDates);

        // Convert barrier to barrier per asset generator
        barrierHVTGen = barSP->convertBarrier(monitoringDates, barSch->lastDate(), 
                                            getToday(), refLevelGen, getMultiFactors(), 0);

        // -----------  Set Up Basket for Final Option / Range Sampling ----------//
        // Tie the pieces together :
        // First creating the per-asset performances, then aggregating them into a single number and
        // finally turning that into an overall perf
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
        overallOption = IDoubleArrayModifierSP(inst->matPerformance->getModifier(&basket));


        //default koRule for all legs.
        koRule = KOStubRuleSP(new KOStubRule("B",false,"AsSchedule"));    //koStubRule,isAccrueUpToSettle,payTimingRule
        
        // -----------  Make EqCpnKOLeg class ----------//
        eqKOLeg = EqCpnKOLegSP(new EqCpnKOLeg(inst->eqCpn,perfForELC,inst->valueDate,koRule,inst->inAdvance));
        eqObsMap = eqKOLeg->setMonStep(simSeries->getAllDates());
        eqKOLeg->setPayDatesAndDiscFacts(monitoringDates, discount, Today); // discount factor is to today, not final date.
        // make a IndexArray to know the position of iObs against iStep.
        int iObs = 0;
        obsIndex.resize(eqObsMap.size()-1);
        for (int iStep = 0; iStep < eqObsMap.size()-1;  iStep++){
            obsIndex[iStep] = iObs;
            if (eqObsMap[iStep] == 0)
                iObs++;
        }

        // -----------  Set up pvFactors for Rebate ----------//
        // rebate payment date is always same to eqLink Coupon Dates!!
        payDates = eqKOLeg->getPayDates();
        couponDfGen = SVGenDiscFactorSP(new SVGenDiscFactor(inst->valueDate,
                                                    inst->discount.getSP(),
                                                    payDates));

        // -----------  Set up Libor Leg ----------//
        if (inst->hasFloater){
            // Build the SV Generator
            // for SV and KOLiborLegSP will be generated at UpDate.
            liborLegGen = LiborLeg::LiborLegSVGenSP(inst->floater->createLiborLegSVGen(inst->discount.getSP()));
        }

    };
   
public:
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form
        barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("RainbowYenmanSVMC::payoff");
        try {

            // Same begin & end for all assets, so read from the first
            const SVPath& path = spotSV->path(0/*iAsset*/);
            int    beginIdx = path.begin();
            int    endIdx   = path.end();

            int    iAsset, iStep;

            // 1. determine raw performances for each asset
            // 2. establish KO multiplier ("prob") for the basket (to be)
            // 3. (new). calculate Range Acrrue Payment..
            // 4. aggregate into a basket
            // 5. apply any modifier (call/put/etc) to basket
            // 6. form "expected payoff" combining overall option and rebate
            sum = sumSoFar;
        
            for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
                refLevel[iAsset] = refLevelSV->refLevel(iAsset);
            }

            // 1a. determine raw performances for each asset
            iStep = beginIdx + eqObsMap[beginIdx];
            if (iStep<endIdx){
                int iObs = obsIndex[iStep];
                for (; iStep<endIdx; iStep += eqObsMap[iStep+1]+1, iObs ++) {
                    for(iAsset=0; iAsset<nbAssets; iAsset++) {
                        const SVPath& path = spotSV->path(iAsset);            
                        assetComps[iAsset] = path[iStep]/refLevel[iAsset];
                    }
                    perfForELC[iObs] = assetBasket->aggregate();            
                }
            }
            if (inst->inAdvance && endIdx > 0){
                // need to calc the final performance.
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetComps[iAsset] = spotSV->path(iAsset)[endIdx -1]/refLevel[iAsset];
                }
            }

            // 2. establish KO multiplier
            // Note this is outside the "if past" clause, which allows
            // the barrier to handle its own past.
            double koFactor;
            bool isConditionMet;
            int metAtStep;
            isConditionMet = barrierHVTSV->hitValueTime(koFactor, metAtStep);

            int atStep = isConditionMet ? metAtStep : endIdx;
            atStep += inst->inAdvance ? 0 : 1;      // inAdvance, the coupon, whose sample is KO date, is cancelled.
            //capped so as to avoid ABR error.
            //"atStep" is used for only eqCashFlow
            //isAdvance = true casee, it won't be > endIdx.
            //isAdvance = false case, the last KO event is not sensitive to eqLeg results
            // as stub Rule = "B". 
            atStep = atStep > endIdx ? endIdx : atStep; 

            if (doingPast()){ 
                sumSoFar = sum;
                // prepare the KNOWN_CASH_FLOW by looking at the past values.
                eqKOLeg->makeKnownCashFlow(beginIdx, atStep, nbAssets, hasFuture());
                // re-set libor KNOWN_CASH_FLOW when it's already breached.  
                if (inst->hasFloater && isConditionMet)
                    koFltSV->makeKnownCashFlows(monitoringDates[metAtStep], isConditionMet);
                if (isConditionMet){
                    DateTime payDt = eqKOLeg->payDateOnKO(metAtStep);
                    knownRebate = CashFlowSP(new CashFlow(payDt,rebate->getLevel(metAtStep)));
                }
            }


            if (!pathGen->doingPast() || !hasFuture()) {
                // Compute a payoff, but only when we have a "complete" situation : either 
                // doingPast() and all is past, or !doingPast().

                // 3(new). calculate eq link coupon.  
                // eqKOLeg class has value as of today.
                double eqlinkPay = eqKOLeg->getValue(beginIdx, atStep, nbAssets, couponDfSV);
            
                // 4. aggregate into a basket
                basket() = assetBasket->aggregate();

                // 5. apply any modifier (call/put/etc) to basket
                overallOption->apply();

                // 6. form "expected payoff" combining overall option and rebate
                double rebateValue = isConditionMet? 
                    (couponDfSV->path()[eqKOLeg->payDateIndex(metAtStep)] * (rebate->getLevel(metAtStep) + redemption)) : 0.0;
                double payoff = koFactor * (basket() + redemption) * matDfSV->firstDF()
                                + (1.0 - koFactor) * rebateValue + eqlinkPay;
                if (inst->hasFloater)
                    payoff +=  isConditionMet ? koFltSV->getKOPV(monitoringDates[metAtStep]) 
                                                : liborLegSV->getPV(Today);
            
                // Done...
                prices.add(inst->notional * payoff); 
            }
        
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = refLevelGen->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);

            // Interpolate at the last date of the barrier level
            double barrierLevel = barrierHVTGen->getBarrierData()->
                getAssetBarrier(iAsset)->barrierLevels.back();
            double interpLevel  = fwdStarting ? 
                barrierLevel: barrierLevel * refLevelSV->refLevel(iAsset);

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));

            return reqarr;
    }

    /** For now am using this instead of implementing the IHandlePaymentEvents interface
        since I'd like to reuse the PAYMENT_DATES and KNOWN_CASHFLOWS features from
        the IMCProduct, but need to add BARRIER_LEVEL support. Should review. XXX */
    void recordExtraOutput(CControl*     control,
                           Results*      results,
                           const IMCPrices& prices) const {
        static const string method = "RainbowYenman::recordExtraOutput";
        try{
            CashFlowArraySP knownCfl = CashFlowArraySP(new CashFlowArray(0));

            // KNOWN CASH FLOW
            // not yed implemented after Maturity...
            OutputRequest* request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished())
            {   // add plain swap as known cash flow.  (Although, could be cancelled....) 
                //CashFlowArraySP rng = CashFlowArraySP(new CashFlowArray(eqKOLeg.getKnownCashFlows()));
                CashFlowArraySP rng = eqKOLeg->getKnownCashFlows();
                if (inst->hasFloater)
                    knownCfl = CashFlow::merge(rng,koFltSV->getKnownCashFlows());
                else 
                    knownCfl = rng;
                // when KO occurs, knownRebate is exiting and it's already include the redemption.
                // If all simulation date are in past, then need to add redemption at final.
                if (!!knownRebate)  
                    knownCfl->append(knownRebate);
                else if (!hasFuture()){  
                    rng->back().amount += redemption;
                }
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
                DateTimeArray listPayDates = eqKOLeg->getPayDates();
                DateTime::removeDuplicates(listPayDates,false);
                DateTimeArray payDateArray = DateTime::merge(listPayDates,cashPayDates);
                OutputRequestUtil::recordPaymentDates(control,results,&payDateArray);
            }

            // BARRIER_LEVEL is delegated to the state variable
            barrierHVTSV->recordBarrierLevels(control, results, getMultiFactors());

        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

};

///////////////////////////////////////////////////////////////////////////////////

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RainbowYenman::createProduct(const MonteCarlo* model) const {    
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    
    DateTimeArray eqObservDates = eqCpn->getObservDates();
    if (inAdvance)
        eqObservDates = DateTime::merge(eqObservDates, matDates);
    simSeries->addDates(eqObservDates);

    // get the barrier sampling date from barrier union
    BarrierSP realBar = BarrierSP(barSimHit.makeBarrierSchedule());
    realBar->validatePop2Object();
    ScheduleSP barSch = realBar->getBarrierSchedule();
    DateTimeArray barDates = barSch->getDateArray();

    // migrate bar dates & eq link coupon observation date
    DateTimeArray monitoringDates = DateTime::merge(eqCpn->getObservDates(), 
                                                    barDates);

    smartPtr<IMCPathConfig> pathConfig = model->getPathConfig();
    
    simSeries->addDates(monitoringDates);
        if(model->stateVarUsed()) {
            // State variables
            return new RainbowYenmanSVMC(this, simSeries);
        } else {
            // Otherwise, use old methodology
            return new RainbowYenmanMC(this, simSeries);
        }
}


CClassConstSP const RainbowYenman::TYPE = CClass::registerClassLoadMethod(
    "RainbowYenman", typeid(RainbowYenman), RainbowYenman::load);

// * for class loading (avoid having header file) */
bool RainbowYenmanLoad() {
    return (RainbowYenman::TYPE != 0);
}

CClassConstSP const RainbowYenman::BarSimHit::TYPE = CClass::registerClassLoadMethod(
    "RainbowYenman::BarSimHit", typeid(RainbowYenman::BarSimHit), RainbowYenman::BarSimHit::load);

DRLIB_END_NAMESPACE
