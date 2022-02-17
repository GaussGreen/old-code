//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EGKnockIn.cpp
//
//   Description : As CallableEquityKOSwap but with limited input for EGK trades:
//
//   Date        : May 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/SwapLegIntFace.hpp"
#include "edginc/Average.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/CallableEquityKOSwap.hpp"
#include "edginc/EGKBond.hpp"

DRLIB_BEGIN_NAMESPACE


/*****************************************************************************/

/// EGKnockIn product - as CKS but with limited interface
class EGKnockIn: public Generic1Factor, 
              virtual public FDModel::IIntoProduct,
              virtual public LastSensDate {

public:
    //---------------------------------------------//
    //  Class for Simple Equity Linked CashFlow   //
    //---------------------------------------------//
    class     MyEqLinkCashFlow: public CObject {
    public:
        static CClassConstSP const TYPE;  

        void validatePop2Object(){
            static const string method("MyEqLinkCashFlow::validatePop2Object");
            try {
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        };
        
        // data
        IDoubleArrayModifierMakerSP  eqLinkPerf;        // performace of each observation d
        DateTimeArray                observDates;      // eq link coupons observation Dates
        DateTimeArray                paymentDates;      // eq link coupons payment Dates.
        //int                          offset;            // offset for maturity of eqCpn.

        MyEqLinkCashFlow():CObject(TYPE){}; 
   
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("EGK Bond's Equity Linked CashFlows");
            REGISTER(EGKnockIn::MyEqLinkCashFlow, clazz);
            SUPERCLASS(CObject);
            FIELD(eqLinkPerf,       "performace of each observation date");
            FIELD(observDates,        "each eq link coupons observation dates");
            FIELD(paymentDates,       "coupon payment dates");
            FIELD_MAKE_OPTIONAL(paymentDates);
            EMPTY_SHELL_METHOD(defaultMyEqLinkCashFlow);
        };
        
        static IObject* defaultMyEqLinkCashFlow(){
            return new EGKnockIn::MyEqLinkCashFlow();
        };
    };



private: 
    /// fields ////////
    EGKnockIn::MyEqLinkCashFlow eqCpnSimple;
    KnockInPerfMakerSP  knockInPerf;         
    double            redemption;    
    LiborLegSP  floater;
    EGKBond::SimpleTrigger  trigger;
    EGKBond::AverageIn     averageIn;

    CashFlowArray samples;
    bool          isSkipCheckSamples;

#ifdef  TREE_THETA_CAP
    // for theta smoothing
    double  thetaSmoothThrehold;    
    bool    isPositiveThetaSmooth;
#endif

    // for internal use  
    int           numOfAvgInSV; // $unregistered

    CallableEquityKOSwapSP       originalCKS;    
    
    InstrumentSettlementSP myFinalSettle; // $unregistered

public:
    static CClassConstSP const TYPE;
    
    // validation
    void validatePop2Object(){
    static const string method = "EGKnockIn::validatePop2Object";
    try{
        // no strong reason to comment out below.  But IMS can prohibit it, so leave the library can take the flag = true.
//        if (fwdStarting)
//            throw (method, "EGKnockIn doesn't use fwdStarting flag.  Need to be false.  Use Average-In for Fwd Starting");
    } catch (exception& e){
        throw ModelException(e, method);
        }
    }

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
    static const string method = "EGKnockIn::GetMarket";
    try{
        int i, j;

        // validation
        if (eqCpnSimple.observDates.size() <=1)
            throw ModelException(method, "Equity Linked coupon size should be >=2");

        DateTimeArray sampleDates = CashFlow::dates(samples);
        DateTimeArray barDates = trigger.barrier->getDates();
        if (!DateTime::isSubset(sampleDates, eqCpnSimple.observDates) && !isSkipCheckSamples)
            throw ModelException(method, "samples schedule should include all equity coupon observetion dates.  Check observation timing, too.");
        if (!DateTime::isSubset(sampleDates, barDates) && !isSkipCheckSamples)
            throw ModelException(method, "samples schedule should include all trigger observetion dates.  Check observation timing, too.");
        // check barrier and equity observation date are same date.
        for (i=j=0; i<barDates.size() && j<eqCpnSimple.observDates.size() ; i++){
            if (barDates[i] < eqCpnSimple.observDates[j]){
                // nothing.  just check next bar date.
            }
            else if (barDates[i] != eqCpnSimple.observDates[j]){
                throw ModelException(method, "equity obaservation date should be 1) same to trigger observetion dates, or 2) any dates later than the last observation date.");
            }
            else
                j++;
        }
        // validation of KnockInPerfMaker
        if (averageIn.numAvgDates() <=0)
            throw ModelException(method, "Empty in average-in schedule!  Needs a date + sampled level (0 for future) at least.");
        DateTimeArray   startDates = CashFlow::dates(averageIn.avgIn);
        DateTime        endOfStart = startDates[startDates.size()-1];
        bool            isCont;
        ScheduleSP      knockInBar = knockInPerf->getBarrier(false,isCont);
        DateTimeArray   kiBarDates = knockInBar->getDates();
        DateTimeArray   matDates = knockInPerf->getMatDates();
        DateTime        lastEqSample = eqCpnSimple.observDates[eqCpnSimple.observDates.size()-1];
        if (kiBarDates[0] < endOfStart){
            throw ModelException(method, "First Date of knock-in monitor["
                                        + kiBarDates[0].toString() +
                                        "] should later \n than or equal to average start date["
                                        + endOfStart.toString() + "].");
        }
        if (kiBarDates[kiBarDates.size()-1] > lastEqSample){
            throw ModelException(method, "Last Date of knock-in monitor["
                                        + kiBarDates[kiBarDates.size()-1].toString() +
                                        "] should be earlier \n than  or equal to the last of eqCpn's observe date["
                                        + lastEqSample.toString() + "].");
        }
        if (matDates.size()!=1)
            throw ModelException(method, "matDates's size in knockInPerf should be 1.");
        // the followig is not necessary, but lazy to test....
        if (matDates[0] != lastEqSample)
            throw ModelException(method, "matDates [" + matDates[0].toString()
                                        + "] should be same to the last of eqCpn's observation["
                                        + lastEqSample.toString() + "].");


        if (!DateTime::isSubset(sampleDates, barDates) && !isSkipCheckSamples)
            throw ModelException(method, "samples schedule should include all trigger observetion dates.  Check observation timing, too.");

        //--------------- Setting up CKS ---------------//
        // dummy null objects
        KOFixedLegMaker nullKOFixed;
        CallableEquityKOSwap::Trigger nullTrigger;
        CallableEquityKOSwap::CallSchedule nullCallSched;

        // koStubRule = "Bond" type
        KOStubRuleSP koRule = KOStubRuleSP(new KOStubRule("B",false,"AsSchedule"));     // this settlemetn timing rule is not used, yet.
        
        ///////////////////////////////////////////
        //making EqCpn by adding koRule.
        //HolidayConstSP hol = AssetUtil::getHoliday(asset.get());
        //DateTimeArray observDates;
        //for (i=0;i<eqCpnSimple.paymentDates.size();i++)
        //    observDates[i] = hol->addBusinessDays(eqCpnSimple.paymentDates[i],-eqCpnSimple.offset);
        DateTime initialAccrueDate = eqCpnSimple.observDates[0].rollDate(-1);   // accrue start just one day before.

        
        EqLinkCashFlowSP eqCpn;
        if (eqCpnSimple.paymentDates.size()>0){
            eqCpn = EqLinkCashFlowSP(new EqLinkCashFlow(eqCpnSimple.eqLinkPerf,
                                                       eqCpnSimple.observDates,
                                                       eqCpnSimple.paymentDates));
        }else{
            eqCpn = EqLinkCashFlowSP(new EqLinkCashFlow(eqCpnSimple.eqLinkPerf,
                                                       eqCpnSimple.observDates));
        }

        EqCpnKOMakerSP koELC = EqCpnKOMakerSP(new EqCpnKOMaker(eqCpn,initialAccrueDate,koRule));

        ///////////////////////////////////////////
        //making koLiborLeg using above (same as eqCpn) koRule.
        bool hasFloater = (!!floater && floater->getSize() > 0);
        KOLiborLegMaker koLiborLeg;
        if (hasFloater)
            koLiborLeg = KOLiborLegMaker(floater, koRule);

        ///////////////////////////////////////////
        //making rebate.  
        // When KO, return redemption.  
        // So barrier schedule should not have any KO date 
        // after the last coupon accrue started!! (in LiborLeg and EqCpn)
        DateTime lastAccrueStart = eqCpnSimple.observDates[eqCpnSimple.observDates.size()-2];
        if (hasFloater){
            if (lastAccrueStart > floater->AccrualDates[floater->getSize()-1])
            lastAccrueStart = floater->AccrualDates[floater->getSize()-1];
        }

        DoubleArray rebValues = trigger.barrier->getValues();
        for (i=0; i<barDates.size(); i++){
            if (barDates[i] > lastAccrueStart)
                throw ModelException(method, "All barrier date should be earlier or eqaul than final accrue start date [" +
                                                lastAccrueStart.toString() + "].");
            rebValues[i] = redemption;
            if (trigger.bonusCoupons.size() > 0){
                rebValues[i] += trigger.bonusCoupons[i];
            }
        }

        ScheduleSP rebate = ScheduleSP(new Schedule(barDates, rebValues, trigger.barrier->getInterp()));
                                                    
        ///////////////////////////////////////////
        //making trigger
        DateTimeArray payDates;
        if (eqCpnSimple.paymentDates.size()>0)
            payDates = eqCpnSimple.paymentDates;
        else
            payDates = DateTimeArray(0);
        CallableEquityKOSwap::Trigger upTrigger(true,    // isUp
                                                trigger.barrier,
                                                trigger.ecoBarrier, 
                                                false,          //IntraDayMonitor
                                                true,           // hasRebate
                                                rebate,
                                                trigger.smoothingType,
                                                trigger.peakDelta,
                                                payDates);
      
        ///////////////////////////////////////////
        // making finalPerf
        Generic1Factor::GetMarket(model, market);   // Need to get market for finalPerf.
        myFinalSettle = instSettle;
        if (eqCpnSimple.paymentDates.size()>0){
            DateTime fpSettleDate = eqCpnSimple.paymentDates[eqCpnSimple.paymentDates.size()-1];            
            DateTime matDate = matDates[matDates.size()-1];
            if (fpSettleDate < matDate){
                throw ModelException(method, "settle date for FinalPerf[" 
                                            + fpSettleDate.toString() +
                                             "] should be later or equal than last of the matDates["
                                            + matDate.toString() + "] .");
            }
            // don't use instrument's instrument rule.
            myFinalSettle = smartPtr<CashSettleDate> (new CashSettleDate(fpSettleDate));
        }


        // Create a TargetRedemptionNote to make use of common validation and for pricing
        originalCKS = CallableEquityKOSwapSP(
                       new     CallableEquityKOSwap(
                                    // Generic1Factor
                                    valueDate, oneContract, notional, initialSpot, 
                                    ccyTreatment,instSettle, premiumSettle, asset, discount,                                    
                                    // for CKS
                                    true,                           //hasEqLeg,
                                    koELC,
                                    knockInPerf,                 // knockInPerf,         
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
                                    averageIn.avgIn,
                                    averageIn.useIniSpotAsRefLvl,
                                    samples,
                                    numOfAvgInSV,
                                    myFinalSettle
#ifdef  TREE_THETA_CAP
                                    ,
                                    thetaSmoothThrehold,
                                    isPositiveThetaSmooth
#endif
                                    ));


        originalCKS->validatePop2Object();

        // parent
        //Generic1Factor::GetMarket(model, market);
        // There is no EGKnockIn specific market data
        // All processing will be done on the originalCKS
        originalCKS->GetMarket(model, market);

    } catch (exception& e){
        throw ModelException(e, method);
        }
    }

    // Pass-through to internal CallableEquityKOSwap.
    // The EGK Bond itself has not acquired market data.
    void Validate(){
    static const string method = "EGKnockIn::Validate";
    try{
        if (!originalCKS.get())
            throw ModelException(method, "missing the object of CallableEquityKOSwap!!");
        originalCKS->Validate();
    } catch (exception& e){
        throw ModelException(e, method);
        }
    }

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const
    {
        return originalCKS->createProduct(model);
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const{
        DateTimeArray dts = knockInPerf->getMatDates();
        DateTime lastDate = dts[dts.size()-1]; 
        DateTime end = instSettle->settles(lastDate, asset.get());
        return end;    
    }

    bool priceDeadInstrument(CControl* control, CResults* results) const{
        return originalCKS->priceDeadInstrument(control, results);
    }
    
private: 
    friend class EGKnockInHelper;
    friend class EGKnockInProd;

    EGKnockIn():Generic1Factor(TYPE), numOfAvgInSV(1), isSkipCheckSamples(false)
#ifdef  TREE_THETA_CAP
                ,thetaSmoothThrehold(0.0),isPositiveThetaSmooth(false)
#endif
                {};// for reflection

    EGKnockIn(const EGKnockIn& rhs); // not implemented
    EGKnockIn& operator=(const EGKnockIn& rhs); // not implemented    
};

class EGKnockInHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("EGKnockIn");
        REGISTER(EGKnockIn, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        FIELD(knockInPerf, "knock-in option at maturity");
        FIELD(redemption, "final redemption amount");
        FIELD(eqCpnSimple, "Equity Linked Cash Flow");
        FIELD(trigger, "Trigger scheduel and levels");        
        FIELD(floater,"Libor Leg");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(averageIn,    "average-in");
        FIELD(isSkipCheckSamples,    "true for checking past samples.");
        FIELD_MAKE_OPTIONAL(isSkipCheckSamples);
        FIELD(samples,  "asset level sampling schedule and results");
        FIELD_MAKE_OPTIONAL(samples);
#ifdef  TREE_THETA_CAP
        FIELD(thetaSmoothThrehold, "The threshold for theta smoothing.");
        FIELD_MAKE_OPTIONAL(thetaSmoothThrehold);
        FIELD(isPositiveThetaSmooth, "true for smoothing so as to FV positively increase.");
        FIELD_MAKE_OPTIONAL(isPositiveThetaSmooth);
#endif
        FIELD(originalCKS, "Internal CallableEquityKOSwap");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(originalCKS);  
        EMPTY_SHELL_METHOD(defaultEGKnockIn);
    }
    static IObject* defaultEGKnockIn(){
        return new EGKnockIn();
    }
};

CClassConstSP const EGKnockIn::TYPE = CClass::registerClassLoadMethod(
    "EGKnockIn", typeid(EGKnockIn), EGKnockInHelper::load);

CClassConstSP const EGKnockIn::MyEqLinkCashFlow::TYPE = CClass::registerClassLoadMethod(
    "EGKnockIn::MyEqLinkCashFlow", typeid(EGKnockIn::MyEqLinkCashFlow), EGKnockIn::MyEqLinkCashFlow::load);

// * for class loading (avoid having header file) */
bool EGKnockInLoad() {
    return (EGKnockIn::TYPE != 0);
}

DRLIB_END_NAMESPACE

