//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EGKBond.cpp
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

// validation
void EGKBond::validatePop2Object(){
    static const string method = "EGKBond::validatePop2Object";
    try{
    //        originalCKS->validatePop2Object();
    //        originalCKS->Validate();

    } catch (exception& e){
        throw ModelException(e, method);
        }
}

/** Get the asset and discount market data */
void EGKBond::GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
    static const string method = "EGKBond::GetMarket";
    try{
        int i, j;

        // validation
        if (fixedLeg->getSize() <=1)
            throw ModelException(method, "fixed Leg size should be >=2");
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
            
        if (!DateTime::isSubset(sampleDates, barDates) && !isSkipCheckSamples)
            throw ModelException(method, "samples schedule should include all trigger observetion dates.  Check observation timing, too.");

        //--------------- Setting up CKS ---------------//
        // dummy null objects
        KOLiborLegMaker nullKOLibor;
        CallableEquityKOSwap::Trigger nullTrigger;
        CallableEquityKOSwap::CallSchedule nullCallSched;

        // koStubRule = "Bond" type
        KOStubRuleSP koRule = KOStubRuleSP(new KOStubRule("B",false,"AsSchedule"));
    
        //making EqCpn by adding koRule.
        //HolidayConstSP hol = AssetUtil::getHoliday(asset.get());
        //DateTimeArray observDates;
        //for (i=0;i<eqCpnSimple.paymentDates.size();i++)
        //    observDates[i] = hol->addBusinessDays(eqCpnSimple.paymentDates[i],-eqCpnSimple.offset);
        DateTime initialAccrueDate = eqCpnSimple.observDates[0].rollDate(-1);   // accrue start just one day before.

        ///////////////////////////////////////////////
        // making payment date for EqLinkCpn...
        j=0;
        DateTimeArray payDates;
        for (i=0; i<eqCpnSimple.observDates.size(); i++){
            while (eqCpnSimple.observDates[i] > fixedLeg->PaymentDatesArray[j] && j<fixedLeg->getSize()){
                j++;
            }
            if (j==fixedLeg->getSize()){
                throw ModelException(method, "fixedLeg payment date should have any date which is later than eqObservDate[" +
                                              eqCpnSimple.observDates[i].toString() +
                                              "].");
            }
            payDates.push_back(fixedLeg->PaymentDatesArray[j]);
        }

        EqLinkCashFlowSP eqCpn;
        if (payDates.size()>0)
            eqCpn = EqLinkCashFlowSP(new EqLinkCashFlow(eqCpnSimple.eqLinkPerf,
                                                       eqCpnSimple.observDates,
                                                       payDates));
        else
            eqCpn = EqLinkCashFlowSP(new EqLinkCashFlow(eqCpnSimple.eqLinkPerf,
                                                       eqCpnSimple.observDates));

        EqCpnKOMakerSP koELC = EqCpnKOMakerSP(new EqCpnKOMaker(eqCpn,initialAccrueDate,koRule));

        //making koFixedLeg using above (same as eqCpn) koRule.
        KOFixedLegMaker koFixedLeg(fixedLeg, koRule);

        ///////////////////////////////////////////////
        //making rebate.  
        // When KO, return redemption.  
        // So barrier schedule should not have any KO date 
        // after the last coupon accrue started!! (in FixedLeg and EqCpn)
        DateTime lastAccrueStart = fixedLeg->AccrueStartDates[fixedLeg->getSize()-1];
        if (lastAccrueStart > eqCpnSimple.observDates[eqCpnSimple.observDates.size()-2])
            lastAccrueStart = eqCpnSimple.observDates[eqCpnSimple.observDates.size()-2];

        DoubleArray rebValues = trigger.barrier->getValues();
        for (i=0; i<barDates.size(); i++){
            if (barDates[i] > lastAccrueStart)
                throw ModelException(method, "All barrier date should be earlier than final accrue start date [" +
                                                lastAccrueStart.toString() + "].");
            rebValues[i] = redemption;
            if (trigger.bonusCoupons.size() > 0){
                rebValues[i] += trigger.bonusCoupons[i];
            }
        }

        ScheduleSP rebate = ScheduleSP(new Schedule(barDates, rebValues, trigger.barrier->getInterp()));
        
        ///////////////////////////////////////////////                                                
        //making trigger
        CallableEquityKOSwap::Trigger upTrigger(true,    // isUp
                                                trigger.barrier,
                                                trigger.ecoBarrier, 
                                                false,          //IntraDayMonitor
                                                true,           // hasRebate
                                                rebate,
                                                trigger.smoothingType,
                                                trigger.peakDelta,
                                                fixedLeg->PaymentDatesArray);
        
        ///////////////////////////////////////////////
        // making finalPerf
        Generic1Factor::GetMarket(model, market);   // Need to get market for finalPerf.
        finalPerfMaker = AveragePerfMakerSP(new AveragePerfMaker(avgOutPerf));
        DateTime fpSettleDate = fixedLeg->PaymentDatesArray[fixedLeg->getSize()-1];
        DateTimeArray matDates = finalPerfMaker->getMatDates();
        DateTime matDate = matDates[matDates.size()-1];
        if (fpSettleDate < matDate){
            throw ModelException(method, "settle date for FinalPerf[" 
                                        + fpSettleDate.toString() +
                                         "] should be later or equal than last of the matDates["
                                        + matDate.toString() + "] .");
        }
        // don't use instrument's instrument rule.
        myFinalSettle = smartPtr<CashSettleDate> (new CashSettleDate(fpSettleDate));

        ///////////////////////////////////////////////
        // Create a CallableEquityKOSwap to make use of common validation and for pricing
        originalCKS = CallableEquityKOSwapSP(
                       new     CallableEquityKOSwap(
                                    // Generic1Factor
                                    valueDate, oneContract, notional, initialSpot, 
                                    ccyTreatment,instSettle, premiumSettle, asset, discount,                                    
                                    // for CKS
                                    true,                           //hasEqLeg,
                                    koELC,
                                    finalPerfMaker,                 // avgOutPerf,         
                                    redemption,    
                                    true,                           // hasFixedLeg,
                                    false,                          // hasFloater,
                                    koFixedLeg,                     // 
                                    nullKOLibor,                           // koFloater
                                    true,                           // hasUpTrigger,
                                    false,                          // hasLowTrigger,
                                    upTrigger,
                                    nullTrigger,                           // lowTrigger
                                    false,                          // hasCallSchedule,
                                    nullCallSched,                        // callSchedule,    
                                    isAvgIn,
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
        // There is no EGKBond specific market data
        // All processing will be done on the originalCKS
        originalCKS->GetMarket(model, market);

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

// Pass-through to internal CallableEquityKOSwap.
// The EGK Bond itself has not acquired market data.
void EGKBond::Validate(){
    originalCKS->Validate();
}


/** create a fd payoff tree - new for all fd/tree state variable interface */
FDProductSP EGKBond::createProduct(FDModel* model) const
{
    return originalCKS->createProduct(model);
}

/** when to stop tweaking */
DateTime EGKBond::endDate(const Sensitivity* sensControl) const{
    DateTime lastDate = avgOutPerf->getLastDate();
    DateTime end = instSettle->settles(lastDate, asset.get());
    return end;    
}

bool EGKBond::priceDeadInstrument(CControl* control, CResults* results) const{
    return originalCKS->priceDeadInstrument(control, results);
}


class EGKBondHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("EGKBond");
        REGISTER(EGKBond, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        FIELD(avgOutPerf, "option at maturity, AvgOutPerfMaker");
        FIELD(redemption, "final redemption amount");
        FIELD(eqCpnSimple, "Equity Linked Cash Flow");
        FIELD(fixedLeg,  "fixed cash flow");
        FIELD(trigger, "Trigger scheduel and levels");        
        FIELD(isAvgIn,  "is average-in?");
        FIELD(averageIn,    "average-in");
        FIELD_MAKE_OPTIONAL(averageIn);
        FIELD(isSkipCheckSamples,    "true for checking past samples.");
        FIELD_MAKE_OPTIONAL(isSkipCheckSamples);
        FIELD(samples,  "asset level sampling schedule and results");
        FIELD_MAKE_OPTIONAL(samples);
//        FIELD(numOfAvgInSV, "number of state variable for Average-In.");
//        FIELD_MAKE_OPTIONAL(numOfAvgInSV);
#ifdef  TREE_THETA_CAP
        FIELD(thetaSmoothThrehold, "The threshold for theta smoothing.");
        FIELD_MAKE_OPTIONAL(thetaSmoothThrehold);
        FIELD(isPositiveThetaSmooth, "true for smoothing so as to FV positively increase.");
        FIELD_MAKE_OPTIONAL(isPositiveThetaSmooth);
#endif
        FIELD(originalCKS, "Internal CallableEquityKOSwap");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(originalCKS);  
        FIELD(finalPerfMaker, "just internal usage, but need to hold during the tweaks");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(finalPerfMaker);
        EMPTY_SHELL_METHOD(defaultEGKBond);
    }
    static IObject* defaultEGKBond(){
        return new EGKBond();
    }
};

CClassConstSP const EGKBond::TYPE = CClass::registerClassLoadMethod(
    "EGKBond", typeid(EGKBond), EGKBondHelper::load);

CClassConstSP const EGKBond::SimpleEqLinkCashFlow::TYPE = CClass::registerClassLoadMethod(
    "EGKBond::SimpleEqLinkCashFlow", typeid(EGKBond::SimpleEqLinkCashFlow), EGKBond::SimpleEqLinkCashFlow::load);

CClassConstSP const EGKBond::SimpleTrigger::TYPE = CClass::registerClassLoadMethod(
    "EGKBond::SimpleTrigger", typeid(EGKBond::SimpleTrigger), EGKBond::SimpleTrigger::load);

CClassConstSP const EGKBond::AverageIn::TYPE = CClass::registerClassLoadMethod(
    "EGKBond::AverageIn", typeid(EGKBond::AverageIn), EGKBond::AverageIn::load);

// * for class loading (avoid having header file) */
bool EGKBondLoad() {
    return (EGKBond::TYPE != 0);
}

DRLIB_END_NAMESPACE

