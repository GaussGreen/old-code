//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KOFixedLeg.cpp
//
//   Description   class for KO Fixed Leg 
//
//
//   $Log: KOFixedLeg.cpp,v $
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/KOFixedLeg.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"

DRLIB_BEGIN_NAMESPACE

//-----------------------------------//
// make class                        //
//-----------------------------------//
KOFixedLegMaker::KOFixedLegMaker():CObject(TYPE) {}; 
KOFixedLegMaker::KOFixedLegMaker(FixedLegSP fix, 
                                 KOStubRuleSP koRule)
                                 :CObject(TYPE),fixedLeg(fix), koRule(koRule) {}; 

void KOFixedLegMaker::validatePop2Object() {
    static const string method("KOFixedLeg::validatePop2Object");
    try {
        if (!!koRule){
            if (koRule->getKOStubRule() != "S" && koRule->getIsAccrueUpToSettle())
                throw ModelException(method, "isAccrueUpToSettle should not be true, if koStubRule is not 'S'.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

KOFixedLeg* KOFixedLegMaker::makeKOFixedLeg(const InstrumentSettlementSP settle) const{
    this->koRule->validatePop2Object(); // set payTimeRule;
    KOFixedLeg* koFixed = new KOFixedLeg(this, settle);
    return koFixed;
}

DateTimeArray KOFixedLegMaker::getPayDates() const
{
    DateTimeArray payDates = fixedLeg->PaymentDatesArray;
	return payDates;
}

/** Invoked when Class is 'loaded' */
void KOFixedLegMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Callable Deposit Call Schedule");
    REGISTER(KOFixedLegMaker, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultKOFixedLegMaker);
    FIELD(koRule, "ko rule for Fixed  Leg");
    FIELD_MAKE_OPTIONAL(koRule);
    FIELD(fixedLeg, "Fixed Leg");
    FIELD_MAKE_OPTIONAL(fixedLeg);
}
    
IObject* KOFixedLegMaker::defaultKOFixedLegMaker(){
    return new KOFixedLegMaker();
}


CClassConstSP const KOFixedLegMaker::TYPE = CClass::registerClassLoadMethod(
    "KOFixedLegMaker", typeid(KOFixedLegMaker), KOFixedLegMaker::load);

/////////////////////
// internal class ///
/////////////////////
KOFixedLeg::KOFixedLeg(const KOFixedLegMaker* maker,
                       const InstrumentSettlementSP settle):maker(copy(maker))
{
    // copy 
    fixedLeg = maker->fixedLeg;
    // validation of accru dates
    for (int i=1; i<fixedLeg->getSize(); i++){
        if (fixedLeg->AccrueStartDates[i] != fixedLeg->AccrueEndDates[i-1])
            throw ModelException("KOFixedLeg::KOFixedLeg", "AccrueStartDates should be previous AccrueEnd Date to use KOFixedLeg");
    }
    DateTimeArray accrual = fixedLeg->AccrueStartDates;
    accrual.push_back(fixedLeg->AccrueEndDates[fixedLeg->getSize()-1]);
    koSettle = KOSettleSP(maker->koRule->makeKOSettle(maker->koRule,
                fixedLeg->PaymentDatesArray,
                accrual,
                settle));
}

// constructor without maker.
KOFixedLeg::KOFixedLeg(const FixedLegSP fix,
                       const KOStubRuleSP koRule,
                       const InstrumentSettlementSP settle):fixedLeg(fix)
{
	DateTimeArray accrual = fixedLeg->AccrueStartDates;
    accrual.push_back(fixedLeg->AccrueEndDates.back());
    string test2 = koRule->getKOStubRule();

    koSettle = KOSettleSP(koRule->makeKOSettle(koRule,
                                    fixedLeg->PaymentDatesArray,
                                    accrual,
                                    settle));

    string test = koSettle->getKOStubRule();
}

// constructor without maker w/o settle
KOFixedLeg::KOFixedLeg(const FixedLegSP fix,
                       const KOStubRuleSP koRule):fixedLeg(fix)
{
	DateTimeArray accrual = fixedLeg->AccrueStartDates;
    accrual.push_back(fixedLeg->AccrueEndDates.back());
    koSettle = KOSettleSP(koRule->makeKOSettle(koRule,
                fixedLeg->PaymentDatesArray,
                accrual));
}

// return the critical dates.
DateTimeArray KOFixedLeg::getCritDates() const
{
    // to do : depend on KO stub rule and settlement rule, crit date should be different.
    DateTimeArray accstart = fixedLeg->AccrueStartDates;
    DateTimeArray payDates = fixedLeg->PaymentDatesArray;
    DateTimeArray critDates = DateTime::merge(accstart,payDates);
	return critDates;
}

// return AccrueStartDates.
DateTimeArray KOFixedLeg::getAccrueDates() const
{
    return fixedLeg->AccrueStartDates;
}

// return the KO value as of hitDate
double KOFixedLeg::getKOValue(const DateTime hitDate, const YieldCurve* yc){
    CashFlow cfl;
    if (getCashFlowOnKO(hitDate, koSettle->getKOStubRule(), &cfl))
        return cfl.amount * yc->pv(hitDate, cfl.date);                
    else 
        return 0.0;
}

// return non-KO (plain) value as of hitDate
double KOFixedLeg::getPV(const DateTime hitDate, const YieldCurve* yc){
    return fixedLeg->getPV(hitDate, yc);
}

// get pv (present = valueDate in disocunt) of cashflows which will be canceled KO on hitDate
// If Bond type (B), remove the full amount of current coupon, 
//                   whose accrue end is same day or in future as of hit date.
// if Swap type (S), remove accrued so far.
// otherwise, return total pv of leg value.
double KOFixedLeg::getKOPV(const DateTime&   hitDate, 
                       const YieldCurve* discount) {
    static const string method("FixedLeg::getKOPV");
    try {
        double pv = getPV(hitDate,discount);
        string koStubRule = koSettle->getKOStubRule();
        if (koStubRule == "B"){
            // substract already accrue started coupon value, but not yet paid.
            // such coupons are not cancelled.
            for(int i=0; i<fixedLeg->AccrueStartDates.size(); i++){
                if (fixedLeg->AccrueStartDates[i] < hitDate && fixedLeg->PaymentDatesArray[i] > hitDate){
                    pv -= fixedLeg->CouponAmounts[i] * discount->pv(hitDate, fixedLeg->PaymentDatesArray[i]);
                }
            }
        }
        else if (koStubRule == "S"){
            int k = Neighbour(hitDate, fixedLeg->AccrueEndDates, 0, fixedLeg->AccrueEndDates.size()-1, 1);
            if (k>=0){
                // substract accrued so far.
                pv -= fixedLeg->getAccrued(hitDate) * discount->pv(hitDate, fixedLeg->PaymentDatesArray[k]);
            }
        }
        return pv;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// return cash flow when KO occur on hitDate
CashFlow KOFixedLeg::getKOCashFlow(const DateTime hitDate) 
{
    CashFlow cfl;
    if (getCashFlowOnKO(hitDate, koSettle->getKOStubRule(), &cfl)){
        // success to get cfl!
    }
    else {
        cfl.date = hitDate;
        cfl.amount = 0.0;
    }
    return cfl;	    
}

CashFlowArrayConstSP KOFixedLeg::getCashFlows() 
{
	CashFlowArrayConstSP cf = fixedLeg->getCashFlowArray();
    return cf;
}

CashFlowArraySP KOFixedLeg::getKnownCashFlows(const DateTime hitDate) 
{
	knownCFLs = fixedLeg->getCashFlowArray(hitDate);
    return knownCFLs;
}


////// private functions /////////////

// return Fixed Leg after KO //
CashFlowArraySP KOFixedLeg::getKOCashFlowArray(const DateTime& hitDate, 
                                           const YieldCurve* discount,
                                           const string koStubRule,
                                           const bool isReverse){
    int i,k;
    CashFlowArraySP fixedflow = fixedLeg->getCashFlowArray(hitDate);
    if (koStubRule == "B"){
        // substract already accrue started coupon value, but not yet paid.
        // such coupons are not cancelled.
        for(int i=0; i<fixedflow->size();i++){
            k = Neighbour((*fixedflow)[i].date,fixedLeg->PaymentDatesArray,0,fixedLeg->PaymentDatesArray.size()-1,0);
            if(k<0)
                break;
            if(fixedLeg->AccrueStartDates[k] < hitDate){
                (*fixedflow)[i].amount = 0.0;                
            }            
        }
    }
    else if (koStubRule == "S"){//Not Tested
        k = Neighbour(hitDate, fixedLeg->AccrueEndDates, 0, fixedLeg->AccrueEndDates.size()-1, 1);
        if (k>=0){
            // substract accrued so far.  Assuming first fixed Leg.
            (*fixedflow)[0].amount -= fixedLeg->getAccrued(hitDate) * discount->pv(hitDate, fixedLeg->PaymentDatesArray[k]);
        }
    }    
    if (isReverse){
        for(i=0;i<(*fixedflow).size();i++){
            (*fixedflow)[i].amount *= -1.0;
        }
    }
    return fixedflow;
}

// return true if there is still cash flow even KO occurs on hitDate.
// Also, return cash flow, which are not cancelled by KO occurs on hitDate.  
// If the accrue end (in KOSettle) is past, then no cash flow returns, as same as "N" case.
bool KOFixedLeg::getCashFlowOnKO(const DateTime hitDate,
                                 const string koStubRule,
                                 CashFlow*  cfl){
    static const string method("FixedLeg::getCashFlowOnKO");
    try {
        int numSize = fixedLeg->getSize();
        if (hitDate >= fixedLeg->AccrueEndDates[numSize-1])
            return false;       // no coupon on this date
        else {
            // Usually, coupon are considered as "from Include / to Exclude". 
            // i.e. today = Accrue Date means it's now accrue start date and previous coupon are in past.
            // Thus, when KO occurs on Accrue End Date, koStub Rule do nothing.
            // When KO occurs on AccrueStartDate, "B" pays full cpn of next.

            DateTime when;            
            if (koSettle->getSettleDate(hitDate, when))
                cfl->date = when;
            else
                cfl->date = hitDate;
            when = koSettle->getAccrueUpToSettle() ? cfl->date : hitDate;

            int i, k;
            if (koStubRule == "S"){
                if (!fixedLeg->isDCC())
                    throw ModelException(method, "when you use koStubRule = S, you need turn on isDCC and specify dcc.");
                else if ( (k = Neighbour(when, fixedLeg->PaymentDatesArray, 0, numSize-1, 0)) >= 0 
                    && hitDate < fixedLeg->PaymentDatesArray[k]){ 
                    // be careful.  getAccrued return 0 if when = payment date.
                        DayCountConventionSP daycount = DayCountConventionSP(DayCountConventionFactory::make(fixedLeg->dcc));
                        cfl->amount = fixedLeg->CouponAmounts[k]*daycount->years(fixedLeg->AccrueStartDates[k], fixedLeg->AccrueEndDates[k]);        
                }
                else
                    cfl->amount = fixedLeg->getAccrued(when);                                          
            }
            else if (koStubRule == "B"){
                for(i=0;i<numSize;i++){
                    if (fixedLeg->AccrueStartDates[i]<hitDate && hitDate<=fixedLeg->AccrueEndDates[i]){
                        cfl->amount = fixedLeg->CouponAmounts[i];                                      
                        if (hitDate == fixedLeg->PaymentDatesArray[i])
                            cfl->amount = 0.0;  // already paid, so exclude from "KO cashflow".
                        if (fixedLeg->isDCC()){
                            DayCountConventionSP daycount = DayCountConventionSP(DayCountConventionFactory::make(fixedLeg->dcc));
                            cfl->amount *= daycount->years(fixedLeg->AccrueStartDates[i], fixedLeg->AccrueEndDates[i]);        
                        }
                    }
                }
            }
            else if (koStubRule =="N"){
                // no coupon 
                return false;
            }
            else{
                throw ModelException(method, " invalid koStubRule request.");
            }                
            return true;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE

