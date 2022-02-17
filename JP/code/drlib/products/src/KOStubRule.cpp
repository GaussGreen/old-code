//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EqLinkCashFlow.hpp
//
//   Description   Coupon stream linked to Equity generalized performance
//
//
//   $Log: EqLinkCashFlow.cpp,v $
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/KOStubRule.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/UtilFuncs.hpp"

DRLIB_BEGIN_NAMESPACE

KOStubRule::KOStubRule():CObject(TYPE),isAccrueUpToSettle(false){
    koStubRule = "N";
    payTimingRule = "AsSchedule";
}

KOStubRule::KOStubRule(const string koStubRule,
                       const bool isAccrueUpToSettle,
                       const string payTimingRule):CObject(TYPE),
                                                   koStubRule(koStubRule),
                                                   isAccrueUpToSettle(isAccrueUpToSettle),
                                                   payTimingRule(payTimingRule){
    validatePop2Object();
};

void KOStubRule::validatePop2Object(){
    static const string method("KOStubRule::validatePop2Object");
    try {
        if (payTimingRule == "AsSchedule")
            payTimeRule = AsSchedule;
        else if (payTimingRule == "AsInstSettle")
            payTimeRule = AsInstSettle;
        else
            throw ModelException(method, "payTiming Rule '" + payTimingRule + "' is unknown type");
        // not necessary, but not tested yet.
        if (koStubRule != "S" && isAccrueUpToSettle)
            throw ModelException(method, "isAccrueUpToSettle should not be true, if koStubRule is not 'S'.");
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

KOSettle* KOStubRule::makeKOSettle(KOStubRuleSP koRule,
                                    const DateTimeArray& payDates,
                                    const DateTimeArray& accrueDates){
    return new KOSettle(koRule, payDates, accrueDates);
}

KOSettle* KOStubRule::makeKOSettle(KOStubRuleSP koRule,
                                    const DateTimeArray& payDates,
                                    const DateTimeArray& accrueDates,
                                    InstrumentSettlementSP settle){
    return new KOSettle(koRule,settle,payDates,accrueDates);
}

string KOStubRule::getKOStubRule(){
    return koStubRule;
}

bool KOStubRule::getIsAccrueUpToSettle(){
	return isAccrueUpToSettle; 
}

string KOStubRule::getpayTimingRule(){
	return payTimingRule; 
}

        //
bool KOStubRule::isInstSettle(){
    return payTimeRule == AsInstSettle;
}
        
        /** Invoked when Class is 'loaded' */
void KOStubRule::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("CallableEquityKOSwap KO Rule");
    REGISTER(KOStubRule, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultKOStubRule);
    FIELD(koStubRule, "B = Bond (full pay) : N = None (No) : S = Swap (Accrue So Far)");
    FIELD_MAKE_OPTIONAL(koStubRule);
    FIELD(isAccrueUpToSettle, "accrue up to payment date for true.  Only active when koStubRule = S.");
    FIELD_MAKE_OPTIONAL(isAccrueUpToSettle);
    FIELD(payTimingRule, "payment timing rule");            
    FIELD_MAKE_OPTIONAL(payTimingRule);
}

IObject* KOStubRule::defaultKOStubRule(){
    return new KOStubRule();
}

CClassConstSP const KOStubRule::TYPE = CClass::registerClassLoadMethod(
    "KOStubRule", typeid(KOStubRule), KOStubRule::load);


// constructor  (is this good?  should take koStubRule, no?)
KOSettle::KOSettle(KOStubRuleSP            koRule,
                   const InstrumentSettlementSP  settle,
                   const DateTimeArray           payDates,
                   const DateTimeArray           accrueDates,
                   const bool                    isRollingSettle):koRule(koRule), settle(settle){

    // koRule class... 
    this->koRule->validatePop2Object();
    this->settle = settle;
    this->payDates = payDates;
    this->accrueDates = accrueDates;
    this->isRollingSettle = isRollingSettle;
    // validation
    if (accrueDates.size() != payDates.size()+1)
        throw ModelException("KOSettle","the size of payDates should be size of accrueDates -1!");
}

// constructor w/o settle
KOSettle::KOSettle(KOStubRuleSP            koRule,
                   const DateTimeArray           payDates,
                   const DateTimeArray           accrueDates):koRule(koRule){
    this->koRule->validatePop2Object(); // set payTimingRule
    this->payDates = payDates;
    this->accrueDates = accrueDates;
    this->isRollingSettle = this->koRule->isInstSettle();
    // validation
    if (accrueDates.size() != payDates.size()+1)
        throw ModelException("KOSettle","the size of payDates should be size of accrueDates -1!");
}

// constructor w/ settle
KOSettle::KOSettle(KOStubRuleSP            koRule,
                       const InstrumentSettlementSP  settle,
                       const DateTimeArray           payDates,
                       const DateTimeArray           accrueDates):koRule(koRule), settle(settle){
    this->koRule->validatePop2Object(); // set payTimingRule
    this->payDates = payDates;
    this->accrueDates = accrueDates;
    this->isRollingSettle = this->koRule->isInstSettle();
    // validation
    if (accrueDates.size() != payDates.size()+1)
        throw ModelException("KOSettle","the size of payDates should be size of accrueDates -1!");
}

bool KOSettle::isEmpty() const {
    return payDates.empty();
}

//////////////////////////////
// KO Settle Class          //
//////////////////////////////
// return the settleDate.  Those data should be set up by product class.
bool KOSettle::getSettleDate(const DateTime hitDate, DateTime& settleDate){
    static const string method = "KOSettle::settleDate";        
    try
    {
        if (hitDate < this->getFirstAccrueDate() || this->getLastAccrueDate() <= hitDate)
            return false;
        else if (isRollingSettle){
            if (CashSettleDate::TYPE->isInstance(settle.get())) {
                throw ModelException(method,"Rebate pay at hit does not allow instrument settlement of CashSettleDate!");
            }
            settleDate = settle->settles(hitDate,0);        
        }
        else{
            int nSize = payDates.size();
            // Basically, From Include / To Exclude.  i.e.  "On Date" means it's already next period.
            // Find out the current period, based on AccrueStart (not Accrue End to achive From Include/ To Exclude)
            int k = Neighbour(hitDate, accrueDates, 0, nSize-1, -1);
            if (k<0 || k >= nSize) 
                throw ModelException(method,"hitDates is not found in front of any accrueDates.");
            settleDate = payDates[k];
        }
        return true;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

bool KOSettle::getAccrueUpToSettle(){
    return koRule->getIsAccrueUpToSettle();
}

DateTimeArray KOSettle::getAccrueDates(){
    return accrueDates;
}

DateTime KOSettle::getLastAccrueDate(){
    return accrueDates[accrueDates.size()-1];
}

DateTime KOSettle::getFirstAccrueDate(){
    return accrueDates[0];
}

string KOSettle::getKOStubRule(){
    return koRule->koStubRule;
}

bool KOSettle::hasSettle(){
    bool hasSettle = (!!settle.get());
    return hasSettle;
}
        
DRLIB_END_NAMESPACE
