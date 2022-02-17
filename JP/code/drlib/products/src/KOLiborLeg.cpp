//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KOLiborLeg.cpp
//
//   Description   class for KO Libor & KO Fixed Leg 
//
//
//   $Log: KOLiborLeg.cpp,v $
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/KOLiborLeg.hpp"

DRLIB_BEGIN_NAMESPACE
   
//-----------------------------------//
// make class                        //
//-----------------------------------//
KOLiborLegMaker::KOLiborLegMaker():CObject(TYPE) {}; 
KOLiborLegMaker::KOLiborLegMaker(const LiborLegSP floater, 
                                 const KOStubRuleSP koRule)
                                 :CObject(TYPE),floater(floater), koRule(koRule){}; 

void KOLiborLegMaker::validatePop2Object() {
    static const string method("KOLiborLegMaker::validatePop2Object");
    try {
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

KOLiborLeg* KOLiborLegMaker::makeKOLiborLeg(const InstrumentSettlementSP settle,
                                            const DateTime valDate) const{
    KOLiborLeg* koLibor = new KOLiborLeg(this,settle,valDate);
    return koLibor;
}

DateTimeArray KOLiborLegMaker::getPayDates() const
{
	return floater->PayDates;
}

void KOLiborLegMaker::getMarket(const IModel*     model,
                                const MarketData* market,
                                const YieldCurveWrapper discount)
{
    floater->setCouponCurve(discount); // in case floater gets coupon curve from inst
    floater->getMarket(model, market); 
}

// tweak theata and feed fixing level....
void KOLiborLegMaker::setFixingforThetaShift(const DateTime& valueDate, 
                                              const YieldCurve* discount,
                                              const DateTime& rollDate)
{
    floater->setFixingforThetaShift(valueDate,
                                    discount,
                                    rollDate);	
}

/** Invoked when Class is 'loaded' */
void KOLiborLegMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Callable Deposit Call Schedule");
    REGISTER(KOLiborLegMaker, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultKOLiborLegMaker);
    FIELD(koRule, "ko rule for Libor Leg");
    FIELD_MAKE_OPTIONAL(koRule);
    FIELD(floater, "Libor Leg");
    FIELD_MAKE_OPTIONAL(floater);
}

IObject* KOLiborLegMaker::defaultKOLiborLegMaker(){
    return new KOLiborLegMaker();
}

CClassConstSP const KOLiborLegMaker::TYPE = CClass::registerClassLoadMethod(
    "KOLiborLegMaker", typeid(KOLiborLegMaker), KOLiborLegMaker::load);

/////////////////////
// internal class ///
/////////////////////

// all constructor allow to useCcyBasis.  
KOLiborLeg::KOLiborLeg(const KOLiborLegMaker* maker,
                       const InstrumentSettlementSP settle,
                       const DateTime valDate): valDate(valDate),floater(maker->floater)
{
    floater->validatePop2Object();
    // NB : setCouponCurve is already called in KOLiborLegMaker class.
    maker->koRule->validatePop2Object(); // set payTimeRule;
    koSettle.reset(maker->koRule->makeKOSettle(maker->koRule,
                                               floater->PayDates,
                                               floater->AccrualDates,
                                               settle));

}

KOLiborLeg::KOLiborLeg(const LiborLegSP flt,
                       const KOStubRuleSP koRule,
                       const InstrumentSettlementSP settle,
                       const DateTime valDate,
                       const YieldCurveWrapper discount):floater(flt),valDate(valDate)
{
    floater->validatePop2Object();
    floater->setCouponCurve(discount);
    koSettle.reset(koRule->makeKOSettle(koRule,
                                        floater->PayDates,
                                        floater->AccrualDates,
                                        settle));
}

// constructor, different settle schedule from floater's schedule.
KOLiborLeg::KOLiborLeg(LiborLegSP flt,
                       KOStubRuleSP koRule,
                       InstrumentSettlementSP settle,
                       const DateTimeArray& payDates,
                       const DateTimeArray& accrueDates,
                       DateTime valDate,
                       const YieldCurve* discount):floater(flt),valDate(valDate)
{
    YieldCurveWrapper disc(discount->getName());
    floater->validatePop2Object();
    floater->setCouponCurve(disc);
    koSettle.reset(koRule->makeKOSettle(koRule,
                                        payDates,
                                        accrueDates,
                                        settle));
}

// no settle version....
KOLiborLeg::KOLiborLeg(const LiborLegSP flt,
                       const KOStubRuleSP koRule,
                       const DateTime valDate,
                       const YieldCurve* discount):floater(flt),valDate(valDate)
{
    YieldCurveWrapper disc(discount->getName());
    floater->validatePop2Object();
    floater->setCouponCurve(disc);
    koSettle.reset(koRule->makeKOSettle(koRule,
                                        floater->PayDates,
                                        floater->AccrualDates));
}

// no settle version....
KOLiborLeg::KOLiborLeg(const LiborLegSP flt,
                       const KOSettleSP koSettle,
                       const DateTime valDate,
                       const YieldCurve* discount):
floater(flt),valDate(valDate), koSettle(koSettle){
    YieldCurveWrapper disc(discount->getName());
    floater->validatePop2Object();
    floater->setCouponCurve(disc);
}

// return the critical dates.
DateTimeArray KOLiborLeg::getCritDates() const
{
    // to do : depend on KO stub rule and settlement rule, crit date should be different.
    DateTimeArray accstart = floater->AccrualDates;
    DateTimeArray payDates = floater->PayDates;
    DateTimeArray critDates = DateTime::merge(accstart,payDates);
	return critDates;
}

// return AccrueDates.
DateTimeArray KOLiborLeg::getAccrueDates() const
{
    return floater->AccrualDates;
}

// return simple cash flow array
CashFlowArrayConstSP KOLiborLeg::getCashFlows(const DateTime valDate, const YieldCurve* yc)
{
    return floater->getCashFlowArray(valDate, yc);	// 
}

double KOLiborLeg::getKOValue(const DateTime hitDate, const YieldCurve* yc){
    double pvValue;
    CashFlowArray cfl;
    if (getCashFlowOnKO(hitDate, koSettle->getKOStubRule(), &pvValue, &cfl))
        return pvValue;
    else 
        return 0.0;
}

// return non-KO (plain) value as of hitDate
double KOLiborLeg::getPV(const DateTime hitDate,const YieldCurve* yc)
{
    DateTime last = floater->PayDates[floater->getSize()-1];
    double value = floater->getPV(valDate, hitDate, last, yc);
    return value;
}

// is this necessary???
double KOLiborLeg::getKOPV(const DateTime&   hitDate)
{
    return getKOPV(valDate, hitDate,  floater->couponCurve.get(), koSettle->getKOStubRule());
}

CashFlowArraySP KOLiborLeg::makeKnownCashFlows(const DateTime   hitDate,            //simulated path
                                              const bool       isHitAlready){        //is Hit in past                                             
    static const string method = "KOLiborLeg::makeKnownCashFlows";
    try{
       
        //initialization.
        knownCFL = CashFlowArraySP(new CashFlowArray(0));

        // in case of already KO.
        if (isHitAlready){
            if (koSettle->isEmpty())
                knownCFL = getKOCashFlowArray(hitDate,floater->couponCurve.get(),koSettle->getKOStubRule(),true); // ????
            else{
                // change the payment date & amount
                CashFlowArray cpn;
                double pvValue;
                if(getCashFlowOnKO(hitDate, koSettle->getKOStubRule(), &pvValue, &cpn)){
                    for (int i=0; i<cpn.size(); i++)
                            knownCFL->push_back(cpn[i]);
                }
            }
        }
        else{
            floater->makeKnownCashFlow();
            knownCFL = floater->getKnownCashFlows();
        }

        return knownCFL;
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
};

CashFlowArraySP KOLiborLeg::getKnownCashFlows(){
    return knownCFL;
}

////////////////////
// private class ///
////////////////////

// Input is KO date time array.  
// Output is the PV values (as of value Date) at correspoinding time point, when the KO occurs.  
// Need to call getMarket & setCouponCurve before.
DoubleArray KOLiborLeg::getPVsAlongKOTimeStep(const DateTimeArray&    koDates,
                                            const double            scaling){
    static const string method("KOLiborLeg::getPVsAlongKOTimeStep");
    try {
        int nbSteps = koDates.size();
        DoubleArray pvValues = DoubleArray(nbSteps, 0.0);  //initiazlize as 0.0
        // check the valueDate has been set already.
        floater->checkMarket();

        // calculate plain cash flow.
        // Then subtract KO Value, which is the amount to be cancelled by KO on hitDate.
        CashFlowArray cpn;   // dummy, not used
        double plainValue = getPV(valDate, floater->couponCurve.get());        
        DateTime when, settleDate;
        DateTime lastAccrue = koSettle->getLastAccrueDate();
        DateTime firstAccrue = koSettle->getFirstAccrueDate();
        double df;
        for (int i = 0; i<nbSteps; i++){
            pvValues[i] = plainValue;
            if (koDates[i] < firstAccrue)//nothing happen as the libor is not started.
                pvValues[i] = 0;
            else if(koSettle->getSettleDate(koDates[i], settleDate)){            
                if (settleDate > valDate){
                    df =  floater->couponCurve->pv(valDate,koDates[i]);
                    pvValues[i] -= df * getKOPV(valDate, koDates[i], floater->couponCurve.get(), koSettle->getKOStubRule());
                } else {
                    pvValues[i] = 0.0;      //already paid everything.
                }
            }
            pvValues[i] *= scaling;
        }
        return pvValues;

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return libor leg after KO //
CashFlowArraySP KOLiborLeg::getKOCashFlowArray(const DateTime& hitDate, 
                                           const YieldCurve* discount,
                                           const string koStubRule,
                                           const bool isReverse){
    int i,k;
    CashFlowArraySP koLibor = floater->getCashFlowArray(hitDate, hitDate, discount);
    if (koSettle->getKOStubRule() == "B"){
        // substract already accrue started coupon value, but not yet paid.
        // such coupons are not cancelled.
        for(int i=0; i<koLibor->size();i++){
            k = Neighbour((*koLibor)[i].date,floater->PayDates,0,floater->PayDates.size()-1,0);
            if(k<0)
                break;
            if(floater->AccrualDates[k] < hitDate){
                (*koLibor)[i].amount = 0.0;                
            }            
        }
    }
    else if (koSettle->getKOStubRule() == "S"){//Not Tested
        k = Neighbour(hitDate, floater->AccrualDates, 0, floater->AccrualDates.size()-1, 1);
        if (k>=0){
            // substract accrued so far.  Assuming first fixed Leg.
            (*koLibor)[0].amount -= floater->getAccrued(hitDate, hitDate, discount)
                                  * discount->pv(hitDate, floater->PayDates[k]);
        }
    }    
    if (isReverse){
        for(i=0;i<(*koLibor).size();i++){
            (*koLibor)[i].amount *= -1.0;
        }
    }
    return koLibor;
}

// return true if there is still cash flow even KO occurs on hitDate.
// return pv of next coupon. (present = hitDate)
// Also, return cash flow, which are not cancelled hitDate KO occurs.  
// If the accrue end (in KOSettle) is past, it pick up as coupon to be paid even "N" case.
bool KOLiborLeg::getCashFlowOnKO(const DateTime hitDate,
                                   const string koStubRule,
                                   double* pvValue,
                                   CashFlowArray*  coupon){
    static const string method("KOLiborLeg::getCashFlowOnKO");
    try {
        //validation
        floater->checkMarket();

        // initializatino
        coupon->clear();

        if (hitDate >= floater->PayDates[floater->PayDates.size()-1])
            return false;

        if (hitDate < koSettle->getFirstAccrueDate() || 
            hitDate > koSettle->getLastAccrueDate())
            return false;
        
        DateTime settleDate = hitDate;
        koSettle->getSettleDate(hitDate, settleDate);

        CashFlowArrayConstSP myLibor = floater->getCashFlowArray();

        // Usually, coupon are considered as "from Include / to Exclude". 
        // i.e. today = Accrue Date means it's now accrue start date and previous coupon are in past.
        // Thus, when KO occurs on Accrue End Date, koStub Rule do nothing,
        // so need to add its payment for any case.
        // "B" case, it pays full cpn of next at later.
        
        int i;
        // accrued coupon.
        CashFlow cfl;            
        if (koStubRule == "S"){
            DateTime acrEnd = koSettle->getAccrueUpToSettle() ? settleDate : hitDate;
            // be careful.  getAccrued return 0 if cfl.date = payment date.
            // So need to add all coupons whose payment date > hitDate and <= settleDate);
            int k = Neighbour(acrEnd, floater->PayDates, 0, floater->getSize()-1, 0);
            for (i=k; i>=0; i--){//find out the first PayDate > hitDte
                if((*myLibor)[i].date<=hitDate){
                    i++;
                    break;
                }
            }
            for (;i<=k && i>=0;i++)
                coupon->push_back((*myLibor)[i]);            
            cfl.amount = floater->getAccrued(valDate, acrEnd, floater->couponCurve.get());              
            cfl.date = settleDate;            
            coupon->push_back(cfl);
        }
        else if (koStubRule == "B"){
            for(i=0;i<floater->AccrualDates.size()-1;i++){
                if (floater->AccrualDates[i]<hitDate && hitDate<=floater->AccrualDates[i+1]){
                    if (hitDate == floater->PayDates[i])
                        cfl.amount = 0.0;      // already paid.
                    else
                        cfl.amount = (*myLibor)[i].amount;
                    cfl.date = settleDate;
                    coupon->push_back(cfl);                
                    break;
                }
            }
        }
        else if (koStubRule =="N"){
            // no coupon 
        }
        else{
            throw ModelException(method, " invalid koStubRule request.");
        }
        
        // calculate PV value.
        *pvValue = 0.0;
        if (coupon->size()>0){
            for (i=0;i<coupon->size();i++)
                *pvValue += (*coupon)[i].amount *  floater->couponCurve->pv(hitDate, (*coupon)[i].date);
            return true;
        }
        else
            return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// with internal value date version
// Need to call getMarket to set valueDate.
double KOLiborLeg::getKOPV(const DateTime&   hitDate, 
                         const string koStubRule) {
    static const string method("KOLiborLeg::getKOPV");
    try {
        floater->checkMarket();
        return getKOPV(valDate, hitDate,  floater->couponCurve.get(), koStubRule);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// get pv of cashflows which will be canceled KO on hitDate.
// here, present = hitDate, not base Date. baseDate is also used in in disocunt.
// If Bond type (B), remove the next coupon fully.
// if Swap type (S), remove accrued so far.  For "S", need to set up KOSettle class.
// otherwise, return total pv of leg value.
double KOLiborLeg::getKOPV(const DateTime&   baseDate, 
                           const DateTime&   hitDate, 
                           const YieldCurve* discount,
                           const string koStubRule) {
    static const string method("KOLiborLeg::getKOPV");
    try {
        
        double pv = 0.0;
        int i;
        // calculate Libor Leg.	(Should be from baseDate, to get Floating Rate)
        // libor contains just cash flow as of baseDate, but not discounted.
        CashFlowArrayConstSP libor = floater->getCashFlowArray(baseDate, discount);        

        // libor include all cash flow include the past. Here, check it this true.
        // if it's not true, memory access error would occurs.
        if (libor->size() != floater->AccrualDates.size()-1)
            throw ModelException(method,"libor size is not same to input array.");

        DateTime settleDate = hitDate;
        if (koSettle->hasSettle()){
            koSettle->getSettleDate(hitDate, settleDate);
        }
        // if the base Date is later than the discount date then the present value is 0
        if (libor->size()==0){
            return 0;
        }
        else {
            for (i = 0; i < libor->size(); i++)
            {// pv the cash flow, include future flows (accrue ended)
                // i.e. including current period coupon
                if (floater->AccrualDates[i+1] > hitDate)
                    pv += (*libor)[i].amount * discount->pv(hitDate, (*libor)[i].date);
            }
            
            if (koStubRule == "B"){
                // substract coupons, which are already accrue started and not paid.
                for(i=0; i<floater->AccrualDates.size()-1; i++){
                    if (floater->AccrualDates[i] < hitDate && hitDate < floater->AccrualDates[i+1]){
                        if (!koSettle->hasSettle()) // no change at settle date.
                            pv -= (*libor)[i].amount * discount->pv(hitDate, (*libor)[i].date);
                        else
                            pv -= (*libor)[i].amount * discount->pv(hitDate, settleDate);
                    }
                }
            }
            else{
                if (koStubRule == "S"){
                    if (koSettle->hasSettle() == false){
                        throw ModelException(method,
                                             "Need to set KOSettle to use koStubRule = 'S'.");
                    }
                    // be careful.  getAccrued return 0 if cfl.date = payment date,
                    // also count only the one period.  If settlement is later than libor accrue end,
                    // it could miss.  So need to sum up.
                    DateTime accrueDate = koSettle->getAccrueUpToSettle() ? settleDate : hitDate;
                    for (i=0; i<floater->AccrualDates.size()-1; i++){
                        if (hitDate < floater->AccrualDates[i+1] && floater->AccrualDates[i+1] <= accrueDate)
                            pv -= (*libor)[i].amount * discount->pv(hitDate, (*libor)[i].date); 
                    }
                    double accrue = floater->getAccrued(valDate, accrueDate, floater->couponCurve.get());              
                    pv -= accrue * discount->pv(hitDate, settleDate);
                }
            }
            return pv;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



// set the settlement info from product class
//void KOLiborLeg::setKOSettle(const InstrumentSettlementSP settle,
//                           const DateTimeArray payDateArray,
//                           const DateTimeArray accrueDates,
//                           const bool isRollingSettle,
//                           const bool isAccrueUpToSettle){
//    koSettle = KOSettleSP(new KOSettle(settle,payDateArray,accrueDates,isRollingSettle,isAccrueUpToSettle));
//}

// tweak spreads.  isTweak = false means roll the spreads back.
void KOLiborLeg::tweakSpreads(bool isTweak)
{
    floater->tweakSpreads(isTweak);
}



////////////////////////
// for StateVariable ///
////////////////////////

// constructor
KOLiborLegSV::KOLiborLegSV(const LiborLeg::LiborLegSVSP flt,
                                        const KOSettleSP koSettle,
                                        const DateTime valDate,
                                        const YieldCurve* discount):
KOLiborLeg(flt->getLibor(),koSettle,valDate,discount), floaterSV(flt) {
}

double KOLiborLegSV::getKOPV(const DateTime&   hitDate){
    static const string method("KOLiborLeg::KOLiborLegSV::getKOPV");
    try {
        double pv = 0.0;
        int i;
        // calculate Libor Leg.	(Should be from baseDate, to get Floating Rate)
        const CashFlowArray* libor = floaterSV->getCashFlowArray(valDate);

        // libor include all cash flow include the past. Here, check it this true.
        // if it's not true, memory access error would occurs.
        if (libor->size() != floater->AccrualDates.size()-1)
            throw ModelException(method,"libor size is not same to input array.");

        DateTime settleDate = hitDate;
        if (koSettle->hasSettle()){
            koSettle->getSettleDate(hitDate, settleDate);
        }
        // if the base Date is later than the discount date then the present value is 0
        if (libor->size()==0){
            return 0;
        }
        else {
            for (i = 0; i < libor->size(); i++)
            {// pv the cash flow, include future flows (accrue ended)
                // i.e. including current period coupon
                if (floater->AccrualDates[i+1] > hitDate)
                    pv += (*libor)[i].amount * floaterSV->payDatePV(i);
            }
            
            if (koSettle->getKOStubRule() == "B"){
                // substract coupons, which are already accrue started and not paid.
                for(i=0; i<floater->AccrualDates.size()-1; i++){
                    if (floater->AccrualDates[i] < hitDate && hitDate < floater->AccrualDates[i+1]){
                        if (!koSettle->hasSettle()) // no change at settle date.
                            pv -= (*libor)[i].amount * floaterSV->payDatePV(i);
                        else
                            pv -= (*libor)[i].amount * floaterSV->payDatePV(i);
                    }
                }
            }
            else{
                if (koSettle->getKOStubRule() == "S"){
                    throw ModelException(method,
                                            "not yet implemented.");
                    if (koSettle->hasSettle() == false){
                        throw ModelException(method,
                                             "Need to set KOSettle to use koStubRule = 'S'.");
                    }
                    // be careful.  getAccrued return 0 if cfl.date = payment date,
                    // also count only the one period.  If settlement is later than libor accrue end,
                    // it could miss.  So need to sum up.
                    DateTime accrueDate = koSettle->getAccrueUpToSettle() ? settleDate : hitDate;
                    for (i=0; i<floater->AccrualDates.size()-1; i++){
                        if (hitDate < floater->AccrualDates[i+1] && floater->AccrualDates[i+1] <= accrueDate)
                            pv -= (*libor)[i].amount * floaterSV->payDatePV(i); 
                    }
                    double accrue = floater->getAccrued(valDate, accrueDate, floater->couponCurve.get());              
//                    pv -= accrue * discount->pv(hitDate, settleDate);
//                  Need to get PV from settleDat!!
                }
            }
            return pv;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
               
DRLIB_END_NAMESPACE
