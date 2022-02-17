//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondCashFlows.cpp
//
//   Description : Bond implementation where you input the cash flows
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : November 16, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BondCashFlows.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RollingSettlement.hpp"

DRLIB_BEGIN_NAMESPACE
/** Bond implementation where you input the cash flows */

BondCashFlows::~BondCashFlows() {
    // empty
}

/** Change the settlement of the bond if applicable */
void BondCashFlows::setSettlement(SettlementSP newSettle) {
    settle = SettlementSP(newSettle.get());
}


BondCashFlows::BondCashFlows(): Bond(TYPE){
    isValid = false;
    errorMessage = "bond not initialized";

    // initialize optional parameters
    dayCountConvString = "30/360";
    redemptionPct = 1.;
    lastCFIncludesRedemption = true;
    exCouponDates = DateTimeArraySP(   );
}

/** Pull out the component assets & correlations from the market data */
void BondCashFlows::getMarket(const IModel* model, const MarketData* market){
    try {
        settle->getMarket(model, market);

    } 
    catch (exception& e){
        throw ModelException(e, "BondCashFlows::getMarket");
    }
}

DateTime BondCashFlows::getMaturityDate() const
{
    throwIfNotValid();
    return cashFlowsProcessed->back().date;
}

DateTime BondCashFlows::getUnadjMaturityDate() const
{
    throwIfNotValid();
    return cashFlowsProcessed->back().date;
}

CashFlowArraySP BondCashFlows::getCashFlows() const
{
    static const string method = "BondCashFlows::getCashFlows";
    try {  
        throwIfNotValid();
        
        return CashFlowArraySP(cashFlowsProcessed.clone());
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

CashFlowArraySP BondCashFlows::getCashFlows(const DateTime &startDate) const
{
    static const string method = "BondCashFlows::getCashFlows_DateTime";

    try { 
        throwIfNotValid();
        CashFlowArraySP myCashFlows(new CashFlowArray());
        int i;
        int numCFs;

        if (startDate <= getMaturityDate()) {

            if (startDate >= exCouponDatesProcessed->back()) {
                // special case. return principal but no coupon
                myCashFlows->resize(1);
                (*myCashFlows)[0].date = getMaturityDate();
                (*myCashFlows)[0].amount = getRedemption();
            } else {
                // count the cashflows after the startDate
                numCFs = 0;
                while(numCFs < cashFlowsProcessed->size() 
                      && (*exCouponDatesProcessed)[cashFlowsProcessed->size() - numCFs - 1] > startDate) {
                    numCFs += 1;
            }
                myCashFlows->resize(numCFs);
                for (i = 0; i<numCFs; i++) {
                    (*myCashFlows)[numCFs-i-1] = (*cashFlowsProcessed)[cashFlowsProcessed->size()-i-1];
                }
            }
        }

        return myCashFlows;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double BondCashFlows::getAccruedAtDate(const DateTime &aiDate) const 
{
    double accrued;
    int idxNextCF;
    
    throwIfNotValid();

    if ( aiDate.getDate() <= datedDate.getDate() || 
         aiDate.getDate() >= cashFlowsProcessed->back().date.getDate() ) {
        accrued = 0.;
    } else {
        // figure out where we fall in the schedule
        idxNextCF = 0;
        while((*cashFlowsProcessed)[idxNextCF].date.getDate() <= aiDate.getDate()) {
            idxNextCF += 1;
        }
        
        if (idxNextCF == 0) {
            if (cashFlowsProcessed->size() == 1) { // special treatment for zeros
                if ((*accruedDCC).days(datedDate, (*cashFlowsProcessed)[idxNextCF].date, false, false)!= 0) {
                    accrued = ((*cashFlowsProcessed)[idxNextCF].amount - faceValue*redemptionPct) * 
                        (*accruedDCC).days(datedDate, aiDate, false, false)/
                        (*accruedDCC).days(datedDate, (*cashFlowsProcessed)[idxNextCF].date, false, false);
                } else {
                    accrued = 0.0;
                }
            
                if (aiDate.getDate() >= (*exCouponDatesProcessed)[idxNextCF].getDate()) {
                    // subtract off the next coupon amount if it's already gone ex
                    accrued -= ((*cashFlowsProcessed)[idxNextCF].amount- faceValue*redemptionPct);
                }
            } else {
                if ((*accruedDCC).days(datedDate, (*cashFlowsProcessed)[idxNextCF].date, false, false) != 0) {
                    accrued = (*cashFlowsProcessed)[0].amount * 
                        (*accruedDCC).days(datedDate, aiDate, false, false)/
                        (*accruedDCC).days(datedDate, (*cashFlowsProcessed)[idxNextCF].date, false, false);
                } else {
                    accrued = 0.0;
                }

                if (aiDate.getDate() >= (*exCouponDatesProcessed)[idxNextCF].getDate()) {
                    // subtract off the next coupon amount if it's already gone ex
                    accrued -= (*cashFlowsProcessed)[idxNextCF].amount;
                }
            }
        } else if (idxNextCF == cashFlowsProcessed->size()-1) {
            if ((*accruedDCC).days((*cashFlowsProcessed)[idxNextCF-1].date, (*cashFlowsProcessed)[idxNextCF].date, false, false) != 0) {
                accrued = ((*cashFlowsProcessed)[idxNextCF].amount - faceValue*redemptionPct) * 
                (*accruedDCC).days((*cashFlowsProcessed)[idxNextCF-1].date, aiDate, false, false)/
                (*accruedDCC).days((*cashFlowsProcessed)[idxNextCF-1].date, (*cashFlowsProcessed)[idxNextCF].date, false, false);
            } else {
                accrued = 0.0;
            }

            if (aiDate.getDate() >= (*exCouponDatesProcessed)[idxNextCF].getDate()) {
                // subtract off the next coupon amount if it's already gone ex
                accrued -= ((*cashFlowsProcessed)[idxNextCF].amount - faceValue*redemptionPct);
            }
        } else {
            if ((*accruedDCC).days((*cashFlowsProcessed)[idxNextCF-1].date, (*cashFlowsProcessed)[idxNextCF].date, false, false)!= 0) {
                accrued = (*cashFlowsProcessed)[idxNextCF].amount * 
                (*accruedDCC).days((*cashFlowsProcessed)[idxNextCF-1].date, aiDate, false, false)/
                (*accruedDCC).days((*cashFlowsProcessed)[idxNextCF-1].date, (*cashFlowsProcessed)[idxNextCF].date, false, false);
            } else {
                accrued = 0.0;
            }

            if (aiDate.getDate() >= (*exCouponDatesProcessed)[idxNextCF].getDate()) {
                // subtract off the next coupon amount if it's already gone ex
                accrued -= (*cashFlowsProcessed)[idxNextCF].amount;
            }
        }
    }
    
    return accrued;
}
    
double BondCashFlows::getMaturityCoupon() const
{
    double matCoupon;

    throwIfNotValid();

    matCoupon = cashFlowsProcessed->back().amount - redemptionPct*faceValue;

    return matCoupon;
}

DateTimeArraySP BondCashFlows::getExCouponDates() const
{
    throwIfNotValid();
 
    return exCouponDatesProcessed;
}

double BondCashFlows::getFaceValue() const {
    return faceValue;
}

double BondCashFlows::getRedemption() const {
    return redemptionPct * faceValue;
}

void BondCashFlows::validatePop2Object()
{
    // try{

    // convert the day count convention string to an object
    accruedDCC = DayCountConventionSP(DayCountConventionFactory::make(dayCountConvString));
      
     // if no settlement (as we only introduced it in May 2004), invent one
    if (!settle.get()) {
        settle = SettlementSP(new RollingSettlement);
    }
 
    validateInputs();
        
    CashFlowArraySP cashFlowsCopy(cashFlows);

    if (lastCFIncludesRedemption == false) {
       cashFlowsCopy->back().amount += redemptionPct*faceValue;
    }

    cashFlowsProcessed = cashFlowsCopy;
        
    if (!exCouponDates || exCouponDates->empty() == true) {
        // exCouponDates were not set
        DateTimeArraySP exCDs(new DateTimeArray);
        int i;
    
        exCDs->resize(cashFlows->size());
        for (i=0; i < cashFlows->size(); i++) {
            (*exCDs)[i] = (*cashFlows)[i].date;
        }
        exCouponDatesProcessed = exCDs;
    } else {
        // they were set
        exCouponDatesProcessed = exCouponDates;
    }
        
    isValid = true;
    // }
    // catch (exception& e) {
    //     errorMessage = ModelException(e).stackTrace();
    //    errorMessage = errorMessage + "\n" + "BondCashFlows::validatePop2Object: Error in object. Continuing but will fail if used";
    //    isValid = false;
    // }
    
    return;
}

void BondCashFlows::throwIfNotValid() const
{
    if (isValid == false) {
        throw ModelException("BondCashFlows::throwIfNotValid", 
                             "attempting to use invalid bond (" + errorMessage + ")");
    }
}

void BondCashFlows::validateInputs() const
{
    static const string method = "BondCashFlows::validateInputs";
    try {
        if (faceValue < 0.) {
            throw ModelException(method, "face value (" + Format::toString(faceValue) + ") cannot be negative");
        } 

        // validate that there are some cashflows
        if (cashFlows->empty()){
            throw ModelException(method, "CashFlow array is empty");
        } 
        // validate that cashFlow dates are increasing 
        int i;
        for (i = 1; i < cashFlows->size(); i++){
            if ((*cashFlows)[i-1].date.isGreater((*cashFlows)[i].date)){
                throw ModelException(method, "CashFlow dates are not increasing");
            }
        }
        // validate that cashFlows are positive
        for (i = 0; i < cashFlows->size(); i++){
            if ((*cashFlows)[i].amount < 0.){
                throw ModelException(method, "CashFlows must be positive");
            }
        }
        if (!(!exCouponDates) && !(exCouponDates->empty())) {
            // validate that number of ex-div dates and cashFlows are the same
            if (exCouponDates->size() != cashFlows->size()) {
                throw ModelException(method, "CashFlow and exCoupon arrays are not the same size");
            }
            // validate that exDivDates are increasing
            for (i = 1; i < exCouponDates->size(); i++){
                if ((*exCouponDates)[i-1].isGreater((*exCouponDates)[i])){
                    throw ModelException(method, "exCoupon dates are not increasing");
                }
            }
            // validate that last exDivDate is on or before last cashFlowDate
             if (exCouponDates->back() > cashFlows->back().date) {
                throw ModelException(method, "last exCouponDate cannot be after the last cashFlow date");
            }
        }
        // validate that datedDate is before first cashFlowDate
        if (datedDate.isGreater((*cashFlows)[0].date)) {
            throw ModelException(method, "datedDate must not be after first cashFlowDate");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


double BondCashFlows::priceFromYield(double yield, bool clean, DateTime yDate) const
{
   static const string method = "CDRBondCashFlows::priceFromYield";
   try {
       double myPrice = 0.;
       if (yDate > getMaturityDate()) {
           throw ModelException("startDate cannot be after maturityDate");
       }
       CashFlowArraySP cfs = getCashFlows(yDate);
       int i;

       // very simple calc. It works for 30/360 where accruedDCC->years will come out to
       // the right fraction, but it will be off a little bit for act/365 and such.
       // It also doesn't do the right thing for long coupons, but you'ld need to take
       // in pseudo coupons for that. As cashflows bonds are rarely used and YTM is not
       // critical, improve only if extremely bored or someone starts complaining.
       // -Tycho
       if (cfs->size() > 1) {

          for (i=cfs->size()-1; i>0; i--) {
                myPrice += (*cfs)[i].amount;
                myPrice /= (1+yield*accruedDCC->years((*cfs)[i-1].date, (*cfs)[i].date)); // should this have eom or bad day adjustments
            }
       }

       if (cashFlowsProcessed->size() > 1) { // we're not a zero
           myPrice += (*cfs)[0].amount;
           DateTime previousDate;
           if (cashFlowsProcessed->size() > cfs->size()) {
               previousDate = (*cashFlowsProcessed)[cashFlowsProcessed->size()-cfs->size()-1].date;
           } else {
               previousDate = datedDate;
           }
           if (cfs->size() > 1) {
               myPrice /= pow((1+yield*accruedDCC->years(previousDate, (*cfs)[0].date)),((*cfs)[0].amount-getAccruedAtDate(yDate))/(*cfs)[0].amount); // not right for long front
           } else {
               myPrice /= pow((1+yield*accruedDCC->years(previousDate, (*cfs)[0].date)),((*cfs)[0].amount-getRedemption()-getAccruedAtDate(yDate))/((*cfs)[0].amount-getRedemption())); // not right for long front
           }
       } else {
           // a zero -- quote semiannually
           myPrice = (*cfs)[0].amount;
           myPrice /= pow((1+yield/2.), (*accruedDCC).years(yDate, (*cfs)[0].date)*2.);
       }

       if (clean == true) {
           myPrice -= getAccruedAtDate(yDate);
       }

       return myPrice;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

Bond* BondCashFlows::getBondToPut(const DateTime &putDate, double putLevel) const
{
   throwIfNotValid();
    static const string method = "BondCashFlows::getBondToPut";
    try {
        BondCashFlowsSP putBond(dynamic_cast<BondCashFlows*>(clone()));
        
        putBond->isValid = false;
        putBond->cashFlowsProcessed = CashFlowArraySP(   );
        putBond->exCouponDatesProcessed = DateTimeArraySP(   );

        if (putDate.getDate() == getMaturityDate().getDate()) {
            putBond->redemptionPct = Maths::max(redemptionPct, putLevel/faceValue);
            (putBond->cashFlows)->back().amount = getMaturityCoupon();
            putBond->lastCFIncludesRedemption = false;
        } else {
            // pop off the end ones while they're after the put date
            int cfIdx;
            for (cfIdx=putBond->cashFlows->size()-1; cfIdx>=0; cfIdx--) {
                if (putDate.getDate() < (*(putBond->cashFlows))[cfIdx].date.getDate()) {
                    putBond->cashFlows->pop_back();
                    if (!(!(putBond->exCouponDates)) && putBond->exCouponDates->empty() == false) { 
                        putBond->exCouponDates->pop_back();
                    }
                } else {
                    break;
                }
            }
        
            if (cfIdx >= 0 && putDate.getDate() == (*(putBond->cashFlows))[cfIdx].date.getDate()) {
                (*(putBond->cashFlows))[cfIdx].amount += putLevel;
            } else {
                putBond->cashFlows->push_back(CashFlow(putDate, putLevel));
                putBond->lastCFIncludesRedemption = true;
                if(!(!(putBond->exCouponDates)) && putBond->exCouponDates->empty() == false) {
                    if ((*exCouponDates)[cfIdx+1].getDate() < putDate.getDate())
                        putBond->exCouponDates->push_back((*exCouponDates)[cfIdx+1]);
                    else {
                        putBond->exCouponDates->push_back(putDate);
                    }
                }
            }
        }

        putBond->validatePop2Object();
        
        return putBond.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DateTime BondCashFlows::getAccrualStartDate() const
{
    return datedDate;
}

// when does bond settle?
DateTime BondCashFlows::settles(const DateTime& tradeDate) const {
    return settle->settles(tradeDate);
}

class BondCashFlowsHelper{
public:
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BondCashFlows, clazz);
        SUPERCLASS(Bond);
        EMPTY_SHELL_METHOD(defaultBondCashFlows);
        FIELD(faceValue, "face value");
        FIELD(cashFlows, "the bond's cash flows");
        FIELD(datedDate, "date interest begins to accrue");
        FIELD(dayCountConvString, "day count convention for accrued interest");
        FIELD_MAKE_OPTIONAL(dayCountConvString);
        FIELD(redemptionPct, "redemption level of the bond as percent of faceValue");
        FIELD_MAKE_OPTIONAL(redemptionPct);
        FIELD(lastCFIncludesRedemption, "Does the last cash flow include the redemption level");
        FIELD_MAKE_OPTIONAL(lastCFIncludesRedemption);
        FIELD(exCouponDates, "the dates the coupons go ex");
        FIELD_MAKE_OPTIONAL(exCouponDates);
        FIELD(settle, "how bond settles");
        FIELD_MAKE_OPTIONAL(settle);
        // private members that need to be copied when bond is copied
        FIELD(isValid, "private member");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(errorMessage, "private member");
        FIELD_MAKE_TRANSIENT(errorMessage);
        FIELD(cashFlowsProcessed, "private member");
        FIELD_MAKE_TRANSIENT(cashFlowsProcessed);
        FIELD(exCouponDatesProcessed, "private member");
        FIELD_MAKE_TRANSIENT(exCouponDatesProcessed);
        FIELD(accruedDCC, "private member");
        FIELD_MAKE_TRANSIENT(accruedDCC);

        // addin constructor
        Addin::registerConstructor("BOND_CASH_FLOWS",
                                   Addin::MARKET,
                                   "Creates a Bond specified via cash flows",
                                   BondCashFlows::TYPE);
    }

    static IObject* defaultBondCashFlows(){
        return new BondCashFlows();
    }
};

CClassConstSP const BondCashFlows::TYPE = CClass::registerClassLoadMethod(
    "BondCashFlows", typeid(BondCashFlows), BondCashFlowsHelper::load);
bool  BondCashFlowsLoad() {
    return (BondCashFlows::TYPE != 0);
   }


DRLIB_END_NAMESPACE

