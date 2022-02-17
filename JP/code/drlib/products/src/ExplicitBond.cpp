//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ExplicitBond.cpp
//
//   Description : Explicit Bond implementation
//
//   Author      : André Segger
//
//   Date        : 10 September 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ExplicitBond.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/ExpiryWindow.hpp"

DRLIB_BEGIN_NAMESPACE
/** Bond implementation where you input the cash flows */

#define SHORT_STUB_THRESHOLD       7

ExplicitBond::ExplicitBond(): Bond(TYPE)
{
}

void ExplicitBond::fieldsUpdated(const CFieldArray& fields) {
    calculateLiborRates();
}


/** Change the settlement of the bond if applicable */
void ExplicitBond::setSettlement(SettlementSP newSettle) {
    settle = SettlementSP(newSettle.get());
}

CashFlowArraySP ExplicitBond::getOutstandingNotionals()
{
    CashFlowArraySP cfs(new CashFlowArray(resetDate->size()));
    for (int i=0; i< cfs->size() ; ++i) {
        (*cfs)[i].date   = (*resetDate)[i];
        (*cfs)[i].amount = (*notional)[i];
    }
    return cfs;
}

void ExplicitBond::getMarket(const IModel* model, const MarketData* market)
{
    try {
        market->GetReferenceDate(valueDate);
        liborReferenceCurve.getData(model, market);
        validatePop2Object();
        settle->getMarket(model, market);

        if (!liborReferenceCurve) {
            throw ModelException("ExplicitBond::getMarket", 
                "The Libor reference curve must not be null for bonds with floating coupons");
        }

        calculateLiborRates(/*false*/);
    } catch (exception& e) {
        throw ModelException(e, "FloatingBond::getMarket");
    }
}

void ExplicitBond::validatePop2Object()
{
    static const string method = "ExplicitBond::validatePop2Object";

    // check that all required arrays exist
    if (!resetDate) {
        throw ModelException(method, "reset date array must not be null");
    }
    if (!initialAccrualDate) {
        throw ModelException(method, "initial accrual date array must not be null");
    }
    if (!finalAccrualDate) {
        throw ModelException(method, "final acrrual date array must not be null");
    }
    if (!exCouponDate) {
        throw ModelException(method, "ex-coupon date array must not be null");
    }
    if (!paymentDate) {
        throw ModelException(method, "payment date array must not be null");
    }
    if (!notional) {
        throw ModelException(method, "notional array must not be null");
    }
    if (!couponRate) {
        throw ModelException(method, "coupon rate array must not be null");
    }
    if (!isFixed) {
        throw ModelException(method, "fixed/floating array must not be null");
    }
    if (!liborSpread) {
        throw ModelException(method, "Libor spread array must not be null");
    }
    if (!redemption) {
        throw ModelException(method, "redemption array must not be null");
    }
    if (!accrualDCC) {
        throw ModelException(method, "acrrual DCC array must not be null");
    }
    if (!paymentDCC) {
        throw ModelException(method, "payment DCC array must not be null");
    }


    // check that all arrays have equal size
    if (resetDate->size() != initialAccrualDate->size() ) {
        throw ModelException(method, "Must have same number of initial accrual dates as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != finalAccrualDate->size() ) {
        throw ModelException(method, "Must have same number of final  accrual dates as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != paymentDate->size() ) {
        throw ModelException(method, "Must have same number of payment dates as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != notional->size() ) {
        throw ModelException(method, "Must have same number of notionals as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != couponRate->size() ) {
        throw ModelException(method, "Must have same number of coupon rates as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != isFixed->size() ) {
        throw ModelException(method, "Must have same number of fixed/floating flags as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != liborSpread->size() ) {
        throw ModelException(method, "Must have same number of libor spreads as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (redemptionDate->size() != redemption->size() ) {
        throw ModelException(method, "Must have same number of redemption amounts as redemption dates (" +
            Format::toString(resetDate->size()) +").");
    }
    if (resetDate->size() != accrualDCC->size() ) {
        throw ModelException(method, "Must have same number of accrual DCCs (" + 
            Format::toString(accrualDCC->size()) + ") as reset dates (" +
            Format::toString(resetDate->size()) +".");
    }
    if (resetDate->size() != paymentDCC->size() ) {
        throw ModelException(method, "Must have same number of payment DCCs (" + 
            Format::toString(paymentDCC->size()) + " as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }
    // Many tests have zero length exCouponDate arrays 
    if (exCouponDate->size()!=0 && resetDate->size() != exCouponDate->size()) {
        throw ModelException(method, "Must have same number of ex coupon dates (" +
            Format::toString(exCouponDate->size()) + " as reset dates (" +
            Format::toString(resetDate->size()) +").");
    }

    if (!!couponCap) {
        if (resetDate->size() != couponCap->size() ) {
            throw ModelException(method, "Must have same number of coupon caps as reset dates (" +
                Format::toString(resetDate->size()) +".");
        }
    }

    if (!!couponFloor) {
        if (resetDate->size() != couponFloor->size() ) {
            throw ModelException(method, "Must have same number of coupon floors as reset dates (" +
                Format::toString(resetDate->size()) +".");
        }
    }

    if (resetDate->size() == 0 && redemptionDate->size() == 0) {
        throw ModelException(method, "Must have at least one coupon or redemption payment");
    }


    // create transient fields
    if (!coupons) {
        coupons = DoubleArraySP(new DoubleArray(resetDate->size()));
    }
    if (!rateToUse) {
        rateToUse = DoubleArraySP(new DoubleArray(resetDate->size()));
    }
    if (!forwards) {
        forwards = DoubleArraySP(new DoubleArray(resetDate->size()));
    }
    if (!dayCountFracs) {
        dayCountFracs = DoubleArraySP(new DoubleArray(resetDate->size()));
    }

    int  i;
    bool hasFloatingCoupons = false;
    for (i=1; i<resetDate->size() ; ++i) {
        // check that reset dates are in ascending order
        if ((*resetDate)[i] < (*resetDate)[i-1]) {
            throw ModelException(method, "Reset dates are not in strictly ascending order.");
        }
        // check that initial accrual dates are in ascending order
        if ((*initialAccrualDate)[i] < (*initialAccrualDate)[i-1]) {
            throw ModelException(method, "Initial accrual dates are not in strictly ascending order.");
        }
        // check that final accrual dates are in ascending order
        if ((*finalAccrualDate)[i] < (*finalAccrualDate)[i-1]) {
            throw ModelException(method, "Final accrual dates are not in strictly ascending order.");
        }
    }

    for (i=0; i<resetDate->size() ; ++i) {
        // check that initial accrual date is on or after reset date
        if ((*resetDate)[i] > (*initialAccrualDate)[i]) {
            throw ModelException(method, "Initial accrual dates must be after reset dates.");
        }
        // check that final accrual date after initial accrual date
        if ((*finalAccrualDate)[i] <= (*initialAccrualDate)[i]) {
            throw ModelException(method, "Final accrual dates must be after initial accrual dates .");
        }
        // check that required data has been set for floating coupons+
        if (!(*isFixed)[i]) {
            hasFloatingCoupons = true;
        }
    }

    if (hasFloatingCoupons) {
        if (!liborPeriod) {
            throw ModelException(method, "The Libor period must not be null for bonds with floating coupons");
        }
    }

    // if no settlement (as we only introduced it in May 2004), invent one
    if (!settle.get()) {
        settle = SettlementSP(new RollingSettlement);
    }

    // convert the day count convention string to an object
    liborDCC  = DayCountConventionSP(DayCountConventionFactory::make(liborDCCString));

    bondAccrualDCC = DayCountConventionArraySP(new DayCountConventionArray(resetDate->size()));
    bondPaymentDCC = DayCountConventionArraySP(new DayCountConventionArray(resetDate->size()));
    for (i=0 ; i<resetDate->size() ; ++i) {
        (*bondAccrualDCC)[i] = DayCountConventionSP(DayCountConventionFactory::make((*accrualDCC)[i]));
        (*bondPaymentDCC)[i] = DayCountConventionSP(DayCountConventionFactory::make((*paymentDCC)[i]));
    }

    if ( initialAccrualDate->size() > 0 ) {
        datedDate = (*initialAccrualDate)[0];
    } else {
        // zero coupon bond
        datedDate = valueDate;
    }
}


void ExplicitBond::calculateLiborRates(/*bool isTheta*/)
{
    // set the floating rates
    DateTime liborEndDate;
    double   libor;
    if ( !valueDate.empty()) {
        for (int idx=0 ; idx<resetDate->size() ; ++idx) {
            if (((*resetDate)[idx] > valueDate) && !(*isFixed)[idx]) {
                if ( idx == 0 || idx == resetDate->size()-1 ) {
                    // check if there's short stub
                    DateTime provisionalEndDate = liborPeriod->toDate((*initialAccrualDate)[idx]);
                    long     dateDiff           = abs(provisionalEndDate.getDate() -
                                                      (*finalAccrualDate)[idx].getDate());
                    
                    if ( dateDiff > SHORT_STUB_THRESHOLD ) {
                        // the libor period does not match the accrual period - calculate the forward 
                        // rate between the accrual dates - this will need to be adjusted for the spot
                        // offset when the spot offset is in the yield curve
                        libor        = liborReferenceCurve->fwd((*initialAccrualDate)[idx],
                                                                (*finalAccrualDate)[idx],
                                                                liborDCC.get(),
                                                                CompoundBasis::SIMPLE);
                    } else {
                        liborEndDate = liborPeriod->toDate((*resetDate)[idx]);
                        libor        = liborReferenceCurve->fwd((*resetDate)[idx],
                                                                liborEndDate,
                                                                liborDCC.get(),
                                                                CompoundBasis::SIMPLE);
                    }
                } else {
                    liborEndDate =  liborPeriod->toDate((*resetDate)[idx]);
                    libor        = liborReferenceCurve->fwd((*resetDate)[idx],
                                                            liborEndDate,
                                                            liborDCC.get(),
                                                            CompoundBasis::SIMPLE);
                    
                }
            } else {
                if ( (*isFixed)[idx] ) {
                    libor        = (*couponRate)[idx];
                } else {
                    libor        = (*couponRate)[idx];
                }
            }

            // only apply libor spread to floating coupons
            if ( !(*isFixed)[idx] ) {
                libor += (*liborSpread)[idx];
            }

            // negative coupons will not be paid by the bondholder - floor coupon at 0.0
            libor = Maths::max(libor, 0.0);

            // apply caps and floors to floating coupons
            if ( !(*isFixed)[idx] ) {
                if ( !!couponCap ) {
                    libor = Maths::min(libor, (*couponCap)[idx]);
                }

                if ( !!couponFloor ) {
                    libor = Maths::max(libor, (*couponFloor)[idx]);
                }
            }


            (*forwards)[idx]      = libor;
            (*dayCountFracs)[idx] = (*bondPaymentDCC)[idx]->accruedFactor((*finalAccrualDate)[idx],
                                                                          (*initialAccrualDate)[idx],
                                                                          (*finalAccrualDate)[idx],
                                                                          false,
                                                                          false);

            if ( (*isFixed)[idx] ) {
                (*coupons)[idx] = libor * (*notional)[idx] * (*dayCountFracs)[idx];
            } else {
                (*coupons)[idx] = libor * (*notional)[idx] * (*dayCountFracs)[idx];
            }

            // (*coupons)[idx] += (*redemption)[idx];
        }
    }

    // calculate cash flows
    DateTimeArray paymentDates(paymentDate->size());
    DateTimeArray redemptionDates(redemptionDate->size());
    int i;
    for (i=0 ; i<paymentDate->size() ; ++i) {
        paymentDates[i] = (*paymentDate)[i];
    }
    for (i=0 ; i<redemptionDate->size() ; ++i) {
        redemptionDates[i] = (*redemptionDate)[i];
    }

    // create an array 
    DateTimeArray mergedDates = DateTime::merge(paymentDates,redemptionDates);
    cashFlows = CashFlowArraySP(new CashFlowArray(mergedDates.size()));
    int couponIdx     = 0;
    int redemptionIdx = 0;
    for(i=0 ; i<mergedDates.size() ; ++i) {
        (*cashFlows)[i].date    = mergedDates[i];
        (*cashFlows)[i].amount = 0.0;
        if (couponIdx < paymentDate->size() && (*paymentDate)[couponIdx] == mergedDates[i]) {
            (*cashFlows)[i].amount += (*coupons)[couponIdx];
            ++couponIdx;
        }
        if (redemptionIdx < redemptionDate->size() && (*redemptionDate)[redemptionIdx] == mergedDates[i]) {
            (*cashFlows)[i].amount += (*redemption)[redemptionIdx];
            ++redemptionIdx;
        }
    }
}

CashFlowArraySP ExplicitBond::getCouponRates()
{
    // create an array 
    DateTimeArray mergedDates = DateTime::merge(*(paymentDate.get()),*(redemptionDate.get()));
    CashFlowArraySP couponRates(new CashFlowArray(mergedDates.size()));
    int couponIdx     = 0;
    for(int i=0 ; i<mergedDates.size() ; ++i) {
        (*cashFlows)[i].date     = mergedDates[i];
        (*couponRates)[i].amount = 0.0;
        if (couponIdx < paymentDate->size() && (*paymentDate)[couponIdx] == mergedDates[i]) {
            (*couponRates)[i].amount = (*forwards)[couponIdx];
            ++couponIdx;
        }
    }
    return couponRates;
}

CashFlowArraySP ExplicitBond::getDayCountFractions()
{
    // create an array 
    DateTimeArray mergedDates = DateTime::merge(*(paymentDate.get()),*(redemptionDate.get()));
    CashFlowArraySP dayCounts(new CashFlowArray(mergedDates.size()));
    int couponIdx     = 0;
    for(int i=0 ; i<mergedDates.size() ; ++i) {
        (*cashFlows)[i].date     = mergedDates[i];
        (*dayCounts)[i].amount = 0.0;
        if (couponIdx < paymentDate->size() && (*paymentDate)[couponIdx] == mergedDates[i]) {
            (*dayCounts)[i].amount = (*dayCountFracs)[couponIdx];
            ++couponIdx;
        }
    }
    return dayCounts;
}


DateTime ExplicitBond::getMaturityDate() const
{
    // the greater of the last coupon and last redemption date
    DateTime matDate;
    if ( paymentDate->size() > 0 && redemptionDate->size() > 0 ) {
        if ( paymentDate->back() > redemptionDate->back()) {
            matDate = paymentDate->back();
        } else {
            matDate = redemptionDate->back();
        }
    } else if ( paymentDate->size() > 0 ) {
        matDate = paymentDate->back();
    } else if ( redemptionDate->size() > 0 ) {
        matDate = redemptionDate->back();
    }
    return matDate;
}

DateTime ExplicitBond::getUnadjMaturityDate() const
{
    return getMaturityDate();
}

// returns all the cashflows
CashFlowArraySP ExplicitBond::getCashFlows() const
{
    static const string method = "ExplicitBond::getCashFlows";
    try {  
        CashFlowArraySP cfs = CashFlowArraySP(cashFlows.clone());
        return cfs;
    } catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

/** returns only the cashFlows that occur after the start date */
CashFlowArraySP ExplicitBond::getCashFlows(const DateTime &startDate) const
{
    static const string method = "ExplicitBond::getCashFlows";

    try { 
        CashFlowArraySP myCashFlows(new CashFlowArray());
        int i;
        int numCFs;

        if (startDate <= getMaturityDate()) {

            // count the cashflows after the startDate
            numCFs = 0;
            while(numCFs < cashFlows->size() 
                  && (*cashFlows)[cashFlows->size() - numCFs - 1].date > startDate) {
                ++numCFs;
            }
            myCashFlows->resize(numCFs);
            for (i = 0; i<numCFs; i++) {
                (*myCashFlows)[numCFs-i-1] = (*cashFlows)[cashFlows->size()-i-1];
            }
        }

        return myCashFlows;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// returns all the coupons
CashFlowArraySP ExplicitBond::getCoupons() const
{
    static const string method = "ExplicitBond::getCoupons";
    try {
        CashFlowArraySP couponCF(new CashFlowArray(coupons->size()));
        for (int i=0 ; i<couponCF->size() ; ++i) {
            (*couponCF)[i].date   = (*paymentDate)[i];
            (*couponCF)[i].amount = (*coupons)[i];
        }
        return couponCF;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

}

/** returns only the coupons that occur after the start date */
CashFlowArraySP ExplicitBond::getCoupons(const DateTime &startDate) const
{
    static const string method = "ExplicitBond::getCoupons";
    try {
        CashFlowArraySP couponCF(new CashFlowArray(0));
        int i, numCoupons = 0;
        
        // count the coupons after the startDate
        while (numCoupons < coupons->size()
            && (*paymentDate)[coupons->size() - numCoupons - 1] > startDate) {
             ++numCoupons;
        }

        couponCF->resize(numCoupons);
        for (i = 0; i<numCoupons; i++) {
            (*couponCF)[numCoupons-i-1].date   = (*paymentDate)[coupons->size()-i-1];
            (*couponCF)[numCoupons-i-1].amount = (*coupons)[coupons->size()-i-1];
        }
        return couponCF;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** get the accrued interest at a date */
double ExplicitBond::getAccruedAtDate(const DateTime &aiDate) const
{
    double accrued   =  0.0;
    int    idx;
    
    if ( resetDate->size() == 0 ) {
        accrued = 0.;
    } else if ( aiDate.getDate() <= datedDate.getDate() || 
         aiDate.getDate() >= finalAccrualDate->back().getDate() ) {
        accrued = 0.;
    } else {
        // figure out where we fall in the schedule
        idx = 0;
        while( idx < finalAccrualDate->size() && 
               (*finalAccrualDate)[idx].getDate() <= aiDate.getDate()) {
            ++idx;
        }

        // Has coupon gone ex? -> negative accrued interest
        bool goneEx = false;
        if (exCouponDate->size() > 0) {
            goneEx = ((aiDate >= (*exCouponDate)[idx])? true: false);
        }

        // Get the coupon rate for the current coupon period
        // If bond is floating, the coupon rate will be populated since aiDate is past the reset date
        // If bond is floating need to add the libor spread to rate fixing - to agree with
        // coupon amount calculation in calculateLiborRates()
        double currentCouponRate = (*couponRate)[idx];
        if ( !(*isFixed)[idx] ) {
            currentCouponRate += (*liborSpread)[idx];
        }

        // Assume all coupon periods are regular. Since we have no coupon frequency there is no
        // way of knowing whether a coupon period is irregular
        // Also assume no end of month rule since not supplied
        accrued = currentCouponRate * (*notional)[idx] * 
                  (*bondAccrualDCC)[idx]->accruedFactor(aiDate, 
                                                 (*initialAccrualDate)[idx],
                                                 (*finalAccrualDate)[idx],
                                                 false,
                                                 goneEx);
    }
    
    return accrued;
}

DateTimeArraySP ExplicitBond::getExCouponDates() const
{
    return exCouponDate;
}

double ExplicitBond::getFaceValue() const
{
    return (*notional)[0];
}

CashFlowArraySP ExplicitBond::getRedemptionPayments() const
{
    CashFlowArraySP redemptionPayments(new CashFlowArray(redemptionDate->size()));
    for(int i=0; i<redemptionPayments->size() ; ++i) {
      (*redemptionPayments)[i].date   = (*redemptionDate)[i];
      (*redemptionPayments)[i].amount = (*redemption)[i];
    }

    return redemptionPayments;
}

double ExplicitBond::getRedemption() const
{
    double redemptionAmount = 0.0;
    if ( redemption->size() > 0) {
        redemptionAmount = redemption->back();
    }
    return redemptionAmount;

}

double ExplicitBond::getMaturityCoupon() const
{
    return (coupons->size()==0)?0.0:coupons->back();
}

DateTime ExplicitBond::getAccrualStartDate() const
{
    return datedDate;
}

double ExplicitBond::getNotional(const DateTime &date) const
{
    double currentNotional;
    if ( resetDate->size() > 0 ) {
        if (date < (*initialAccrualDate)[0] || date > getMaturityDate() ) {
            currentNotional = 0.0;
        } else {
            int idx = 1;
            while (idx < initialAccrualDate->size() && (*initialAccrualDate)[idx] < date ) ++idx;
            if ( idx < initialAccrualDate->size() ) {
                currentNotional = (*notional)[idx-1];
            } else {
                currentNotional = notional->back();
            }
        }
    } else {
        // for zero coupon bonds the notional is the sum of all outstanding
        // redemption payments
        currentNotional = 0.0;
        int idx = redemption->size()-1;
        while ( idx >= 0 && (*redemptionDate)[idx] > date ) {
            currentNotional += (*redemption)[idx];
            --idx;
        }
    }

    return currentNotional;
}

bool ExplicitBond::isZeroBond() const
{
    return resetDate->size() == 0;
}


/** required for Yield to Maturity calculation */
double ExplicitBond::priceFromYield(double yield, bool clean, DateTime yDate) const
{
   static const string method = "ExplicitBond::priceFromYield";
   try {
       double myPrice = 0.;

       int i;
       // very simple calc. It works for 30/360 where accruedDCC->years will come out to
       // the right fraction, but it will be off a little bit for act/365 and such.
       // It also doesn't do the right thing for long coupons, but you'd need to take
       // in pseudo coupons for that. As cashflows bonds are rarely used and YTM is not
       // critical, improve only if extremely bored or someone starts complaining. -André

       if (!isZeroBond() || redemptionDate->size() != 1) {
           if (coupons->size() > 1) {
              for (i=coupons->size()-1; i>=0; i--) {
                  if ( (*initialAccrualDate)[i] > yDate || 
                       ((*initialAccrualDate)[i] <= yDate  && (*finalAccrualDate)[i] <= yDate ) ) {
                      // include full coupon
                      myPrice += (*cashFlows)[i].amount;
                      // AM: Not at all nice - this function needs revisiting 
                      myPrice /= (1+yield*(*bondAccrualDCC)[i]->accruedFactor((*paymentDate)[i],
                                                                              yDate,
                                                                              (*paymentDate)[i],
                                                                              false,
                                                                              false));

                  }
              }
           }
       } else {
           // a single redemption zero bond -- quote semiannually
           if ( yDate <= (*redemptionDate)[0] ) {
               myPrice = (*redemption)[0];
               myPrice /= pow((1+yield/2.), (*bondAccrualDCC)[0]->accruedFactor((*redemptionDate)[0],
                                                                                yDate, 
                                                                                (*redemptionDate)[0],
                                                                                false,
                                                                                false) * 2.0);
           }
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

/** get the equivalent bond which matures on the put date. Used for yield to put. */
Bond* ExplicitBond::getBondToPut(const DateTime &putDate, double putLevel) const
{
    static const string method = "ExplicitBond::getBondToPut";
    try {
       ExplicitBondSP putBond(dynamic_cast<ExplicitBond*>(clone()));

       // pop off the coupons which are paid after the put date
       int cfIdx;
       for (cfIdx=putBond->paymentDate->size()-1; cfIdx>=0; cfIdx--) {
           if (putDate.getDate() < (*(putBond->paymentDate))[cfIdx].getDate()) {
               putBond->resetDate->pop_back();
               putBond->initialAccrualDate->pop_back();
               putBond->finalAccrualDate->pop_back();
               putBond->exCouponDate->pop_back();
               putBond->paymentDate->pop_back();
               putBond->notional->pop_back();
               putBond->couponRate->pop_back();
               putBond->isFixed->pop_back();
               putBond->liborSpread->pop_back();
               putBond->accrualDCC->pop_back();
               putBond->paymentDCC->pop_back();
               if (!!putBond->couponCap && !putBond->couponCap->empty()) { 
                   putBond->couponCap->pop_back();
               }
               if (!!putBond->couponFloor && !putBond->couponFloor->empty()) { 
                   putBond->couponFloor->pop_back();
               }
           } else {
               break;
           }
       }

       // remove all redemption payments after the put date
       for (cfIdx=putBond->redemptionDate->size()-1; cfIdx>=0; cfIdx--) {
           if (putDate.getDate() < (*(putBond->redemptionDate))[cfIdx].getDate()) {
               putBond->redemptionDate->pop_back();
               putBond->redemption->pop_back();
           }
       }

       // add the put payment on the put date
        if (putBond->redemptionDate->size() > 0 &&
            putDate.getDate() == (putBond->redemptionDate)->back().getDate()) {
            (putBond->redemption)->back() = Maths::max((putBond->redemption)->back(), putLevel);
        } else {
            putBond->redemptionDate->push_back(putDate);
            putBond->redemption->push_back(putLevel);
        }

        putBond->validatePop2Object();
        
        return putBond.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

bool ExplicitBond::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);

        if (BootstrappedYieldCurve::TYPE->isInstance(liborReferenceCurve.get())) {
            IObject* obj = dynamic_cast<IObject*>(liborReferenceCurve.get());
            BootstrappedYieldCurve* yieldCurve = dynamic_cast<BootstrappedYieldCurve*>(obj);
            yieldCurve->sensShift(shift);

            // check whether there is a fixing today
            for (int i=0; i<resetDate->size() ;++i) {
                if ( !(*isFixed)[i] && (*resetDate)[i] == valueDate && 
                     Maths::isZero((*couponRate)[i])) {

                    if ( i == 0 || i == resetDate->size()-1 ) {
                        DateTime provisionalEndDate = liborPeriod->toDate((*initialAccrualDate)[i]);
                        long     dateDiff           = abs(provisionalEndDate.getDate() -
                                                          (*finalAccrualDate)[i].getDate());

                        if ( dateDiff > SHORT_STUB_THRESHOLD ) {
                            // the libor period does not match the accrual period - calculate the forward 
                            // rate between the accrual dates - this will need to be adjusted for the spot
                            // offset when the spot offset is in the yield curve
                            (*couponRate)[i] = liborReferenceCurve->fwd((*initialAccrualDate)[i],
                                                                    (*finalAccrualDate)[i],
                                                                    liborDCC.get(),
                                                                    CompoundBasis::SIMPLE);
                        } else {
                            DateTime liborEndDate = liborPeriod->toDate((*resetDate)[i]);
                            (*couponRate)[i]      = liborReferenceCurve->fwd((*resetDate)[i],
                                                                    liborEndDate,
                                                                    liborDCC.get(),
                                                                    CompoundBasis::SIMPLE);
                        }
                    } else {
                        // set the libor fixing
                        DateTime liborEndDate = liborPeriod->toDate((*resetDate)[i]);
                        (*couponRate)[i]      = liborReferenceCurve->fwd((*resetDate)[i],
                                                                    liborEndDate,
                                                                    liborDCC.get(),
                                                                    CompoundBasis::SIMPLE);
                    }
                }
            }
            calculateLiborRates(/*true*/);
        } else {
            throw ModelException("ExplicitBond::sensShift", "Reference Libor Curve only supports BootstrappedYieldCurve at the moment.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "ExplicitBond::sensShift(theta)");
    }    
    return false; // our components don't have theta type sensitivity
}

// when does bond settle?
DateTime ExplicitBond::settles(const DateTime& tradeDate) const {
    return settle->settles(tradeDate);
}


class ExplicitBondHelper{
public:
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ExplicitBond, clazz);
        SUPERCLASS(Bond);
        EMPTY_SHELL_METHOD(defaultExplicitBond);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate,             "value date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(resetDate,                    "reset dates");
        FIELD(initialAccrualDate,           "initial accrual dates");
        FIELD(finalAccrualDate,             "final accrual dates");
        FIELD(exCouponDate,                 "ex-coupon dates");
        FIELD(paymentDate,                  "payment dates");
        FIELD(notional,                     "notionals");
        FIELD(redemptionDate,               "redemption dates");
        FIELD(redemption,                   "redemption amounts");
        FIELD(couponRate,                   "coupon rates");
        FIELD(isFixed,                      "fixed = true");
        FIELD(liborSpread,                  "Libor spreads");
        FIELD(couponCap,                    "coupon caps");
        FIELD_MAKE_OPTIONAL(couponCap);
        FIELD(couponFloor,                  "coupon floors");
        FIELD_MAKE_OPTIONAL(couponFloor);
        FIELD(accrualDCC,                   "accrual day count conventions");
        FIELD(paymentDCC,                   "payment day count conventions");
        FIELD(liborPeriod,                  "Libor period");
        FIELD_MAKE_OPTIONAL(liborPeriod);
        FIELD(liborDCCString,        "Libor day count convention");
        FIELD_MAKE_OPTIONAL(liborDCCString);
        FIELD(liborReferenceCurve,   "Libor reference curve");
        FIELD_MAKE_OPTIONAL(liborReferenceCurve);
        FIELD(bondFrequency,         "bond frequency - required for Act/Act ISMA bonds");
        FIELD_MAKE_OPTIONAL(bondFrequency);
        FIELD(settle, "how bond settles");
        FIELD_MAKE_OPTIONAL(settle);

        // transient fields
        FIELD(cashFlows,                            "internal field");
        FIELD_MAKE_TRANSIENT(cashFlows);
        FIELD(coupons,                              "internal field");
        FIELD_MAKE_TRANSIENT(coupons);
        FIELD(forwards,                             "internal field");
        FIELD_MAKE_TRANSIENT(forwards);
        FIELD(dayCountFracs,                        "internal field");
        FIELD_MAKE_TRANSIENT(dayCountFracs);
        FIELD(rateToUse,                            "internal field");
        FIELD_MAKE_TRANSIENT(rateToUse);
        FIELD(liborDCC,                             "internal field");
        FIELD_MAKE_TRANSIENT(liborDCC);
        FIELD(datedDate,                     "internal field");
        FIELD_MAKE_TRANSIENT(datedDate);
        FIELD(bondPaymentDCC,                       "internal field");
        FIELD_MAKE_TRANSIENT(bondPaymentDCC);
        FIELD(bondAccrualDCC,                       "internal field");
        FIELD_MAKE_TRANSIENT(bondAccrualDCC);

        // addin constructor
        Addin::registerConstructor("EXPLICIT_BOND_STREAM",
                                   Addin::MARKET,
                                   "Creates an explicit bond stream",
                                   ExplicitBond::TYPE);
        }

        static IObject* defaultExplicitBond(){
            return new ExplicitBond();
        }
};

CClassConstSP const ExplicitBond::TYPE = CClass::registerClassLoadMethod(
    "ExplicitBond", typeid(ExplicitBond), ExplicitBondHelper::load);
bool   ExplicitBondLoad() {
    return (ExplicitBond::TYPE != 0);
}


DRLIB_END_NAMESPACE

