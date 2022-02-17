//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FloatingBond.cpp
//
//   Description : Floating Bond implementation
//
//   Author      : André Segger
//
//   Date        : 30 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FloatingBond.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/BondCashFlows.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE
/** Bond implementation where you input the cash flows */

FloatingBond::FloatingBondPayment::FloatingBondPayment(
                                    const DateTime& refixDate,
                                    const DateTime& paymentDate,
                                    const double&   fixing,
                                    const bool      isFixed) : CObject(TYPE)
{
    static const string method = "FloatingBond::FloatingBondPayment::FloatingBondPayment";
    
    try {
        this->refixDate    = refixDate;
        this->paymentDate  = paymentDate;
        this->fixing       = fixing;
        this->isFixed      = isFixed;
    }
    catch (exception& e) 
    {
        throw ModelException(&e, method);
    }
}

void FloatingBond::fieldsUpdated(const CFieldArray& fields){
    calculateLiborRates(false);
}


/** Change the settlement of the bond if applicable */
void FloatingBond::setSettlement(SettlementSP newSettle) {
    settle = SettlementSP(newSettle.get());
}

void FloatingBond::FloatingBondPayment::validatePop2Object()
{
    if ( refixDate > paymentDate ) {
        throw ModelException("FloatingBond::FloatingBondPayment::validatePop2Object", 
            "Refix date (" + refixDate.toString() + ") is after payment date (" +
            paymentDate.toString() + ")");
    }
}

/* for reflection */
FloatingBond::FloatingBondPayment::FloatingBondPayment(): CObject(TYPE)
{
    fixing                  = 0.0;
    isFixed                 = false;
}

FloatingBond::~FloatingBond() {
    // empty
}

FloatingBond::FloatingBond(const double                      faceValue,
                           const DateTime                    datedDate,
                           const double                      redemptionPct,
                           const MaturityPeriodSP            liborPeriod,
                           const string                      liborDCCString,
                           const double                      liborSpread,
                           const YieldCurveWrapper&          liborReferenceCurve,
                           const DateTime                    valueDate,
                           const FloatingBondPaymentArray&   fixings):
            Bond(TYPE), faceValue(faceValue), datedDate(datedDate), redemptionPct(redemptionPct),
            liborPeriod(liborPeriod), liborDCCString(liborDCCString), liborSpread(liborSpread),
            liborReferenceCurve(liborReferenceCurve), valueDate(valueDate), deferredCoupon(false),
            couponDeferralThreshold(0.0), capFloatingCoupons(false), floatingCouponCap(0.0),
            floorFloatingCoupons(false), floatingCouponFloor(0.0)
{           
    this->fixings = fixings;
}

FloatingBond::FloatingBond(const double                      faceValue,
                           const DateTime                    datedDate,
                           const double                      redemptionPct,
                           const MaturityPeriodSP            liborPeriod,
                           const string                      liborDCCString,
                           const double                      liborSpread,
                           const YieldCurveWrapper&          liborReferenceCurve,
                           const DateTime                    valueDate,
                           const FloatingBondPaymentArray&   fixings,
                           const bool                        deferredCoupon,
                           const CAssetWrapper               asset,
                           const DateTime                    couponDeferralEndDate,
                           const double                      couponDeferralThreshold):
            Bond(TYPE), faceValue(faceValue), datedDate(datedDate), redemptionPct(redemptionPct),
            liborPeriod(liborPeriod), liborDCCString(liborDCCString), liborSpread(liborSpread),
            liborReferenceCurve(liborReferenceCurve), valueDate(valueDate), deferredCoupon(deferredCoupon),
            asset(asset), couponDeferralEndDate(couponDeferralEndDate),
            couponDeferralThreshold(couponDeferralThreshold), capFloatingCoupons(false), floatingCouponCap(0.0),
            floorFloatingCoupons(false), floatingCouponFloor(0.0)
{           
    this->fixings = fixings;
}

FloatingBond::FloatingBond(): Bond(TYPE), deferredCoupon(false), couponDeferralThreshold(0.0),
                              capFloatingCoupons(false), floatingCouponCap(0.0), 
                              floorFloatingCoupons(false), floatingCouponFloor(0.0)
{
    isValid = false;

    // initialize optional parameters
    redemptionPct = 1.;
    liborDCCString = "ACT/360";
}

DateTime FloatingBond::getMaturityDate() const
{
    throwIfNotValid();
    return cashFlowsProcessed->back().date;
}

DateTime FloatingBond::getUnadjMaturityDate() const
{
    throwIfNotValid();
    return cashFlowsProcessed->back().date;
}

CashFlowArraySP FloatingBond::getCashFlows() const
{
    static const string method = "FloatingBond::getCashFlows";
    try {  
        throwIfNotValid();
        
        return CashFlowArraySP(cashFlowsProcessed.clone());
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

CashFlowArraySP FloatingBond::getCashFlows(const DateTime &startDate) const
{
    static const string method = "FloatingBond::getCashFlows";

    try { 
        throwIfNotValid();
        CashFlowArraySP myCashFlows(new CashFlowArray(0));
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

double FloatingBond::getAccruedAtDate(const DateTime &aiDate) const 
{
    static const string method = "FloatingBond::getAccruedAtDate";
    try {

        double accrued=0.0;
        int idxNextCF;
    
        throwIfNotValid();

        if ( aiDate.getDate() <= datedDate.getDate() || 
             aiDate.getDate() >= cashFlowsProcessed->back().date.getDate() ) {
            accrued = 0.0;
        } else {
            // figure out where we fall in the schedule
            idxNextCF = 0;
            while((*cashFlowsProcessed)[idxNextCF].date.getDate() <= aiDate.getDate()) {
                idxNextCF += 1;
            }
        
            // Get the coupon rate for the current coupon period. Normally the fixing for the current coupon period
            // is known because valueDate is generally past the reset date. However the CVB calls getAccruedAtDate()
            // for all future time points so need to estimate future coupon rates
            double currentCouponRate = 0.0;
            if (!!(fixings[idxNextCF]->isFixed) || 
                 fixings[idxNextCF]->refixDate < valueDate) {
                // Fixing is known
                currentCouponRate = fixings[idxNextCF]->fixing;
            } else {
                // Fixing not yet known
                currentCouponRate = liborReferenceCurve->fwd(fixings[idxNextCF]->refixDate,
                                                             liborPeriod->toDate(fixings[idxNextCF]->refixDate),
                                                             liborDCC.get(),
                                                             CompoundBasis::SIMPLE);
            }
            
            // Add the libor spread if floating to agree with coupon amount calculation in calculateLiborRates()
            if ( !fixings[idxNextCF]->isFixed ) {
                currentCouponRate += liborSpread;
            }

            // Determine parameters needed by accruedFactor()
            DateTime prevDate, nextDate;
            bool eomAdjSec = false;         // Default to false since not supplied
            bool goneEx = false;

            if (idxNextCF == 0) {
                prevDate = datedDate;
            }
            else { 
                prevDate = (*cashFlowsProcessed)[idxNextCF-1].date;
            }    
            nextDate = (*cashFlowsProcessed)[idxNextCF].date;

            if (aiDate >= (*exCouponDatesProcessed)[idxNextCF]) {
                // coupon has gone ex
                goneEx = true;
            }
    
            if (prevDate == nextDate) {
                accrued = 0.0;
            }
            else {
                accrued = currentCouponRate * faceValue * liborDCC->accruedFactor(aiDate,
                                                                                  prevDate,
                                                                                  nextDate,
                                                                                  eomAdjSec,
                                                                                  goneEx);
            }
        }
    
        return accrued;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
    
double FloatingBond::getMaturityCoupon() const
{
    double matCoupon;

    throwIfNotValid();

    matCoupon = cashFlowsProcessed->back().amount - redemptionPct*faceValue;

    return matCoupon;
}

DateTimeArraySP FloatingBond::getExCouponDates() const
{
    throwIfNotValid();
 
    return exCouponDatesProcessed;
}

double FloatingBond::getFaceValue() const {
    return faceValue;
}

double FloatingBond::getRedemption() const {
    return redemptionPct * faceValue;
}


void FloatingBond::validatePop2Object()
{
    // convert the day count convention string to an object
    liborDCC = DayCountConventionSP(DayCountConventionFactory::make(liborDCCString));
        
    validateInputs();

    cashFlowsProcessed = CashFlowArraySP( new CashFlowArray(fixings.size()));
    DateTimeArraySP exCDs(new DateTimeArray(fixings.size()));

    for (int i=0; i < fixings.size(); i++) {
        (*cashFlowsProcessed)[i].date   = fixings[i]->paymentDate;
        (*cashFlowsProcessed)[i].amount = fixings[i]->fixing;
        // ex-coupon dates are currently assumed to be identical to the payment dates
        (*exCDs)[i] = fixings[i]->paymentDate;
    }
    exCouponDatesProcessed = exCDs;
        
    // if no settlement (as we only introduced it in May 2004), invent one
    if (!settle.get()) {
        settle = SettlementSP(new RollingSettlement);
    }
    isValid = true;
    
    return;
}

void FloatingBond::throwIfNotValid() const
{
    if (isValid == false) {
        throw ModelException("FloatingBond::throwIfNotValid", 
                             "attempting to use invalid bond");
    }
}

void FloatingBond::validateInputs() const
{
    static const string method = "FloatingBond::validateInputs";
    try {

        if (deferredCoupon ) {
            /*
            if (!asset) {
                throw ModelException(method, "Need valid asset to use coupon deferral");
            }
            */
            if ( couponDeferralEndDate.empty()) {
                throw ModelException(method, "No coupon deferral end date has been given");
            }
        }

        if (faceValue < 0.) {
            throw ModelException(method, "face value (" + Format::toString(faceValue) + ") cannot be negative");
        } 

        // validate that there are some cashflows
        if (fixings.empty()){
            throw ModelException(method, "CashFlow array is empty");
        } 
        // validate that cashFlow dates are increasing 
        int i;
        for (i = 1; i < fixings.size(); i++) {
            if (fixings[i-1]->refixDate.isGreater(fixings[i]->refixDate)) {
                throw ModelException(method, "Fixing dates are not increasing(" +
                    fixings[i-1]->refixDate.toString() + " is after " +
                    fixings[i-1]->refixDate.toString() + ")");
            }
        }
        for (i = 0; i < fixings.size(); i++){
            if (fixings[i]->refixDate.isGreater(fixings[i]->paymentDate)) {
                throw ModelException(method, "Fixing dates cannot be after payment dates(" +
                    fixings[i]->refixDate.toString() + " is after " +
                    fixings[i]->paymentDate.toString() + ")");

            }
        }

        if (fixings.size() > 0 && fixings[0]->refixDate != datedDate ) {
            throw ModelException(method, "The first fixing date ("+ fixings[0]->refixDate.toString() +
                ") must be equal to the initial accrual date (" + datedDate.toString() +
                ").");
        }

        // validate that cashFlows are positive
        for (i = 0; i < fixings.size(); i++) {
            if (fixings[i]->refixDate > valueDate) {
                if (Maths::isNegative(fixings[i]->fixing)){
                    throw ModelException(method, "Past fixings must be positive");
                }
            } else {
                if ( (fixings[i]->refixDate > valueDate && !Maths::isPositive(fixings[i]->fixing)) ||
                     (fixings[i]->refixDate > valueDate && Maths::isNegative(fixings[i]->fixing))) {
                    string errorMsg;
                    errorMsg = "A non-positive fixing has been found for a past date:\n" + 
                        fixings[i]->refixDate.toString() + ":  " + Format::toString(fixings[i]->fixing);
                    throw ModelException(method, errorMsg);
                }
            }
        }

        // validate that datedDate is before first cashFlowDate
        if (datedDate.isGreater(fixings[0]->paymentDate)) {
            throw ModelException(method, "datedDate must not be after first cashFlowDate");
        }

        // validate caps and floors
        if ( capFloatingCoupons && floatingCouponCap <= 0.0 ) {
            throw ModelException(method, "floating coupon cap must be strictly positive");
        }

        if ( capFloatingCoupons && floorFloatingCoupons && floatingCouponCap <= floatingCouponFloor ) {
            throw ModelException(method, "floating coupon cap must be strictly greater than floating coupon floor");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


double FloatingBond::priceFromYield(double yield, bool clean, DateTime yDate) const
{
   static const string method = "FloatingBond::priceFromYield";
   try {
       double myPrice = 0.;
       if (yDate > getMaturityDate()) {
           throw ModelException("startDate cannot be after maturityDate");
       }
       CashFlowArraySP cfs = getCashFlows(yDate);
       int i;

       // very simple calc. It works for 30/360 where liborDCC->years will come out to
       // the right fraction, but it will be off a little bit for act/365 and such.
       // It also doesn't do the right thing for long coupons, but you'ld need to take
       // in pseudo coupons for that. As cashflows bonds are rarely used and YTM is not
       // critical, improve only if extremely bored or someone starts complaining.
       // -Tycho
       if (cfs->size() > 1) {

          for (i=cfs->size()-1; i>0; i--) {
                myPrice += (*cfs)[i].amount;
                myPrice /= (1+yield*liborDCC->years((*cfs)[i-1].date, (*cfs)[i].date)); // should this have eom or bad day adjustments
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
               myPrice /= pow((1+yield*liborDCC->years(previousDate, (*cfs)[0].date)),((*cfs)[0].amount-getAccruedAtDate(yDate))/(*cfs)[0].amount); // not right for long front
           } else {
               myPrice /= pow((1+yield*liborDCC->years(previousDate, (*cfs)[0].date)),((*cfs)[0].amount-getRedemption()-getAccruedAtDate(yDate))/((*cfs)[0].amount-getRedemption())); // not right for long front
           }
       } else {
           // a zero -- quote semiannually
           myPrice = (*cfs)[0].amount;
           myPrice /= pow((1+yield/2.), (*liborDCC).years(yDate, (*cfs)[0].date)*2.);
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


Bond* FloatingBond::getBondToPut(const DateTime &putDate, double putLevel) const
{
    static const string method = "FloatingBond::getBondToPut";

    throwIfNotValid();
    try {
        BondCashFlowsSP putBond(new BondCashFlows());

        putBond->faceValue                  = faceValue;
        putBond->datedDate                  = fixings[0]->refixDate;
        putBond->dayCountConvString         = liborDCCString;
        putBond->redemptionPct              = putLevel/faceValue;
        putBond->lastCFIncludesRedemption   = true;
        putBond->cashFlows                  = CashFlowArraySP(new CashFlowArray(0));

        if( fixings.size() == 0 )
        {
            throw ModelException(method,
                "There are no fixings for the floating bond");
        }

        if( putDate < fixings[0]->paymentDate )
        {
            throw ModelException(method, 
                "Next put date (" + putDate.toString() +
                ") is before the first payment date (" + fixings[0]->paymentDate.toString() + ")");
        }

        int i=0;
        while( i < fixings.size() && fixings[i]->paymentDate < putDate ) {
            putBond->cashFlows->push_back((*cashFlowsProcessed)[i]);
            ++i;
        }

        if ( putBond->cashFlows->back().date == putDate) {
            putBond->cashFlows->back().amount += 
                putBond->faceValue * putBond->redemptionPct;
        } else {
            CashFlow redemption(putDate, putLevel);
            putBond->cashFlows->push_back(redemption);
        }

        putBond->validatePop2Object();
        putBond->validateInputs();

        return putBond.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void FloatingBond::exportLiborRates(const string& header)
{
    string cashFlow;
    ErrorHandler::writeMsg(header);

    for (int idx=0 ; idx<cashFlowsProcessed->size() ; ++idx) {
        cashFlow = "Date: " + (*cashFlowsProcessed)[idx].date.toString() +
            "    Amount: " + Format::toString((*cashFlowsProcessed)[idx].amount);
        ErrorHandler::writeMsg(cashFlow);
    }
    ErrorHandler::flush();
}

void FloatingBond::getMarket(const IModel* model, const MarketData* market)
{
    try {
        market->GetReferenceDate(valueDate);
        liborReferenceCurve.getData(model, market);
        settle->getMarket(model, market);

        if ( deferredCoupon == true ) {
            CAsset::getAssetMarketData(model, market, "N", liborReferenceCurve, asset);
        }

        validateInputs();
        calculateLiborRates(false);
    } catch (exception& e) {
        throw ModelException(e, "FloatingBond::getMarket");
    }
}

void FloatingBond::calculateLiborRates(bool isTheta)
{
    // set the floating rates
    DateTime liborEndDate;
    double   libor;
    double   deferral = 0.0;
    if ( !valueDate.empty()) {
        bool     pastDeferral = (valueDate <= couponDeferralEndDate)?false:true;
        for (int idx=0 ; idx<cashFlowsProcessed->size() ; ++idx) {
            if ((fixings[idx]->refixDate > valueDate ||
                ( isTheta && fixings[idx]->refixDate == valueDate )) && !fixings[idx]->isFixed) {
                liborEndDate =  liborPeriod->toDate(fixings[idx]->refixDate);
                if ( isTheta && fixings[idx]->refixDate == valueDate && !Maths::isZero(fixings[idx]->fixing )) {
                    libor = fixings[idx]->fixing;
                } else {
                    libor = liborReferenceCurve->fwd(fixings[idx]->refixDate,
                                                     liborEndDate,
                                                     liborDCC.get(),
                                                     CompoundBasis::SIMPLE);
                }
            } else {
                if ( fixings[idx]->isFixed ) {
                    liborEndDate =  fixings[idx]->paymentDate;
                    libor        = fixings[idx]->fixing;
                } else {
                    liborEndDate =  liborPeriod->toDate(fixings[idx]->refixDate);
                    libor        = fixings[idx]->fixing;
                }
            }

            // only apply libor spread to floating coupons
            if ( !fixings[idx]->isFixed ) {
                libor += liborSpread;
            }

            // negative coupons will not be paid by the bondholder - floor coupon at 0.0
            libor = Maths::max(libor, 0.0);

            // apply caps and floors to floating coupons
            if ( !fixings[idx]->isFixed ) {
                if ( capFloatingCoupons ) {
                    libor = Maths::min(libor, floatingCouponCap);
                }

                if ( floorFloatingCoupons ) {
                    libor = Maths::max(libor, floatingCouponFloor);
                }
            }

            // Cashflow is calculated using year fraction between refixDate and paymentDate (not the
            // Libor period)
            (*cashFlowsProcessed)[idx].amount = 
                libor * faceValue * liborDCC->accruedFactor(fixings[idx]->paymentDate,               
                                                            fixings[idx]->refixDate,
                                                            fixings[idx]->paymentDate,
                                                            false,
                                                            false);

            if ( deferredCoupon == true && !pastDeferral && valueDate < (*cashFlowsProcessed)[idx].date) {
                DateTime paymentDate    = (*cashFlowsProcessed)[idx].date;
                double fwdPrice         = asset->fwdValue(paymentDate);

                // choose how to interpolate the vol - go for traditional route for now
                LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(
                                                           asset->getSpot(), 
                                                           valueDate, 
                                                           paymentDate,
                                                           false));
                // interpolate the vol using our LN request
                CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                // calculate the variance
                double variance = volBS->CalcVar(valueDate, paymentDate);

                double z = (log(Maths::max(DBL_EPSILON,couponDeferralThreshold) / fwdPrice) + .5 * variance) / ( sqrt(variance));
                double prob = N1(z);
                double coupon = 0.0;

                // add deferred portion and defer forward
                if ( ( idx == cashFlowsProcessed->size()-1    && 
                       (*cashFlowsProcessed)[idx].date   <= couponDeferralEndDate )||
                   ( (*cashFlowsProcessed)[idx].date   <= couponDeferralEndDate && 
                   (*cashFlowsProcessed)[idx+1].date >  couponDeferralEndDate )) {
                    coupon       = (*cashFlowsProcessed)[idx].amount + deferral;
                    pastDeferral = true;
                } else {
                    coupon   = (1 - prob) * ((*cashFlowsProcessed)[idx].amount + deferral);
                    deferral =      prob  * ((*cashFlowsProcessed)[idx].amount + deferral);
                }
                (*cashFlowsProcessed)[idx].amount = coupon;
            }
        }
        cashFlowsProcessed->back().amount += redemptionPct*faceValue;
    }
}

// get all cashflows known today
CashFlowArraySP FloatingBond::getKnownCashFlows() const
{
    static const string method = "FloatingBond::getKnownCashFlows";
    try {  
        throwIfNotValid();

        if ( fixings.size() != cashFlowsProcessed->size() ) {
            throw ModelException("FloatingBond::sensShift", "Internal error: number of cash flows different from number of fixings");
        }

        CashFlowArraySP knownCashFlows(new CashFlowArray(0));

        for (int i=0 ; i<fixings.size() ; ++i) {
            if ( (*cashFlowsProcessed)[i].date <= valueDate || 
                 fixings[i]->isFixed                        || 
                 fixings[i]->refixDate <= valueDate         ) {
                knownCashFlows->push_back((*cashFlowsProcessed)[i]);
            }
        }

        return knownCashFlows;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}


bool FloatingBond::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);

        if (BootstrappedYieldCurve::TYPE->isInstance(liborReferenceCurve.get())) {
            IObject* obj = dynamic_cast<IObject*>(liborReferenceCurve.get());
            BootstrappedYieldCurve* yieldCurve = dynamic_cast<BootstrappedYieldCurve*>(obj);
            yieldCurve->sensShift(shift);

            // check whether there is a fixing today
            for (int i=0; i<fixings.size() ;++i) {
                if ( !fixings[i]->isFixed &&
                     fixings[i]->refixDate == valueDate && 
                     Maths::isZero(fixings[i]->fixing) ) {
                    // set the libor fixing
                    DateTime liborEndDate  = liborPeriod->toDate(fixings[i]->refixDate);
                    fixings[i]->fixing = liborReferenceCurve->fwd(fixings[i]->refixDate,
                                                 liborEndDate,
                                                 liborDCC.get(),
                                                 CompoundBasis::SIMPLE);

                }
            }


            calculateLiborRates(true);
        } else {
            throw ModelException("FloatingBond::sensShift", "Reference Libor Curve only supports BootstrappedYieldCurve at the moment.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "FloatingBond::sensShift(theta)");
    }    
    return false; // our components don't have theta type sensitivity
}

DateTime FloatingBond::getAccrualStartDate() const
{
    return datedDate;
}

// when does bond settle?
DateTime FloatingBond::settles(const DateTime& tradeDate) const {
    return settle->settles(tradeDate);
}

class FloatingBondHelper{
public:
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FloatingBond, clazz);
        SUPERCLASS(Bond);
        EMPTY_SHELL_METHOD(defaultFloatingBond);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(faceValue, "face value");
        FIELD(fixings, "the bond's fixings");
        FIELD(datedDate, "date interest begins to accrue");
        FIELD(redemptionPct, "redemption level of the bond as percent of faceValue");
        FIELD_MAKE_OPTIONAL(redemptionPct);
        FIELD(liborReferenceCurve, "identifies the Libor reference curve");
        FIELD(liborPeriod, "the Libor Period");
        FIELD(liborDCCString, "Libor day count convention");
        FIELD(liborSpread, "the spread over Libor");
        // new fields for coupon deferral
        FIELD(deferredCoupon,                    "true for deferred coupon");
        FIELD_MAKE_OPTIONAL(deferredCoupon);
        FIELD(asset,                             "underlying asset for coupon deferral");
        FIELD_MAKE_OPTIONAL(asset);
        FIELD(couponDeferralEndDate,             "last date for coupon deferral decision");
        FIELD_MAKE_OPTIONAL(couponDeferralEndDate);
        FIELD(couponDeferralThreshold,           "thereshold for coupon deferral decision");
        FIELD_MAKE_OPTIONAL(couponDeferralThreshold);
        // new fields for coupon cap/floor
        FIELD(capFloatingCoupons,                "apply cap to floating coupons?");
        FIELD_MAKE_OPTIONAL(capFloatingCoupons);
        FIELD(floatingCouponCap,                 "floating coupon cap");
        FIELD_MAKE_OPTIONAL(floatingCouponCap);
        FIELD(floorFloatingCoupons,              "apply floor to floating coupons?");
        FIELD_MAKE_OPTIONAL(floorFloatingCoupons);
        FIELD(floatingCouponFloor,               "floating coupon floor");
        FIELD_MAKE_OPTIONAL(floatingCouponFloor);
        FIELD(settle, "how bond settles");
        FIELD_MAKE_OPTIONAL(settle);

        // private members that need to be copied when bond is copied
        FIELD(isValid, "private member");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(cashFlowsProcessed, "private member");
        FIELD_MAKE_TRANSIENT(cashFlowsProcessed);
        FIELD(exCouponDatesProcessed, "private member");
        FIELD_MAKE_TRANSIENT(exCouponDatesProcessed);
        FIELD(liborDCC, "private member");
        FIELD_MAKE_TRANSIENT(liborDCC);

        // addin constructor
        Addin::registerConstructor("FLOATING_BOND",
                                   Addin::MARKET,
                                   "Creates a Bond specified via cash flows",
                                   FloatingBond::TYPE);
    }

    static IObject* defaultFloatingBond(){
        return new FloatingBond();
    }
};

CClassConstSP const FloatingBond::TYPE = CClass::registerClassLoadMethod(
    "FloatingBond", typeid(FloatingBond), FloatingBondHelper::load);
bool  FloatingBondLoad() {
    return (FloatingBond::TYPE != 0);
   }



class FloatingBondPaymentHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(FloatingBond::FloatingBondPayment, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPayment);

        FIELD(refixDate,     "the fixing date");
        FIELD(paymentDate,   "the payment date");
        FIELD(fixing,        "the libor fixing");
        FIELD(isFixed,       "true = coupon is fixed, false = coupon is libor + spread");
    }

    static IObject* defaultPayment(){
        return new FloatingBond::FloatingBondPayment();
    }
};


CClassConstSP const FloatingBond::FloatingBondPayment::TYPE = CClass::registerClassLoadMethod(
    "FloatingBondPayment", typeid(FloatingBond::FloatingBondPayment), FloatingBondPaymentHelper::load);

typedef FloatingBond::FloatingBondPaymentArray FloatingBondFloatingBondPaymentArray; // msvc 7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("FloatingBondPaymentArray", FloatingBondFloatingBondPaymentArray);

class FloatBondCreateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    double                      faceValue;
    DateTime                    datedDate;
    double                      redemptionPct;            // Optional
    MaturityPeriodSP            liborPeriod;
    string                      liborDCCString;
    double                      liborSpread;
    YieldCurveWrapper           liborReferenceCurve;
    DateTime                    valueDate;

    DateTimeArraySP             refixDates;             // refix dates
    DateTimeArraySP             paymentDates;           // payment dates
    CDoubleArraySP              fixings;                // past fixings
    CBoolArraySP                isFixed;                // whether the payment is fixed or floating

    // parameters for coupon deferral
    bool                        deferredCoupon;
    CAssetWrapper               asset;
    DateTime                    couponDeferralEndDate;
    double                      couponDeferralThreshold;

    // create a dividend list
    static IObjectSP create(FloatBondCreateAddin *params) {
        static const string method = "FloatBondCreateAddin::create";

        // do some validation on the input parameters
        if ((params->refixDates->size() != params->paymentDates->size())  ||
            (params->refixDates->size() > params->fixings->size())   ||
            (params->refixDates->size() > params->isFixed->size()) )
        {
            throw ModelException(method,
                                 "Mismatch between number of fixing dates, "
                                 "payment dates, fixing levels, and fixed/float flags.");
        }
        FloatingBond::FloatingBondPaymentArray array(0);
        FloatingBond::FloatingBondPaymentArraySP newFloatPayments(new FloatingBond::FloatingBondPaymentArray(array));

        // loop through the paramters creating new dividends
        for (int i = 0; i < params->refixDates->size(); i++) {
            FloatingBond::FloatingBondPaymentSP tmpPayment( new FloatingBond::FloatingBondPayment((*params->refixDates)[i],
                                            (*params->paymentDates)[i],
                                            (*params->fixings)[i],
                                            (*params->isFixed)[i]));
            newFloatPayments->push_back(tmpPayment);
        }

        // assemble the floating bond
        FloatingBondSP floatBond(new FloatingBond(params->faceValue,
                                                  params->datedDate,
                                                  params->redemptionPct,
                                                  params->liborPeriod,
                                                  params->liborDCCString,
                                                  params->liborSpread,
                                                  params->liborReferenceCurve,
                                                  params->valueDate,
                                                  *(newFloatPayments.get())));

        // check that we have created a decent floating bond
        floatBond->validatePop2Object();
        return floatBond;
    }

    // create a dividend list
    static IObjectSP create2(FloatBondCreateAddin *params) {
        static const string method = "FloatBondCreateAddin::create2";

        // do some validation on the input parameters
        if ((params->refixDates->size() != params->paymentDates->size())  ||
            (params->refixDates->size() > params->fixings->size())   ||
            (params->refixDates->size() > params->isFixed->size()) )
        {
            throw ModelException(method,
                                 "Mismatch between number of fixing dates, "
                                 "payment dates, fixing levels, and fixed/float flags.");
        }
        FloatingBond::FloatingBondPaymentArray array(0);
        FloatingBond::FloatingBondPaymentArraySP newFloatPayments(new FloatingBond::FloatingBondPaymentArray(array));

        // loop through the paramters creating new dividends
        for (int i = 0; i < params->refixDates->size(); i++) {
            FloatingBond::FloatingBondPaymentSP tmpPayment( new FloatingBond::FloatingBondPayment((*params->refixDates)[i],
                                            (*params->paymentDates)[i],
                                            (*params->fixings)[i],
                                            (*params->isFixed)[i]));
            newFloatPayments->push_back(tmpPayment);
        }

        // assemble the floating bond
        FloatingBondSP floatBond(new FloatingBond(params->faceValue,
                                                  params->datedDate,
                                                  params->redemptionPct,
                                                  params->liborPeriod,
                                                  params->liborDCCString,
                                                  params->liborSpread,
                                                  params->liborReferenceCurve,
                                                  params->valueDate,
                                                  *(newFloatPayments.get()),
                                                  params->deferredCoupon,
                                                  params->asset,
                                                  params->couponDeferralEndDate,
                                                  params->couponDeferralThreshold));

        // check that we have created a decent floating bond
        floatBond->validatePop2Object();
        return floatBond;
    }

    FloatBondCreateAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FloatBondCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFloatBondCreateAddin);
        // order of registration effects order of parameters in addin function
        FIELD(valueDate,            "value date");
        FIELD(faceValue,            "face value of the floating bond");
        FIELD(datedDate,            "accrual start date");
        FIELD(liborPeriod,                 "Libor Reset Reference");
        FIELD(liborDCCString,       "Floating Rate Day Count Convention");
        FIELD(liborSpread,          "Spread over LIBOR");
        FIELD(liborReferenceCurve,  "underlying yield curve (wrapper) for LIBOR fixings");
        FIELD(refixDates,                  "refix dates");
        FIELD(paymentDates,                "payment dates");
        FIELD(fixings,                     "past fixings (0.0 for future fixings");
        FIELD(isFixed,                     "array of fixed/float flags (true = payment is fixed)");
        FIELD(redemptionPct,        "Redemption as percentage of face value (default = 100%)");
        FIELD_MAKE_OPTIONAL(redemptionPct);
        FIELD(deferredCoupon,               "true for deferred coupon");
        FIELD_MAKE_OPTIONAL(deferredCoupon);
        FIELD(asset,                        "underlying asset for coupon deferral");
        FIELD_MAKE_OPTIONAL(asset);
        FIELD(couponDeferralEndDate,        "last date for coupon deferral decision");
        FIELD_MAKE_OPTIONAL(couponDeferralEndDate);
        FIELD(couponDeferralThreshold,      "thereshold for coupon deferral decision");
        FIELD_MAKE_OPTIONAL(couponDeferralThreshold);

        Addin::registerClassObjectMethod("FLOATING_COUPON_BOND",
                                         Addin::MARKET,
                                         "Creates a handle to a floating bond object",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);

        Addin::registerClassObjectMethod("FLOATING_COUPON_BOND2",
                                         Addin::MARKET,
                                         "Creates a handle to a floating bond object",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create2);
    }
    
    static IObject* defaultFloatBondCreateAddin(){
        return new FloatBondCreateAddin();
    }
 
};

CClassConstSP const FloatBondCreateAddin::TYPE = CClass::registerClassLoadMethod(
    "FloatBondCreateAddin", typeid(FloatBondCreateAddin), load);


DRLIB_END_NAMESPACE

