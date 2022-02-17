//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondfloatNtl.cpp
//
//   Description : Bond with notional that floats
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : June 18, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BondFloatNtl.hpp"
#include "edginc/BondParams.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/RollingSettlement.hpp"

DRLIB_BEGIN_NAMESPACE
/** Bond implementation where you input the cash flows */

BondFloatNtl::~BondFloatNtl() {
    // empty
}

BondFloatNtl::BondFloatNtl(): Bond(TYPE), faceValue(100.0) {
    isValid = false;

    // initialize optional parameters
    liborDCCString = "ACT/360";
    redemptionPct = 1.;

    capFloatingCoupons      = false;
    floatingCouponCap       = 0.0;
    floorFloatingCoupons    = false;
    floatingCouponFloor     = 0.0;
}

void BondFloatNtl::fieldsUpdated(const CFieldArray& fields){
    initialize();
}

/** Change the settlement of the bond if applicable */
void BondFloatNtl::setSettlement(SettlementSP newSettle) {
    settle = SettlementSP(newSettle.get());
}

DateTime BondFloatNtl::getMaturityDate() const
{
    throwIfNotValid();
    return maturityDate;
}

double BondFloatNtl::getFaceValue() const {
    return faceValue;
}

double BondFloatNtl::getRedemption() const
{
    throwIfNotValid();
    return redemptionPct*notionalsProcessed->back();
}

double BondFloatNtl::getNotional(const DateTime &ntlDate) const
{
    throwIfNotValid();
    
    double notional=0;
    int i;

    if (ntlDate >= maturityDate) {
        notional = notionalsProcessed->back();
    } else if (ntlDate < (*fixingDates)[0]) {
        notional = 0; // shouldn't matter what we set it to here
    } else {
        for (i = 1; i < fixingDates->size(); i++) {
            if (ntlDate < ((*fixingDates)[i])){
                notional = (*notionalsProcessed)[i-1] * 
                    (1 + (*accretionRates)[i-1] * liborDCC->years((*fixingDates)[i-1], ntlDate));
                break;
            }
        }
    }

    return notional;
}


CashFlowArraySP BondFloatNtl::getCashFlows() const
{
    static const string method = "BondFloatNtl::getCashFlows";
    try {  
        throwIfNotValid();
        CashFlowArraySP myCashFlows(new CashFlowArray());

        myCashFlows->resize(1);
        (*myCashFlows)[0].date = getMaturityDate();
        // why not scaled by redemptionPct??
        (*myCashFlows)[0].amount = notionalsProcessed->back();

        return myCashFlows;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

CashFlowArraySP BondFloatNtl::getCashFlows(const DateTime &startDate) const
{
    static const string method = "BondFloatNtl::getCashFlows_DateTime";

    try { 
        throwIfNotValid();

        if (startDate > getMaturityDate()) {
            CashFlowArraySP noCashFlows(new CashFlowArray(0));
            return noCashFlows;
        }
        else {
            return getCashFlows();
        }

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double BondFloatNtl::getAccruedAtDate(const DateTime &aiDate) const 
{
    throwIfNotValid();
    
    return 0.;
}

DateTimeArraySP BondFloatNtl::getExCouponDates() const
{
    static const string method = "BondFloatNtl::getExCouponDates";
    try {
        
        DateTimeArraySP exDates(new DateTimeArray(0));
        return exDates;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CashFlowArraySP BondFloatNtl::getCoupons() const
{
    static const string method = "BondFloatNtl::getCoupons";
    try {
        
        CashFlowArraySP myCoupons(new CashFlowArray(0));
        return myCoupons;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CashFlowArraySP BondFloatNtl::getCoupons(const DateTime &startDate) const
{
    static const string method = "BondFloatNtl::getCoupons";
    try {
        
        if (startDate > getMaturityDate()) {
            CashFlowArraySP noCoupons(new CashFlowArray(0));
            return noCoupons;
        }
        else {
            return getCoupons();
        }

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


CashFlowArraySP BondFloatNtl::getExAdjCashFlows(const DateTime &startDate, YieldCurveConstSP discount) const
{
    static const string method = "BondFloatNtl::getExAdjCashFlows";
    try {  

        // There are no coupons or ex-dates so return unadjusted maturity cashflow
        CashFlowArraySP cashFlows = getCashFlows();
        
        return cashFlows;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}
    
double BondFloatNtl::getMaturityCoupon() const
{
    return 0.;
}

void BondFloatNtl::validatePop2Object()
{
    if ( fixingDates->size() <= 1 ) {
        throw ModelException("BondFloatNtl::validatePop2Object", 
                             "Need at least two fixing dates for floating notional bonds");
    }
    // convert the day count convention string to an object
    liborDCC = DayCountConventionSP(DayCountConventionFactory::make(liborDCCString));
    if (datedDate.getDate() == 0)
        datedDate          = (*fixingDates)[0];
    else if (datedDate != (*fixingDates)[0])
        throw ModelException("BondFloatNtl::validatePop2Object", 
                             "dated date should be equal to the first fixing date");

    // validate caps and floors
    if ( capFloatingCoupons && floatingCouponCap <= 0.0 ) {
        throw ModelException("BondFloatNtl::validatePop2Object", 
                             "floating coupon cap must be strictly positive");
    }

    if ( capFloatingCoupons && floorFloatingCoupons && floatingCouponCap <= floatingCouponFloor ) {
        throw ModelException("BondFloatNtl::validatePop2Object", 
                             "floating coupon cap must be strictly greater than floating coupon floor");
    }

    // if no settlement (as we only introduced it in May 2004), invent one
    if (!settle.get()) {
        settle = SettlementSP(new RollingSettlement);
    }
    return;
}

void BondFloatNtl::throwIfNotValid() const
{
    if (isValid == false) {
        throw ModelException("BondFloatNtl::throwIfNotValid", 
                             "attempting to use invalid bond");
    }
}

void BondFloatNtl::validateInputs() const
{
    static const string method = "BondFloatNtl::validateInputs";
    try {
        int i;

        if (faceValue < 0.) {
            throw ModelException(method, "face value (" + Format::toString(faceValue) + ") cannot be negative");
        } 

        // validate that fixing dates are increasing
        for (i = 1; i < fixingDates->size(); i++) {
            if ((*fixingDates)[i-1].isGreater((*fixingDates)[i])){
                throw ModelException(method, "fixingDates dates are not increasing");
            }
        }

        // this could be relaxed to price before issue
        if ((*fixingDates)[0].isGreater(valueDate)) {
            throw ModelException(method, "cannot price before first fixing date");
        }
        
        if (fixingDates->back() != maturityDate) {
            throw ModelException(method, "there must be a fixing date on the maturity date");
        }
        
        // validate that the number of fixing dates, rates, and notionals are the same
        if (fixingDates->size() != rates->size()) {
            throw ModelException(method, "fixingDates is not the same length as rates");
        }

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void BondFloatNtl::getMarket(const IModel* model, const MarketData* market)
{
    market->GetReferenceDate(valueDate);
    liborReferenceCurve.getData(model, market);
    settle->getMarket(model, market);

    validateInputs();
    initialize();
    isValid = true;
}

void BondFloatNtl::initialize() 
{
    int i;
    int currentIdx;

    ratesProcessed = DoubleArraySP(rates.clone());
    accretionRates = DoubleArraySP(new DoubleArray());
    accretionRates->resize(rates->size());
    notionalsProcessed = DoubleArraySP(new DoubleArray(rates->size()));
    
    // treat maturity specially
    if (fixingDates->back() == valueDate) {
        // to do
    }

    // figure out where we are today
    for (i = 1; i <fixingDates->size(); i++) {
        if ((*fixingDates)[i].isGreater(valueDate)){
            currentIdx = i-1;
            break;
        }
    }
    
    // set the future rates
    calculateLiborRates();

    // set the accretion rates
    for (i = 0; i <fixingDates->size()-1; i++) {
        // floor at zero, unless there's an explicit floor
        if ( !floorFloatingCoupons) {
            (*accretionRates)[i] = Maths::max((*ratesProcessed)[i] + liborSpread, 0.);
        } else {
            (*accretionRates)[i] = Maths::max((*ratesProcessed)[i] + liborSpread, floatingCouponFloor);
        }

        // cap, if applicable
        if ( capFloatingCoupons ) {
            (*accretionRates)[i] = Maths::min((*accretionRates)[i] , floatingCouponCap);
        }
    }

    (*notionalsProcessed)[0] = faceValue;

    // set the notionals
    for (i = 0; i <notionalsProcessed->size()-1; i++) {
         (*notionalsProcessed)[i+1] = (*notionalsProcessed)[i]*
             (1+(*accretionRates)[i]*liborDCC->years((*fixingDates)[i], (*fixingDates)[i+1]));
    }
}

void BondFloatNtl::calculateLiborRates()
{
    // set the floating rates
    DateTime liborEndDate;
    double   libor;
    if ( !valueDate.empty()) {
        for (int idx=0 ; idx<fixingDates->size()-1 ; ++idx) {
            if ((*fixingDates)[idx] > valueDate ||
                (*fixingDates)[idx] == valueDate ) {
                liborEndDate =  liborPeriod->toDate((*fixingDates)[idx]);
                if ((*fixingDates)[idx] == valueDate && !Maths::isZero((*rates)[idx])) {
                    libor = (*rates)[idx];
                } 
                else {
                    libor = liborReferenceCurve->fwd((*fixingDates)[idx],
                                                     liborEndDate,
                                                     liborDCC.get(),
                                                     CompoundBasis::SIMPLE);
                }
                
                (*ratesProcessed)[idx] = libor;
            }
        }
    }
}

// get all cashflows known today - default implementation is all cashflows are known
CashFlowArraySP BondFloatNtl::getKnownCashFlows() const
{
    CashFlowArraySP knownCashflows;
    if ( valueDate >= maturityDate || 
         (fixingDates->size() > 1 && valueDate >= fixingDates->back())) {
        knownCashflows = getCashFlows();
    } else {
        knownCashflows = CashFlowArraySP(new CashFlowArray(0));
    }
    return knownCashflows;

}

bool BondFloatNtl::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
        
        if (BootstrappedYieldCurve::TYPE->isInstance(liborReferenceCurve.get())) {
            IObject* obj = dynamic_cast<IObject*>(liborReferenceCurve.get());
            BootstrappedYieldCurve* yieldCurve = dynamic_cast<BootstrappedYieldCurve*>(obj);
            yieldCurve->sensShift(shift);
            
            // check whether there is a fixing today
            for (int i=0; i<rates->size() ;++i) {
                if ( (*fixingDates)[i] == valueDate && 
                     Maths::isZero((*rates)[i])) {
                    // set the libor fixing
                    DateTime liborEndDate  = liborPeriod->toDate((*fixingDates)[i]);
                    (*rates)[i] = liborReferenceCurve->fwd((*fixingDates)[i],
                                                           liborEndDate,
                                                           liborDCC.get(),
                                                           CompoundBasis::SIMPLE);

                }
            }
            
            initialize();
        } else {
            throw ModelException("BondFloatNtl::sensShift", "Reference Libor Curve only supports BootstrappedYieldCurve at the moment.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "BondFloatNtl::sensShift(theta)");
    }    
    return false; // our components don't have theta type sensitivity
}



double BondFloatNtl::priceFromYield(double yield, bool clean, DateTime yDate) const
{
   static const string method = "BondFloatNtl::priceFromYield";
   try {
       double myPrice = 0.;
       if (yDate > getMaturityDate()) {
           throw ModelException("startDate cannot be after maturityDate");
       }
       CashFlowArraySP cfs = getCashFlows(yDate);

       // a zero -- quote semiannually for now. Should probably quote on libor period
       myPrice = (*cfs)[0].amount;
       myPrice /= pow((1+yield/2.), (*liborDCC).years(yDate, (*cfs)[0].date)*2.);
 

       if (clean == true) {
           myPrice -= getAccruedAtDate(yDate);
       }

       return myPrice;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

Bond* BondFloatNtl::getBondToPut(const DateTime &putDate, double putLevel) const
{
   throwIfNotValid();
    static const string method = "BondFloatNtl::getBondToPut";
    try {
        BondParamsSP putBond(new BondParams());

        putBond->faceValue          = putLevel;
        putBond->redemptionPct      = redemptionPct;
        putBond->couponPct          = 0;
        putBond->frequency          = 2;
        putBond->maturityDate       = putDate;
        putBond->datedDate          = (*fixingDates)[0];

        putBond->validatePop2Object();
        putBond->initialize();

        return putBond.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DateTime BondFloatNtl::getAccrualStartDate() const
{
    return (*fixingDates)[0];
}

// when does bond settle?
DateTime BondFloatNtl::settles(const DateTime& tradeDate) const {
    return settle->settles(tradeDate);
}

class BondFloatNtlHelper{
public:
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BondFloatNtl, clazz);
        SUPERCLASS(Bond);
        EMPTY_SHELL_METHOD(defaultBondFloatNtl);
        // IMPLEMENTS(ITweakableWithRespectTo<RateParallel>);
        //IMPLEMENTS(ITweakableWithRespectTo<RatePointwise>);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(maturityDate, "maturity date");
        FIELD(faceValue, "face value");
        FIELD(fixingDates, "the Libor fixing dates");
        FIELD(rates, "the Libor rates on the fixing dates");
        FIELD(liborPeriod, "the Libor Period");
        FIELD(liborDCCString, "Libor day count convention");
        FIELD(liborSpread, "the spread over Libor");
        FIELD(liborReferenceCurve, "identifies the Libor reference curve");
        FIELD(redemptionPct, "redemption level of the bond as percent of floating notional");
        FIELD_MAKE_OPTIONAL(redemptionPct);
        FIELD(datedDate, "date notional begins to accrue");
        FIELD_MAKE_OPTIONAL(datedDate);

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
        FIELD(ratesProcessed, "private member");
        FIELD_MAKE_TRANSIENT(ratesProcessed);
        FIELD(accretionRates, "private member");
        FIELD_MAKE_TRANSIENT(accretionRates);
        FIELD(notionalsProcessed, "private member");
        FIELD_MAKE_TRANSIENT(notionalsProcessed);
        FIELD(liborDCC, "private member");
        FIELD_MAKE_TRANSIENT(liborDCC);
    }

    static IObject* defaultBondFloatNtl(){
        return new BondFloatNtl();
    }
};

CClassConstSP const BondFloatNtl::TYPE = CClass::registerClassLoadMethod(
    "BondFloatNtl", typeid(BondFloatNtl), BondFloatNtlHelper::load);
bool  BondFloatNtlLoad() {
    return (BondFloatNtl::TYPE != 0);
   }


class FloatNtlAccreted: public CObject {
    static CClassConstSP const TYPE;

    // addin takes two parameters - the yield curve and dates to 
    //  get pv factor between
    BondFloatNtlSP      bond;
    DateTime            startDate;
    CMarketDataSP       market;

    static IObjectSP floatNtlAccreted(FloatNtlAccreted* params){
        static const string routine = "FloatNtlAccreted::floatNtlAccreted";
        try {

            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            BondFloatNtlSP bond(copy(params->bond.get()));

            bond->getMarket(&model, params->market.get());

            /* roll value date forward by 0 days - this is to populate
             any samples which should be set now, but are not
             populated yet, for example when running overnight grids
             for instruments which have a SOD sample */
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(bond);

            bond->initialize();

            CashFlowArraySP flows(new CashFlowArray(bond->notionalsProcessed->size()));

            // set the notionals
            for (int i = 0; i <bond->notionalsProcessed->size(); i++) {
                (*flows)[i].date   = (*(*(bond.get())).fixingDates.get())[i];
                (*flows)[i].amount = (*(*(bond.get())).notionalsProcessed.get())[i];
            }

            return IObjectSP(flows);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    FloatNtlAccreted():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FloatNtlAccreted, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFloatNtlAccreted);
        FIELD(bond, "bond object");
        FIELD(startDate, "return coupons that fall after this date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(market, "market object");

        Addin::registerClassObjectMethod("GET_FLOAT_NTL_ACCRETED",
                                         Addin::RISK,
                                         "Returns the accreted values on fixing dates for floating notional bonds",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)floatNtlAccreted);
    }

    static IObject* defaultFloatNtlAccreted(){
        return new FloatNtlAccreted();
    }
    
};

CClassConstSP const FloatNtlAccreted::TYPE = CClass::registerClassLoadMethod(
    "FloatNtlAccreted", typeid(FloatNtlAccreted), load);

DRLIB_END_NAMESPACE

