//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Bond.cpp
//
//   Description : Virtual base class for defining bond cashflows and accrued
//                 interest
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 28, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_BOND_CPP
#include "edginc/Bond.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Maths.hpp"
#include "edginc/BadDayFollowing.hpp"

DRLIB_BEGIN_NAMESPACE

const double Bond::ONE_BASIS_POINT = 0.0001;

// defines an interface to be implemented by concrete classes
/** A Bond defines the cashflows and accrued interest for a bond. */    

typedef struct _TEqYTM
{
    Bond*    bond;
    double   yPrice;
    bool     clean;
    DateTime yDate;
} TEqYTM;

Bond::~Bond() {
    // empty
}

Bond::Bond(CClassConstSP clazz): CObject(clazz){}

/** change the settlement of the bond if applicable */
void Bond::setSettlement(SettlementSP newSettle) {
    // nothing to do
}

static void bondLoad(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Bond, clazz);
    SUPERCLASS(CObject);
}

CClassConstSP const Bond::TYPE = CClass::registerClassLoadMethod(
    "Bond", typeid(Bond), bondLoad);
bool  BondLoad() {
    return (Bond::TYPE != 0);
   }


DateTime Bond::getUnadjMaturityDate() const {
    return getMaturityDate();
}

double Bond::getNotional(const DateTime &date) const {
    return getFaceValue();
}


CashFlowArraySP Bond::getRedemptionPayments() const
{
    throw ModelException("Bond::getRedemptionPayments", "not implemented for this bond type");
}

double Bond::getRedemption(bool asPercent) const
{
    double redemp;
    
    if (asPercent == true) {
        redemp = getRedemption()/getFaceValue();
    } else {
        redemp = getRedemption();
    }

    return redemp;
}

DoubleArraySP Bond::getNotionalsAtDates(DateTimeArrayConstSP dates) const
{
    DoubleArraySP notionals(new DoubleArray(dates->size()));
    for (int i =0; i<notionals->size() ; ++i) {
        (*notionals)[i] = getNotional((*dates)[i]);
    }
    return notionals;
}

DateTimeArraySP Bond::getExCouponDates() const
{
    static const string method = "Bond::getExCouponDates";
    try {  
        CashFlowArraySP cashFlows = getCashFlows();
        DateTimeArraySP exCouponDates(new DateTimeArray);
        int i;
    
        exCouponDates->resize(cashFlows->size());
        for (i=0; i < cashFlows->size(); i++) {
            (*exCouponDates)[i] = (*cashFlows)[i].date;
        }
         
        return exCouponDates;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

double Bond::presentValue(const DateTime &valueDate, YieldCurveConstSP discount) const
{
 
    static const string method = "Bond::presentValue";
    try {  
        int i;
        CashFlowArraySP myCFs = getCashFlows(valueDate);
        double value = 0;
    
        value = 0;
        for (i=0; i<myCFs->size(); i++) {
            value += (*myCFs)[i].amount * discount->pv(valueDate, (*myCFs)[i].date);
        }
    
    //    value -= getAccruedAtDate(valueDate);

            
        return value;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

double Bond::redemptionPV(const DateTime &baseDate, YieldCurveConstSP discount) const
{
    static const string method = "Bond::redemptionPV";
    try {  
        double value = 0;

        if ( baseDate <= getMaturityDate() ) {
            value = getRedemption() * discount->pv(baseDate, getMaturityDate());
        }

        return value;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}



double Bond::duration(const DateTime &valueDate, YieldCurveConstSP discount) const
{
 
    static const string method = "Bond::duration";
    try {  
        int i;
        CashFlowArraySP myCFs = getCashFlows(valueDate);
        double value = 0;
        double duration = 0;
    
        value = 0;
		double time;
        for (i=0; i<myCFs->size(); i++) {
			time = valueDate.yearFrac((*myCFs)[i].date);
            duration += time * (*myCFs)[i].amount * discount->pv(valueDate, (*myCFs)[i].date);
            value += (*myCFs)[i].amount * discount->pv(valueDate, (*myCFs)[i].date);
        }
        if (!Maths::isZero(value)) {
            duration /= value;
        } else {
            duration = 0.0;
        }
    //    value -= getAccruedAtDate(valueDate);

            
        return duration;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

double Bond::convexity(const DateTime &valueDate, YieldCurveConstSP discount) const
{
 
    static const string method = "Bond::convexity";
    try {  
        int i;
        CashFlowArraySP myCFs = getCashFlows(valueDate);
        double value = 0;
        double convexity = 0;
    
        value = 0;
		double time;
        for (i=0; i<myCFs->size(); i++) {
			time = valueDate.yearFrac((*myCFs)[i].date);
            convexity += time*time * (*myCFs)[i].amount * discount->pv(valueDate, (*myCFs)[i].date);
            value += (*myCFs)[i].amount * discount->pv(valueDate, (*myCFs)[i].date);
        }
        if (!Maths::isZero(value)) {
            convexity /= value;
        } else {
            convexity = 0.0;
        }
    //    value -= getAccruedAtDate(valueDate);
            
        return convexity;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

/********************************************************************************
 * Calculates value on fromDate of coupons paid after fromDate and on or before 
 * toDate. Returns zero if fromDate and toDate are equal. 
 ********************************************************************************/
double Bond::couponsPV(const DateTime&   fromDate,
                       const DateTime&   toDate,
                       YieldCurveConstSP discount)
{
    static const string method = "Bond::couponsPV";
    try {
        CashFlowArrayConstSP allCoupons = getCoupons();  
 
        DateTimeArraySP allExDates = getExCouponDates();

        if (allCoupons->size() != allExDates->size()) {
            throw ModelException("Number of ex-coupon dates, " + 
                                    Format::toString(allExDates->size()) +
                                 " does not match number of coupons, "+ 
                                    Format::toString(allCoupons->size()));
        }   

        int i;
        double couponPV = 0.;

        if (allCoupons->size() > 0) {
            for (i=0; i<allCoupons->size(); i++) {
                if ((*allExDates)[i] > fromDate && (*allExDates)[i] <= toDate) {
                    if (!(!discount)) {
                        couponPV += (*allCoupons)[i].amount * discount->pv(fromDate, (*allCoupons)[i].date);
                    } else {
                        couponPV += (*allCoupons)[i].amount;
                    }
                }
            }
        }

        return couponPV;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

// Returns all coupons after startDate
CashFlowArraySP Bond::getCoupons(const DateTime &startDate) const
{
    static const string method = "Bond::getCoupons";
    try {
        CashFlowArraySP myCoupons = getCashFlows(startDate);

        if (myCoupons->size() > 0) {
            myCoupons->back().amount -= getRedemption();
        }
            
        return myCoupons;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Returns all coupons
CashFlowArraySP Bond::getCoupons() const
{
    static const string method = "Bond::getCoupons";
    try {
        CashFlowArraySP myCoupons = getCashFlows();

        if (myCoupons->size() > 0) {
            myCoupons->back().amount -= getRedemption();
        }
            
        return myCoupons;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CashFlowArraySP Bond::getExAdjCashFlows(const DateTime &startDate, YieldCurveConstSP discount) const
{
    static const string method = "Bond::getExAdjCashFlows";
    try {  
        CashFlowArraySP cashFlows = getCashFlows();
        DateTimeArraySP exCouponDates = getExCouponDates();
        CashFlowArraySP exAdjCashFlows(new CashFlowArray);
        int numCoupons;
        int i;

        if (startDate > getMaturityDate()) {          
            throw ModelException(method, "startDate (" + startDate.toString() + 
                    ") must not be greater than maturity (" + getMaturityDate().toString() + ")");
        } else if (startDate >= exCouponDates->back()) {
            // special case. return principal but no coupon
            exAdjCashFlows->resize(1);
            (*exAdjCashFlows)[0].date = getMaturityDate();
            (*exAdjCashFlows)[0].amount = getRedemption();
        } else {
            // count the cashflows after the startDate
            numCoupons = 0;
            while(numCoupons < cashFlows->size() 
                  && (*exCouponDates)[cashFlows->size() - numCoupons - 1] > startDate) {
                numCoupons++;
            }
            if (exCouponDates->back() != cashFlows->back().date) {
                exAdjCashFlows->resize(numCoupons+1); // extra 1 for redemption
                (*exAdjCashFlows)[numCoupons].date = getMaturityDate();
                (*exAdjCashFlows)[numCoupons].amount = getRedemption();
                (*exAdjCashFlows)[numCoupons-1].date = exCouponDates->back();
                (*exAdjCashFlows)[numCoupons-1].amount = (cashFlows->back().amount-getRedemption()) *
                    discount->pv(exCouponDates->back(), cashFlows->back().date);
            } else {
                exAdjCashFlows->resize(numCoupons);
                (*exAdjCashFlows)[numCoupons-1].date = exCouponDates->back();
                (*exAdjCashFlows)[numCoupons-1].amount = (cashFlows->back().amount) *
                    discount->pv(exCouponDates->back(), cashFlows->back().date);
            }
  
            for (i = 1; i<numCoupons; i++) {
                (*exAdjCashFlows)[numCoupons-i-1].date = (*exCouponDates)[cashFlows->size()-i-1];
                (*exAdjCashFlows)[numCoupons-i-1].amount = ((*cashFlows)[cashFlows->size()-i-1].amount) *
                    discount->pv((*exCouponDates)[cashFlows->size()-i-1], (*cashFlows)[cashFlows->size()-i-1].date);
            }
        }
        
        return exAdjCashFlows;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

double YTMFunc(double yield, void *inYTM)
{
    TEqYTM *ytmStruct = (TEqYTM *) inYTM;
    double price;

    price = ytmStruct->bond->priceFromYield(yield, ytmStruct->clean, ytmStruct->yDate);

    return price - ytmStruct->yPrice;
}

double Bond::yieldToMaturity(double yPrice, bool clean, DateTime yDate)
{
    static const string method = "Bond::yieldToMaturity";
    try {
        double yield;
        double priceTolPct = 1.e-10;
        TEqYTM ytmStruct;
        
        if (yDate > getMaturityDate()) {
            throw ModelException("startDate cannot be after maturityDate");
        }
        ytmStruct.bond = this;
        ytmStruct.yPrice = yPrice;
        ytmStruct.clean = clean;
        ytmStruct.yDate = yDate;

        yield = zbrentUseful(
                &YTMFunc,   /* (I) The function to find the root of */
                &ytmStruct,           /* (I) Parameter block */
                -.99,            /* (I) Lowvalue for x */ // be careful that yield/freq != -1
                4.,               /* (I) High value for x */
                yPrice*priceTolPct);           /* (I) Tolerance */

        return yield;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// returns the sensitivity of the dirty price to a 1bp move in YTM
double Bond::pvbp(double mktPrice, bool priceIsClean, DateTime fromDate)
{
    static const string method = "Bond::pvbp";
    try {
        double ytm = yieldToMaturity(mktPrice, priceIsClean, fromDate);
        double loPrice = priceFromYield(ytm - ONE_BASIS_POINT/2.0, 
                                        false,          // dirty
                                        fromDate);     
        double hiPrice = priceFromYield(ytm + ONE_BASIS_POINT/2.0, 
                                        false,          // dirty
                                        fromDate);

        return hiPrice - loPrice;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double Bond::yieldToFirstPut(double yPrice, bool clean, DateTime putDate, double putLevel, DateTime yDate)
{
    static const string method = "Bond::yieldToFirstPut";
    try {

        BondSP putBond = BondSP(getBondToPut(putDate, putLevel));

        return putBond->yieldToMaturity(yPrice, clean, yDate);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double Bond::getBondCouponRate() const
{
   return  0.0;
}

/** get coupon details at base date */
void Bond::getCouponsDetails(const DateTime& baseDate, DateTime& startDate, DateTime& endDate,
                             double& couponRate, double& accruedInterest) const
{
    static const string method = "Bond::getCouponsDetails";
    try {
        DateTime        accrualStartDate;
        // DateTime        accrualStartDate = getAccrualStartDate();
        CashFlowArraySP myCoupons        = getCashFlows();
        if ( myCoupons->size() > 0 ) {
            myCoupons->back().amount -= getRedemption();

            if ( baseDate < accrualStartDate || 
                 baseDate > myCoupons->back().date ) {
                startDate       = DateTime();
                endDate         = DateTime(); 
                couponRate      = getBondCouponRate();
                accruedInterest = 0.0;
            } else {
                int i=0;
                while (i<myCoupons->size() && (*myCoupons)[i].date < baseDate)
                    ++i;

                if ( i == 0 ) {
                    startDate = accrualStartDate;
                } else {
                    startDate = (*myCoupons)[i-1].date;
                }

                endDate         = (*myCoupons)[i].date;
                couponRate      = 0.0;
                accruedInterest = getAccruedAtDate(baseDate);
            }
        } else {
            throw ModelException(method, "No cash flows in bond");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


CashFlowArrayConstSP Bond::getCashFlowRef() const
{
    return getCashFlows();
}

/** returns the coupon stream of the bond from which the fixed leg of the asset swap bond floor is calculated */
DateTimeArraySP Bond::getPaymentDates() const
{
    CashFlowArraySP cashFlows = getCashFlows();

    DateTimeArraySP paymentDates(new DateTimeArray(cashFlows->size()));

    for (int i=0; i<paymentDates->size(); ++i) {
        (*paymentDates)[i] = (*cashFlows)[i].date;
    }

    return paymentDates;
}

// get all cashflows known today - default implementation is all cashflows are known
CashFlowArraySP Bond::getKnownCashFlows() const
{
    return getCashFlows();
}

// return any cashflow due in the next CPND_DATE_WINDOW calendar days
// together with coupon due date and the due date adjusted to next valid business day
AccrualCalendarArraySP Bond::couponDue(const DateTime& fromDate) const {
    static const string method("Bond::couponDue");
    try {
        AccrualCalendarArraySP couponsDue(new AccrualCalendarArray(0));

        // go forward 7 calendar days
        DateTime toDate = fromDate.rollDate(AccrualCalendar::CPND_DATE_WINDOW);

        // Retrieve all known cashflows (redemption included)
        CashFlowArraySP cfl(getKnownCashFlows());

        // Count number of known cashflows falling within window
        int numCFs = 0, i=0, startIdx=0;
        for (i=0; i<cfl->size() && ((*cfl)[i].date < toDate); i++) {
            if ((*cfl)[i].date >= fromDate) {
                if (0==numCFs) {
                    startIdx=i;
                }
                numCFs ++;
            }
        }

        DateTime couponDueDate;
        HolidaySP hols(Holiday::weekendsOnly());
        BadDayFollowing bdf;

        for (i = startIdx; i<startIdx+numCFs; i++) {

            couponDueDate = (*cfl)[i].date;
            AccrualCalendar couponData(couponDueDate, 
                                       bdf.adjust(couponDueDate, hols.get()), 
                                       (*cfl)[i].amount);

            couponsDue->push_back(couponData);
        }
        
        return couponsDue;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return ACCA_DATE_WINDOW calendar days of accrued interest data starting with fromDate 
AccrualCalendarArraySP Bond::accrualCalendar(const DateTime& date) const {
    static const string method("Bond::accrualCalendar");
    try {

        // Accrual calcs are independent of time of day so use SOD
        DateTime fromDate(date.getDate(), DateTime::timeConvert(DateTime::START_OF_DAY));

        AccrualCalendarArraySP accrualData(new AccrualCalendarArray(0));

        // Only concerned with cashflows after passed date
        CashFlowArraySP nextCoupons(getCoupons(fromDate));

        int i, j;
        DateTime accDate;
        bool foundNext = false;
        for (i=0; i<AccrualCalendar::ACCA_DATE_WINDOW; i++) {

            foundNext = false;
            accDate = fromDate.rollDate(i);

            // Need the strictly next cashflow date after accDate
            for (j=0; j<nextCoupons->size(); j++) {
                if (accDate < (*nextCoupons)[j].date) {
                    foundNext = true;
                    break;
                }
            }

            // Populate provided there is a coupon to be paid after accDate and accDate falls after datedDate 
            if (foundNext && (getAccrualStartDate() <= accDate)) { 
                AccrualCalendar accCal(accDate, (*nextCoupons)[j].date, getAccruedAtDate(accDate));
                accrualData->push_back(accCal);
            }
        }

        return accrualData;
        
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

class CFTest: public CObject{
    typedef array<CFTest> CFTestArray;
public:
    static CClassConstSP const TYPE;
    DoubleArray doubArray;
    DateTimeArray dtArray;

    /** for reflection */
    CFTest():  CObject(TYPE){}
    CFTest(int i):  CObject(TYPE) {
        doubArray.resize(i);
        dtArray.resize(i);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CFTest, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCFTest);
        FIELD(dtArray, "date time array")
        FIELD(doubArray, "double array");
        clazz->setPrivate(); // hide this class
    }

    static IObject* defaultCFTest(){
        return new CFTest();
    }
    
};

CashFlowArraySP Bond::getCouponRates()
{
    throw ModelException("Bond::getCouponRates", "This method has not been implemented yet!");
    CashFlowArraySP couponRates;
    return couponRates;
}

CashFlowArraySP Bond::getDayCountFractions()
{
    throw ModelException("Bond::getDayCountFractions", "This method has not been implemented yet!");
    CashFlowArraySP dayCountFractions;
    return dayCountFractions;
}


typedef smartConstPtr<CFTest> CFTestConstSP;
typedef smartPtr<CFTest> CFTestSP;

CClassConstSP const CFTest::TYPE = CClass::registerClassLoadMethod(
    "CFTest", typeid(CFTest), load);

class CFTest2: public CObject{
    typedef array<CFTest2> CFTest2Array;
public:
    static CClassConstSP const TYPE;
    DoubleArray doubArray;
    DateTimeArray dtArray1;
    DateTimeArray dtArray2;

    /** for reflection */
    CFTest2():  CObject(TYPE){}
    CFTest2(int i):  CObject(TYPE) {
        doubArray.resize(i);
        dtArray1.resize(i);
        dtArray2.resize(i);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CFTest2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCFTest2);
        FIELD(dtArray1, "date time array")
        FIELD(doubArray, "double array");
        FIELD(dtArray2, "date time array")
        clazz->setPrivate(); // hide this class
    }

    static IObject* defaultCFTest2(){
        return new CFTest2();
    }
    
};

typedef smartConstPtr<CFTest2> CFTest2ConstSP;
typedef smartPtr<CFTest2> CFTest2SP;

CClassConstSP const CFTest2::TYPE = CClass::registerClassLoadMethod(
    "CFTest2", typeid(CFTest2), load);



class CFAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes two parameters - the yield curve and dates to 
        get pv factor between
    */
    BondSP       bond;
    DateTime        startDate;
    CMarketDataSP   market;

    static IObjectSP addinGetCashFlows(CFAddin* params){
        static const string routine = "CFAddin::addinGetCashFlows";
        try {

            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            BondSP bond(copy(params->bond.get()));

            bond->getMarket(&model, params->market.get());

            /* roll value date forward by 0 days - this is to populate
             any samples which should be set now, but are not
             populated yet, for example when running overnight grids
             for instruments which have a SOD sample */
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(bond);

            CashFlowArraySP flows(bond->getCashFlows(params->startDate));
            CFTestSP output(new CFTest(flows->size()));
            

            for (int i = 0; i < flows->size(); i++) {
                output->dtArray[i] = (*flows)[i].date;
                output->doubArray[i] = (*flows)[i].amount;
            }

            return IObjectSP(output);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP addinGetCashFlowsWithExDiv(CFAddin* params){
        static const string routine = "CFAddin::addinGetCashFlowsWithExDiv";
        try {

            CClosedFormLN model("VolSurface");
            params->bond->getMarket(&model, params->market.get());

            CashFlowArraySP flows(params->bond->getCashFlows());
            DateTimeArraySP exDates(params->bond->getExCouponDates());
            int i;
            // count the flows after the start date
            int nFlows;
            if (params->startDate.empty()) {
                nFlows = flows->size();
            } else {
                if (params->startDate > params->bond->getMaturityDate()) {
                    throw ModelException(routine, "startDate cannot be after maturityDate");
                }

                nFlows = 0;
                while(nFlows < flows->size() 
                    && (*exDates)[flows->size() - nFlows - 1] > params->startDate) {
                    nFlows++;
                }
                if (nFlows < flows->size()) {
                    nFlows++; // show one historic
                }
            }

            CFTest2SP output(new CFTest2(nFlows+1)); // extra one for redemption
            
            
            for (i=nFlows-2; i>=0; i--) {
                output->dtArray1[i] = (*flows)[flows->size() - nFlows + i].date;
                output->doubArray[i] = (*flows)[flows->size() - nFlows + i].amount;
                output->dtArray2[i] = (*exDates)[flows->size() - nFlows + i];
            }
            output->dtArray1[nFlows - 1] = flows->back().date;
            output->doubArray[nFlows - 1] = flows->back().amount - params->bond->getRedemption();
            output->dtArray2[nFlows - 1] = exDates->back();
            output->dtArray1[nFlows] = flows->back().date;
            output->doubArray[nFlows] = params->bond->getRedemption();
            output->dtArray2[nFlows] = flows->back().date;

            return IObjectSP(output);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }
    /** for reflection */
    CFAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CFAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCFAddin);
        FIELD(bond, "bond object");
        FIELD(startDate, "return coupons that fall after this date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(market, "market object");

        Addin::registerClassObjectMethod("GET_COUPONS",
                                         Addin::RISK,
                                         "Returns the cash flows",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)addinGetCashFlows);
        Addin::registerClassObjectMethod("GET_COUPONS_WITH_EX_DIV",
                                         Addin::RISK,
                                         "Returns the cash flows and ex div dates",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)addinGetCashFlowsWithExDiv);
    }

    static IObject* defaultCFAddin(){
        return new CFAddin();
    }
    
};

CClassConstSP const CFAddin::TYPE = CClass::registerClassLoadMethod(
    "CFAddin", typeid(CFAddin), load);



class BondAddin2: public CObject{
    static CClassConstSP const TYPE;

    BondSP bond;
    DateTime aiDate;
    

    static double addinGetAccrued(BondAddin2* params){
        static const string routine = "BondAddin2::addinGetAccrued";
        try {
            return params->bond->getAccruedAtDate(params->aiDate);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    BondAddin2():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BondAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBondAddin2);
        FIELD(bond, "bond object");
        FIELD(aiDate, "return the accrued interest for this date");
        Addin::registerInstanceDoubleMethod("GET_ACCRUED",
                                            Addin::RISK,
                                            "Returns the accrued interest for a bond object",
                                            TYPE,
                                            (Addin::DoubleMethod*)addinGetAccrued);
    }

    static IObject* defaultBondAddin2(){
        return new BondAddin2();
    }
    
};

CClassConstSP const BondAddin2::TYPE = CClass::registerClassLoadMethod(
    "BondAddin2", typeid(BondAddin2), load);

// -------------------------------------------------------------------------

class BondAddin3: public CObject {
    static CClassConstSP const TYPE;

    BondSP          bond;
    DateTime        valueDate;
    double          fairValue;
    bool            fvIsClean;
    CMarketDataSP   market;

    static double addinGetYTM(BondAddin3* params) {
        static const string routine = "BondAddin3::addinGetYTM";
        try {

            CClosedFormLN model("VolSurface");
            CMarketDataSP market;

            if (!!market) {
                market = CMarketDataSP(new MarketData(params->valueDate));
            } else {
                market = params->market;
            }

            params->bond->getMarket(&model, market.get());

            return params->bond->yieldToMaturity(params->fairValue, 
                                                 params->fvIsClean, 
                                                 params->valueDate);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    BondAddin3():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BondAddin3, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBondAddin3);
        FIELD(bond,             "bond object");
        FIELD(valueDate, "value date for which YTM is to be calculated");
        FIELD(fairValue, "fair value of the bond");
        FIELD(fvIsClean, "true if fair value is the clean price");
        FIELD(market,           "market data");

        Addin::registerInstanceDoubleMethod("CVB_GET_YTM",
                                            Addin::CONV_BOND,
                                            "Returns the yield-to-maturity for a bond object",
                                            TYPE,
                                            (Addin::DoubleMethod*)addinGetYTM);
    }

    static IObject* defaultBondAddin3(){
        return new BondAddin3();
    }
    
};

CClassConstSP const BondAddin3::TYPE = CClass::registerClassLoadMethod(
    "BondAddin3", typeid(BondAddin3), load);


class BondAddin4: public CObject {
    static CClassConstSP const TYPE;

    BondSP          bond;

    DateTime        valueDate;
    double          fairValue;
    bool            fvIsClean;

    DateTime        putDate;
    double          putLevel;
    bool            putAdjustForAccrued;
    CMarketDataSP   market;

    static double addinGetYTP(BondAddin4* params) {
        static const string routine = "BondAddin4::addinGetYTP";
        try {
            CClosedFormLN model("VolSurface");
            CMarketDataSP market;

            if (!!market) {
                market = CMarketDataSP(new MarketData(params->valueDate));
            } else {
                market = params->market;
            }

            params->bond->getMarket(&model, market.get());

            double firstPutLevel = params->putLevel;
            if ( params->putAdjustForAccrued ) {
                firstPutLevel += params->bond->getAccruedAtDate(params->valueDate);
            }

            return params->bond->yieldToFirstPut(params->fairValue, params->fvIsClean, 
                                                 params->putDate,   firstPutLevel, 
                                                 params->valueDate);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    BondAddin4():  CObject(TYPE){}

    static void load(CClassSP& clazz){
        REGISTER(BondAddin4, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBondAddin4);
        FIELD(bond,                       "bond object");

        FIELD(valueDate,           "value date for which YTM is to be calculated");
        FIELD(fairValue,           "fair value of the bond");
        FIELD(fvIsClean,           "true if fair value is the clean price");

        FIELD(putDate,             "first put date");
        FIELD(putLevel,            "the (unadjusted) put level");
        FIELD(putAdjustForAccrued, "true to adjust the put level for accrued interest");
        FIELD(market,                     "market object");

        Addin::registerInstanceDoubleMethod("CVB_GET_YTP",
                                            Addin::CONV_BOND,
                                            "Returns the yield-to-first-put for a bond object",
                                            TYPE,
                                            (Addin::DoubleMethod*)addinGetYTP);
    }

    static IObject* defaultBondAddin4(){
        return new BondAddin4();
    }
    
};

CClassConstSP const BondAddin4::TYPE = CClass::registerClassLoadMethod(
    "BondAddin4", typeid(BondAddin4), load);

DRLIB_END_NAMESPACE

