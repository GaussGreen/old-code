//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CreditLegConvention.cpp
//
//   Description : Base class representing "credit legs" conventions (see below)
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditLegConvention.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/FixedCashFlow.hpp"
#include "edginc/AdhocCashFlow.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/FullCreditFeeLeg.hpp"
#include "edginc/CDOContingentLeg.hpp"

DRLIB_BEGIN_NAMESPACE

/** Returns the name of this object */
string CreditLegConvention::getName() const
{
    if(name != "")
    {
        return name;
    }
    else
    {
        throw ModelException("Name is empty");
    }
    
}

/**
 * Returns the default start date given the trade date
 * [Implements ICreditLegConvention]
 * */
DateTime CreditLegConvention::startDate(const DateTime& tradeDate) const {
    return tradeDate.rollDate(startDateOffset);
}

/**
 * Generates the fee leg given the end date
 * [Implements ICreditLegConvention]
 * */
ICreditFeeLegSP CreditLegConvention::generateFeeLeg(
    const DateTime& startDate,
    const DateTime& endDate,
    double coupon,
    double upfrontPayment,
    double notional) const
{
    const string method = "CreditLegConvention::generateFeeLeg";
    int i;
    
    // 1. Build "theoretical dates" i.e. non adjusted
    //    accrual dates, starting at "startDate" and ending at
    //    "endDate"
    
    // previousAccrualStart can be used in case of a 
    // stub payment of type BOND or NONE
    DateTime previousAccrualStart = startDate;
    
    // decompose payment period
    int count;
    string interval;
    paymentPeriod->decompose(count, interval);
    
    DateTimeArray theoreticalDates;
    if (stubType == FRONT) {
        int nbDates = 0;
        DateTime currentDate;
        for (currentDate = endDate; currentDate > startDate;) {
            // using "push_back" means that theoreticalDates will
            // be sorted by decreasing dates (end date first)
    		theoreticalDates.push_back(currentDate);
            nbDates++;
            currentDate = MaturityPeriod::toDate(-count * nbDates, interval, endDate);
		}
        if (nbDates == 0) {
            throw ModelException(method,
                "Unable to generate fee leg: start date >= end date");
        }
        
        // update previousAccrualStart
        if (currentDate != startDate) {
            previousAccrualStart = paymentBDCSP->adjust(currentDate, holidays.get());
        }

        // retrieve stub length ("long" / "short" / "trigger")
        string stLength = stubLength;
        if (stLength == TRIGGER) {
            // use stubLengthTrigger
            DateTime triggerDate = stubLengthTrigger->toDate(startDate);
            if (theoreticalDates.back() < triggerDate) {
                stLength = LONG;
            } else {
                stLength = SHORT;
            }
        }
        
        // remove first date (in time) if stub length is LONG, unless
        // we have only one date in theoreticalDates
        if (stLength == LONG && nbDates > 1) {
            theoreticalDates.pop_back();
        }
        
        // add start date (so after this step, theoreticalDates will 
        // have more than 2 elements)
        theoreticalDates.push_back(startDate);
        
        // put theoreticalDates in increasing order i.e. revert the array
        int n = theoreticalDates.size();
        DateTime tempDate;
        for (i = 0; i < (n / 2); ++i) {
			tempDate = theoreticalDates[i];
            theoreticalDates[i] = theoreticalDates[n - 1 - i];
            theoreticalDates[n - 1 - i] = tempDate;
		}
    } else if (stubType == BACK) {
        int nbDates = 0;
        for (DateTime currentDate = startDate; currentDate < endDate;) {
            // using "push_back" means that theoreticalDates will
            // be sorted by increasing dates (start date first)
            theoreticalDates.push_back(currentDate);
            nbDates++;
            currentDate = MaturityPeriod::toDate(count * nbDates, interval, startDate);
        }
        if (nbDates == 0) {
            throw ModelException(method,
                "Unable to generate fee leg: start date >= end date");
        }

        // retrieve stub length ("long" / "short")
        string stLength = stubLength;
        if (stLength == TRIGGER) {
            // use stubLengthTrigger
            int countSLT;
            string intervalSLT;
            stubLengthTrigger->decompose(countSLT, intervalSLT);
            DateTime triggerDate = MaturityPeriod::toDate(-countSLT, intervalSLT, endDate);
            if (theoreticalDates.back() > triggerDate) {
                stLength = LONG;
            } else {
                stLength = SHORT;
            }
        }
        
        // remove last date (in time) if stub length is LONG, unless
        // we have only one date in theoreticalDates
        if (stLength == LONG && nbDates > 1) {
            theoreticalDates.pop_back();
        }
        
        // add end date (after this step, theoreticalDates will 
        // have more than 2 elements and will be sorted in increasing order)
        theoreticalDates.push_back(endDate);
    } else {
        throw ModelException(method,
            "Stub type must be " + FRONT + " or " + BACK);
    }
    
    // 2. Build accrual dates
    
    DateTimeArray accrualDates(theoreticalDates.size() - 1);
    // use paymentDCC to adjust accrual dates (except last one which
    // remains unadjusted)
    for (i = 0; i < accrualDates.size() - 1; ++i) {
		accrualDates[i] = paymentBDCSP->adjust(theoreticalDates[i+1], holidays.get());
	}
    accrualDates[accrualDates.size() - 1] = theoreticalDates[accrualDates.size()];
    
    // 3. Build payment dates - same as accrual dates except for the last
    //    one which is adjusted using lastPaymentBDC
    
    DateTimeArray paymentDates(accrualDates);
    paymentDates[paymentDates.size() - 1] =
        lastPaymentBDCSP->adjust(accrualDates[paymentDates.size() - 1], holidays.get());
    
    // 4. Build cash flows
    
    // take "stubPayment" into account
    double upfrontForBondStub = 0.0;
    if (stubPayment == BOND) {
        upfrontForBondStub =
            coupon * paymentDCCSP->years(previousAccrualStart, startDate);
    }
    
    AbstractCashFlowArraySP cashFlows(
        new AbstractCashFlowArray(paymentDates.size() + 1));

    // upfront cash flow
    (*cashFlows)[0] = AdhocCashFlowSP(new AdhocCashFlow(
        DateTime(),                             // value date - irrelevant for an AdhocCashFlow
        startDate,                              // pay date
        notional * (upfrontPayment - upfrontForBondStub)));  // amount
    
    // "coupon" cash flows
    FixedCashFlowSP fixedCashFlow;
    DateTime prevDate = startDate;
    if (stubPayment == BOND || stubPayment == NONE) {
        prevDate = previousAccrualStart;   
    }
    CouponNotionalTypes couponNotionalType;
    if (notionalObservationType == OBS_START || notionalObservationType == OBS_END) {
        couponNotionalType = OBSERVATION_DATE;
    } else if (notionalObservationType == OBS_MIDDLE) {
        couponNotionalType = AVERAGE;
    } else {
        throw ModelException(method, "Cannot assign coupon notional type");
    }

    DateTimeSP observationDate;
    for (i = 0; i < paymentDates.size(); ++i) {
        fixedCashFlow = FixedCashFlowSP(new FixedCashFlow(
            DateTime(),     // value date - irrelevant for a FixedCashFlow
            paymentDates[i],// pay date
            notional,       // notional
            prevDate,       // accrual start
            accrualDates[i],// accrual end
            paymentDCCSP,   // DCC
            coupon));       // fixing

        if (notionalObservationType == OBS_START) {
            observationDate.reset(new DateTime(prevDate));
        } else if (notionalObservationType == OBS_END) {
            observationDate.reset(new DateTime(accrualDates[i]));
        }
        (*cashFlows)[i+1] = CreditCashFlowSP(new CreditCashFlow(
            fixedCashFlow,        // risk free cash flow
            couponNotionalType,   // coupon notional type
            observationDate));    // observation date
    
        prevDate = accrualDates[i];
    }

    // 5. Create fee leg
    return ICreditFeeLegSP(new FullCreditFeeLeg(cashFlows));
}

/**
 * Generates the contingent leg given the end date
 * [Implements ICreditLegConvention]
 * */
ICreditContingentLegSP CreditLegConvention::generateContingentLeg(
    const DateTime& startDate,
    const DateTime& endDate,
    double coupon,
    double upfrontPayment,
    double notional) const
{
    DateTimeArray obsStartDates(1);
    DateTimeArray obsEndDates(1);
    DateTimeArray payDates(1);
    DoubleArray notionals(1);
    
    obsStartDates[0] = startDate;
    obsEndDates[0] = endDate;
    payDates[0] = lastPaymentBDCSP->adjust(endDate, holidays.get());
    notionals[0] = notional;
    
    return ICreditContingentLegSP(new CDOContingentLeg(
        defaultPaymentDelayDays,
        obsStartDates,
        obsEndDates,
        payDates,
        notionals,
        paymentUponDefault));
}

/** Virtual destructor */
CreditLegConvention::~CreditLegConvention() {}

/** Constructor (used by reflection) */
CreditLegConvention::CreditLegConvention(CClassConstSP clazz):
    MarketObject(clazz),
    name(""),
    startDateOffset(1),
    paymentUponDefault(true),
    defaultPaymentDelayDays(7),
    paymentDCC("Act/360"),
    paymentDCCSP(   ),
    paymentBDC("M"),
    paymentBDCSP(   ),
    lastPaymentBDC("F"),
    lastPaymentBDCSP(   ),
    paymentPeriod(new MaturityPeriod("3M")),
    notionalObservationType(OBS_MIDDLE),
    recoverNotional(false),
    stubType(FRONT),
    stubLength(TRIGGER),
    stubLengthTrigger(new MaturityPeriod("1M")),
    stubPayment(SIMPLE) {}

/** Default constructor */
IObject* CreditLegConvention::defaultConstructor() {
    return new CreditLegConvention(TYPE);
}    

/** Checks parameters immediately after object is constructed */
void CreditLegConvention::validatePop2Object() {
    const string method = "CreditLegConvention::validatePop2Object";
    try {
        // Check stub type: "FRONT" or "BACK"
        if (stubType != FRONT && stubType != BACK) {
            throw ModelException(method,
                "Stub type must be " + FRONT + " or " + BACK);
        }

        // Check stub length: "LONG", "SHORT" or "TRIGGER"
        if (stubLength != LONG && stubLength != SHORT && stubLength != TRIGGER) {
            throw ModelException(method,
                "Stub length must be " + LONG + ", " + SHORT + " or " + TRIGGER);
        }
    
        // Check stub payment: "BOND", "SIMPLE" or "NONE"
        if (stubPayment != BOND && stubPayment != SIMPLE && stubPayment != NONE) {
            throw ModelException(method,
                "Stub payment must be " + BOND + ", " + SIMPLE + " or " + NONE);
        }
        
        // Check payment period (payment period with length 0 could
        // create infinite loops
        int count;
        string interval;
        paymentPeriod->decompose(count, interval);
        if (count <= 0) {
            throw ModelException(method,
                "Invalid payment period (length should be > 0): " + paymentPeriod->toString());
        }
        
        // Check notionalObservationType
        if (notionalObservationType != OBS_START &&
            notionalObservationType != OBS_END &&
            notionalObservationType != OBS_MIDDLE)
        {
            throw ModelException(method,
                "Notional observation type must be " +
                OBS_START + ", " + OBS_END + " or " + OBS_MIDDLE);
        }
        
        // Populates day count and bad day conventions
        paymentDCCSP.reset(DayCountConventionFactory::make(paymentDCC)); 
        paymentBDCSP.reset(BadDayConventionFactory::make(paymentBDC));
        lastPaymentBDCSP.reset(BadDayConventionFactory::make(lastPaymentBDC));
        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Populate from market cache */
void CreditLegConvention::getMarket(const IModel* model, const MarketData* market) {
    try {
        holidays.getData(model, market);
        discount.getData(model, market);
    } catch (exception& e){
        throw ModelException(e, "CreditLegConvention::getMarket");
    }
}

/**
 * Access to discount curve name
 * [Implements ICreditLegConvention]
 * */
string CreditLegConvention::getDiscountName() const {
    return discount.getName();
}

/**
 * Access to discount curve
 * [Implements ICreditLegConvention]
 * */
YieldCurveConstSP CreditLegConvention::getDiscount() const {
    return discount.getSP();
}

/** access to recoverNotional flag */
bool CreditLegConvention::getRecoverNotional() const
{
    return recoverNotional;
}

/** Invoked when Class is 'loaded' */
void CreditLegConvention::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditLegConvention, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(ICreditLegConvention);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name, "object name");
    FIELD_MAKE_OPTIONAL(name);
    // contingent leg fields
    FIELD(startDateOffset,
        "Number of calendar days to add to "
        "the trade date to get the protection start date [default is 1]");
    FIELD_MAKE_OPTIONAL(startDateOffset);
    FIELD(paymentUponDefault,
        "Pay upon default (TRUE) or at payment date (FALSE) [default is TRUE]");
    FIELD_MAKE_OPTIONAL(paymentUponDefault);
    FIELD(defaultPaymentDelayDays,
        "Number of delay days for default payment when "
        "paymentUponDefault = TRUE [default is 7] ");
    FIELD_MAKE_OPTIONAL(defaultPaymentDelayDays);

    // fee leg fields
    FIELD(paymentDCC,
        "Payment day count convention [default is Act/360]");
    FIELD_MAKE_OPTIONAL(paymentDCC);
    FIELD_NO_DESC(paymentDCCSP);
    FIELD_MAKE_TRANSIENT(paymentDCCSP);
    FIELD(paymentBDC,
        "Bad day convention used to adjust all payment dates, "
        "except the last one [default is 'Modified']");
    FIELD_MAKE_OPTIONAL(paymentBDC);
    FIELD_NO_DESC(paymentBDCSP);
    FIELD_MAKE_TRANSIENT(paymentBDCSP);
    FIELD(lastPaymentBDC,
        "Bad day convention used to adjust the last payment date "
        "[default is 'Following']");
    FIELD_MAKE_OPTIONAL(lastPaymentBDC);
    FIELD_NO_DESC(lastPaymentBDCSP);
    FIELD_MAKE_TRANSIENT(lastPaymentBDCSP);
    FIELD(paymentPeriod, "Payment period [default is 3M]");
    FIELD_MAKE_OPTIONAL(paymentPeriod);
    FIELD(holidays, "Holidays");
    FIELD(notionalObservationType,
        "Describes when we observe the coupon notional over an "
        "accrual fee period [T0, T1]: " + OBS_START + "=T0, " +
        OBS_END + "=T1, " + OBS_MIDDLE + "=(T0+T1)/2 [default is " +
        OBS_MIDDLE + "]");
    FIELD_MAKE_OPTIONAL(notionalObservationType);

    FIELD(recoverNotional, "Is tranche reduced by recovered notional (index convention). [default = false]");
    FIELD_MAKE_OPTIONAL(recoverNotional);

    // stub fields
    FIELD(stubType,
        "Stub type (" + FRONT + " or " + BACK +
        ") [default is " + FRONT + "]");
    FIELD_MAKE_OPTIONAL(stubType);
    FIELD(stubLength,
        "Stub length (" + LONG + ", " + SHORT + " or " + TRIGGER +
        ") [default is " + TRIGGER + " i.e. use stubLengthTrigger]");
    FIELD_MAKE_OPTIONAL(stubLength);
    FIELD(stubLengthTrigger,
        "Stub length trigger [default is 1M]");
    FIELD_MAKE_OPTIONAL(stubLengthTrigger);    
    FIELD(stubPayment,
        "Stub payment (" + BOND + ", " + SIMPLE + " or " + NONE +
        ") [default is " + SIMPLE + "]");
    FIELD_MAKE_OPTIONAL(stubPayment);

    // discount curve
    FIELD(discount, "Discount curve");
}

// strings describing when we observe the coupon notional
const string CreditLegConvention::OBS_START  = "OBS_START";
const string CreditLegConvention::OBS_END    = "OBS_END";
const string CreditLegConvention::OBS_MIDDLE = "OBS_MIDDLE";

// strings describing the stub conventions
const string CreditLegConvention::FRONT   = "FRONT";
const string CreditLegConvention::BACK    = "BACK";
const string CreditLegConvention::LONG    = "LONG";
const string CreditLegConvention::SHORT   = "SHORT";
const string CreditLegConvention::TRIGGER = "TRIGGER";
const string CreditLegConvention::BOND    = "BOND";
const string CreditLegConvention::SIMPLE  = "SIMPLE";
const string CreditLegConvention::NONE    = "NONE";
        
/** TYPE (for reflection) */
CClassConstSP const CreditLegConvention::TYPE =
    CClass::registerClassLoadMethod(
        "CreditLegConvention",
        typeid(CreditLegConvention),
        load);

/****************/

bool CreditLegConventionLoad() { 
    return (CreditLegConvention::TYPE != 0); 
}

DRLIB_END_NAMESPACE
