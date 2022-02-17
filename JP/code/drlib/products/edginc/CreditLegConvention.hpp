//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CreditLegConvention.hpp
//
//   Description : Base class representing "credit legs" conventions (see below)
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef CREDIT_LEG_CONVENTION_BASE_HPP
#define CREDIT_LEG_CONVENTION_BASE_HPP

#include "edginc/ICreditLegConvention.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Base class representing "credit legs" conventions.
 * In particular, contains the "stub conventions" (in the following,
 * we define T0=start date, T1=first forthcoming coupon date, T2=coupon date
 * following T1, Te=end date):
 * 1. stubType ("front"/"back"):
 *    - "front" (or "first") means that you construct the coupon dates
 *      backward from the end date Te
 *    - "back" means that you construct the coupon dates onward from the
 *      start date T0 (rarely used)
 *
 * 2a. stubLength ("long"/"short"): determines what happens to the first
 *    (forthcoming) coupon T1
 *    - "long": the first coupon payment T1 is combined with the next one T2
 *    - "short": the first coupon payment T1 occurs normally
 *
 * 2b. stubLengthTrigger (eg: 1 month): determines which stub length to use
 *    - if (T1-T0)<stubLengthTrigger: uses "long"
 *    - otherwise: uses "short"
 *
 * 3. stubPayment ("bond"/"simple"/"none"):
 *    - "bond" means that you receive an upfront payment at T0 and then pay a full coupon at T1
 *    - "simple" means that there is no upfront payment and you pay at T1 the coupon
 *      accrued between T0 and T1
 *    - "none" is equivalent to "bond" with an upfront payment of 0.0
 * */
class PRODUCTS_DLL CreditLegConvention:
    public MarketObject,
    public virtual ICreditLegConvention {
public:
    /**
     * Returns the default start date given the trade date
     * [Implements ICreditLegConvention]
     * */
    virtual DateTime startDate(const DateTime& tradeDate) const;

    /**
     * Generates the fee leg given the end date
     * [Implements ICreditLegConvention]
     * */
    virtual ICreditFeeLegSP generateFeeLeg(
        const DateTime& startDate,
        const DateTime& endDate,
        double coupon,
        double upfrontPayment,
        double notional = 1.0) const;

    /**
     * Generates the contingent leg given the end date
     * [Implements ICreditLegConvention]
     * */
    virtual ICreditContingentLegSP generateContingentLeg(
        const DateTime& startDate,
        const DateTime& endDate,
        double coupon,
        double upfrontPayment,
        double notional = 1.0) const;

    /**
     * Access to discount curve name
     * [Implements ICreditLegConvention]
     * */
    virtual string getDiscountName() const;

    /**
     * Access to discount curve
     * [Implements ICreditLegConvention]
     * */
    virtual YieldCurveConstSP getDiscount() const;

    /** access to recoverNotional flag */
    virtual bool getRecoverNotional() const;

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Virtual destructor */
    virtual ~CreditLegConvention();

    /** Checks parameters immediately after object is constructed */
    virtual void validatePop2Object();

    /** Populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the name of this object */
    virtual string getName() const;

    // strings describing when we observe the coupon notional
    static const string OBS_START;
    static const string OBS_END;
    static const string OBS_MIDDLE;

    // strings describing the stub conventions
    static const string FRONT;
    static const string BACK;
    static const string LONG;
    static const string SHORT;
    static const string TRIGGER;
    static const string BOND;
    static const string SIMPLE;
    static const string NONE;

protected:
    /** Constructor (used by reflection) */
    CreditLegConvention(CClassConstSP clazz = TYPE);

    // --------------------------
    // Contingent leg conventions
    // --------------------------

    /** name of object */
    string name;

    /**
     * Number of calendar days you need to add to
     * the trade date to get the protection start date
     * [optional, default is 1]
     * */
    int startDateOffset;

    /**
     * Pay upon default (TRUE) or at payment date (FALSE)
     * [optional, default is TRUE]
     * */
    bool paymentUponDefault;

    /**
     * Number of delay days for default payment when paymentUponDefault = TRUE
     * [optional, default is (empirically) 7]
     * */
    int defaultPaymentDelayDays;

    // -------------------
    // Fee leg conventions
    // -------------------

    /** Payment day count convention [optional, default is "Act/360"] */
    string paymentDCC;

    /** Payment day count convention object [transient] */
    DayCountConventionSP paymentDCCSP;

    /**
     * Bad day convention used to adjust all payment dates, except the last one
     * (for the last payment date, we use "lastPaymentBDC"
     * [optional, default is "Modified"]
     * */
    string paymentBDC;

    /** Payment bad day convention object [transient] */
    BadDayConventionSP paymentBDCSP;

    /**
     * Bad day convention used to adjust the last payment date
     * [optional, default is "Following"]
     * */
    string lastPaymentBDC;

    /** Last payment bad day convention object [transient] */
    BadDayConventionSP lastPaymentBDCSP;

    /** Payment period [optional, default is 3M] */
    MaturityPeriodSP paymentPeriod;

    /** Holidays [mandatory] */
    HolidayWrapper holidays;

    /**
     * Describes when we observe the coupon notional over an
     * accrual fee period [T0, T1]:
     * - OBS_START  : at the beginning of the period (i.e. T0, quite exotic)
     * - OBS_END    : at the end of the period (i.e. T1)
     * - OBS_MIDDLE : in the middle of the period (i.e. (T0+T1)/2)
     * [optional, default is OBS_MIDDLE]
     * */
     //TODO: we should also have a "CONTINUOUS" observation type that
     //      corresponds to "payAccruedFee=true" in CDSPricer.
     //      Note that "PERIOD_END"  corresponds to "payAccruedFee=false" in CDSPricer.
     //      -> see CDS fee leg pricing in CDSPricer
    string notionalObservationType;

    /** is notional reduced by recovery amount (index convention) [Optional, default = false] */
    bool recoverNotional;

    // ----------------
    // Stub conventions
    // ----------------

    /** Stub type: "FRONT" or "BACK" [optional, default is "FRONT"] */
    string stubType;

    /**
     * Stub length ("LONG", "SHORT" or "TRIGGER")
     * [optional, default is "TRIGGER" i.e. use stubLengthTrigger]
     * */
    string stubLength;

    /** Stub length trigger [optional, default is "1M"] */
    MaturityPeriodSP stubLengthTrigger;

    /** Stub payment: "BOND", "SIMPLE" or "NONE" [optional, default is "SIMPLE"] */
    string stubPayment;

    // --------------
    // Discount curve
    // --------------

    /** Discount curve [mandatory] */
    YieldCurveWrapper discount;

private:
    /** Default constructor */
    static IObject* defaultConstructor();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif /*CREDIT_LEG_CONVENTION_BASE_HPP*/
