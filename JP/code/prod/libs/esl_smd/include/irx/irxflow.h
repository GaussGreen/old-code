/*
***************************************************************************
** HEADER FILE: irxflow.h
**
** Defines data structures and enums used in the irxflow library.
***************************************************************************
*/

#ifndef _IRX_IRXFLOW_H
#define _IRX_IRXFLOW_H


#include <irx/irxutils.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** Defines the bootstrap method for building a zero curve from swap rates. */
typedef enum
{
    IRX_BOOTSTRAP_FLAT,                      /* Flat */
    IRX_BOOTSTRAP_SMOOTH                     /* Smooth */
} IrxTBootstrapMethod;

/** Defines the type of rates. In the formulae, r = rate and t = time. */
typedef enum
{
/** Simple interest rate - the discount factor is
    \begin{equation*}
    \frac{1}{1+rt}
    \end{equation*} */
    IRX_RT_SIMPLE,                           /* Simple */
/** Discount rate - the discount factor is
    \begin{equation*}
    1-rt
    \end{equation*} */
    IRX_RT_DISCOUNT,                         /* Discount */
/** Continuously compounded - the discount factor is
    \begin{equation*}
    e^{-rt}
    \end{equation*} */
    IRX_RT_CONTINUOUS,                       /* Continuous */
/** Compounded N periods per year - the discount factor is
    \begin{equation*}
    (1+\frac{r}{N})^{-Nt}
    \end{equation*} */
    IRX_RT_FREQUENCY                         /* Frequency */
} IrxTRateTypeType;

typedef enum
{
/** Generates one cash flow - partial coupon at first coupon date. */
    IRX_STUB_SIMPLE = 1,                     /* Simple */
/** Generates two cash flows - full coupon at first coupon date and accrued
    interest (of opposite sign) at value date. */
    IRX_STUB_BOND = 2,                       /* Bond */
/** Generates one cash flow - full coupon at first coupon date. */
    IRX_STUB_NONE = 4                        /* None */
} IrxTStubPayment;

typedef struct _IrxTCashFlowList
{
    int             numItems;
    /** Array of size numItems. */
    IrxTDate*       dates;
    /** Array of size numItems. */
    double*         amounts;
} IrxTCashFlowList;

/** Defines a zero curve structure which uses either flat forwards or
    parabolic forwards. */
typedef struct _IrxTZeroCurve
{
    /** Base date of the curve. The discount factor at this date should be 1. */
    IrxTDate        baseDate;
    /** The maximum date that can be used for extrapolating beyond the last date
        in the curve. This can be 0 in which case extrapolation is not allowed. */
    IrxTDate        maxDate;
    /** Number of points in the curve */
    int             numItems;
    /** Array of size numItems.
        Start dates for the discount factor. One of the dates in the curve must be
        the baseDate */
    IrxTDate*       startDates;
    /** Array of size numItems.
        Zero coupon prices corresponding to startDates.
        The zero coupon price corresponding to baseDate must be 1.
        a.k.a. discountFactors */
    double*         prices;
    /** Array of size numItems.
        If using flat forwards, then a0 is the forward rate
        from the corresponding startDate to the next date in the curve.
        Otherwise it is the constaint component in the parabola. */
    double*         a0;
    /** Array of size numItems.
        Linear component in the parabola.
        For flat forwards, a1=0. */
    double*         a1;
    /** Array of size numItems.
        Quadratic component in the parabola.
        For flat forwards, a2=0. */
    double*         a2;



    /************************************************************************
     *  
     * The followings are for backward compatibility with T_CURVE, mostly
     * requirments by wrapper I/O routines.  These should only stay here
     * for the transition purpose only, and in principle we should use 
     * IrxTSwapZeroCurve structure throughout.
     *
     ************************************************************************/

    /* Base date information */
    IrxTDate    Today;               /**< Today's date                        */
    int         SpotDays;            /**< Spot days                           */
    IrxTDate    ValueDate;           /**< Value date, = baseDate              */

    /* Underlying yield curve conventions */
    char        SwapFreq;            /**< Benchmark swap frequency            */
    char        SwapDCC[4];          /**< Benchmark swap day count convention */
    int         MMB;                 /**< Money market basis (360 or 365)     */

} IrxTZeroCurve;

/** Defines the interest rate market parameters. You will be able to coerce
    these parameters from currency codes for which the values are hard-coded
    within the library. */
typedef struct _IrxTMarketConv
{
    /** Date interval for the fixed leg */
    IrxTDateInterval fixedIvl;
    /** Date interval for the floating leg */
    IrxTDateInterval floatIvl;
    /** Date interval for currency basis swaps (into USD) */
    IrxTDateInterval cbsIvl;
    /** Day count convention for the fixed leg */
    IrxTDayCountConv fixedDcc;
    /** Day count convention for the floating leg */
    IrxTDayCountConv floatDcc;
    /** Day count convention for currency basis swaps */
    IrxTDayCountConv cbsDcc;
    /** Day count convention for the money market payment */
    IrxTDayCountConv mmDcc;
    /** Bad day convention for payments */
    IrxTBadDayConv  paymentBdc;
    /** Bad day convention for accruals */
    IrxTBadDayConv  accrualBdc;
    /** Bad day convention for reset of floating leg */
    IrxTBadDayConv  resetBdc;
    /** Bad day convention for money market payments.
        
        Note that if this results in a payment adjustment to on or before the
        start date, then it will automatically switch to use following.
        (Prime example - the 1D rate when the start date is the last business day
        of the month, but not the last calendar day of the month). */
    IrxTBadDayConv  mmBdc;
    /** Number of days to spot */
    int             daysToSpot;
} IrxTMarketConv;

/** Defines the basis of a rate - where the basis is like simple basis,
    continuously compouneded basis etc. */
typedef struct _IrxTRateType
{
    /** Type of the rate - e.g. simple, continously compounded, compounded N periods
        per year etc. */
    IrxTRateTypeType type;
    /** Frequency - only used for rates with N periods per year. This is N. */
    int             frequency;
} IrxTRateType;

typedef struct _IrxTSwapZeroCurve
{
    /** Today - normally before baseDate of the curves. */
    IrxTDate        today;
    /** Market conventions */
    IrxTMarketConv* marketConv;
    /** Holiday */
    IrxTCalendar*   calendar;
    /** Bootstrapped zero curve for discounting */
    IrxTZeroCurve*  discountCurve;
    /** Bootstrapped zero curve for estimation.
        If undefined, then we use the discount curve. */
    IrxTZeroCurve*  indexCurve;
} IrxTSwapZeroCurve;

/** Defines a swap object - which incorporates a fixed leg and a floating
    leg, i.e. a vanilla interest rate swap.
   
    Regarding the first coupon of the floating leg - there is the question
    on how we handle the first coupon when it is a broken dated period.
    Currently we assume that we use the broken period instead of the full
    period. */
typedef struct _IrxTSwap
{
    /** Market conventions - defines things such as DCC, BDC, date intervals */
    IrxTMarketConv* marketConv;
    /** Coupon rate for fixed payments */
    double          couponRate;
    /** Spread above index curve for floating payments */
    double          spread;
    /** Start date of the swap */
    IrxTDate        startDate;
    /** Maturity date of the swap */
    IrxTDate        maturityDate;
    /** Roll date for cash flow dates - typically the startDate */
    IrxTDate        rollDate;
    /** Location of the stub - typically at the front, but might be short or
        long. */
    IrxTStubLocation stubLocation;
    /** Payment for the first coupon. Whether we get a full coupon (or not),
        and if we get a full coupon whether we also get accrued interest. */
    IrxTStubPayment fixedStubPayment;
    /** Payment for the first coupon. Whether we get a full coupon (or not),
        and if we get a full coupon whether we also get accrued interest.
        Note that for a floating payment a BOND stub is only possible when
        there is no requirement to estimate the first coupon payment. */
    IrxTStubPayment floatStubPayment;
} IrxTSwap;

typedef struct _IrxTCouponPayment
{
    /** Payment is multiplied by the notional */
    double          notional;
    /** The coupon rate. Note that when we use a floating payment, this number
        is meaningless. */
    double          couponRate;
    /** The spread. This is added to the coupon rate to compute the payment. */
    double          spread;
    /** Index curve number. Value of 0 indicates that this is a fixed payment.
        Otherwise the index curve number is the curve that we use to compute
        the coupon. */
    int             indexCurveNumber;
    /** Rate is calculated as SIMPLE using this day count convention */
    IrxTDayCountConv dcc;
    /** Interest accrues from accrueStartDate */
    IrxTDate        accrueStartDate;
    /** Interest accrues to accrueEndDate */
    IrxTDate        accrueEndDate;
    /** Rate is observed effective at rateStartDate */
    IrxTDate        rateStartDate;
    /** Rate matures at rateEndDate */
    IrxTDate        rateEndDate;
} IrxTCouponPayment;

typedef struct _IrxTSwapPayments
{
    int             numFlows;
    /** Array of size numFlows.
        Payment dates */
    IrxTDate*       dates;
    /** Array of size numFlows.
        Amounts includes notional payments only. */
    double*         amounts;
    /** Array of size numFlows.
        Includes fixed and floating coupon payments and accrued interest. */
    IrxTCouponPayment* couponPayments;
} IrxTSwapPayments;

/**
***************************************************************************
** Converts BootstrapMethod to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxBootstrapMethodToString (IrxTBootstrapMethod value);

/**
***************************************************************************
** Converts BootstrapMethod from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxBootstrapMethodFromString (const char* str, IrxTBootstrapMethod *val);

/**
***************************************************************************
** Converts RateTypeType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxRateTypeTypeToString (IrxTRateTypeType value);

/**
***************************************************************************
** Converts RateTypeType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxRateTypeTypeFromString (const char* str, IrxTRateTypeType *val);

/**
***************************************************************************
** Converts StubPayment to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxStubPaymentToString (IrxTStubPayment value);

/**
***************************************************************************
** Converts StubPayment from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxStubPaymentFromString (const char* str, IrxTStubPayment *val);

/**
***************************************************************************
** Constructor for IrxTCashFlowList
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMake(
int             numItems,
/** Array of size numItems. */
IrxTDate const* dates,
/** Array of size numItems. */
double const*   amounts
);

/**
***************************************************************************
** Memory allocator for IrxTCashFlowList
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMakeEmpty(
int             numItems
);

/**
***************************************************************************
** Copy constructor for IrxTCashFlowList
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListCopy(IrxTCashFlowList const* src);

/**
***************************************************************************
** Destructor for IrxTCashFlowList
***************************************************************************
*/
void irxCashFlowListFree(IrxTCashFlowList *p);

/**
***************************************************************************
** Constructor for IrxTZeroCurve
***************************************************************************
*/
IrxTZeroCurve* irxZeroCurveMake(
/** Base date of the curve. The discount factor at this date should be 1. */
IrxTDate        baseDate,
/** The maximum date that can be used for extrapolating beyond the last date
    in the curve. This can be 0 in which case extrapolation is not allowed. */
IrxTDate        maxDate,
/** Number of points in the curve */
int             numItems,
/** Array of size numItems.
    Start dates for the discount factor. One of the dates in the curve must be
    the baseDate */
IrxTDate const* startDates,
/** Array of size numItems.
    Zero coupon prices corresponding to startDates.
    The zero coupon price corresponding to baseDate must be 1.
    a.k.a. discountFactors */
double const*   prices,
/** Array of size numItems.
    If using flat forwards, then a0 is the forward rate
    from the corresponding startDate to the next date in the curve.
    Otherwise it is the constaint component in the parabola. */
double const*   a0,
/** Array of size numItems.
    Linear component in the parabola.
    For flat forwards, a1=0. */
double const*   a1,
/** Array of size numItems.
    Quadratic component in the parabola.
    For flat forwards, a2=0. */
double const*   a2
);



/**
***************************************************************************
** Set legacy default parameters
***************************************************************************
*/
int IrxTZeroCurveSetDefault(IrxTZeroCurve *zc);


/**
***************************************************************************
** Set legacy values for IrxTZeroCurve from market
***************************************************************************
*/
int IrxTZeroCurveSetFromMarket(IrxTZeroCurve            *zc,
                               IrxTDate                 today,
                               const IrxTMarketConv     *marketConv);

/**
***************************************************************************
** Memory allocator for IrxTZeroCurve
***************************************************************************
*/
IrxTZeroCurve* irxZeroCurveMakeEmpty(
/** Number of points in the curve */
int             numItems
);

/**
***************************************************************************
** Memory allocator for IrxTZeroCurve
***************************************************************************
*/
int irxZeroCurveConstructEmpty(
IrxTZeroCurve* crv,
int            numItems
);

/**
***************************************************************************
** Copy constructor for IrxTZeroCurve
***************************************************************************
*/
IrxTZeroCurve* irxZeroCurveCopy(IrxTZeroCurve const* src);

/**
***************************************************************************
** Destructor + Delete for IrxTZeroCurve
***************************************************************************
*/
void irxZeroCurveFree(IrxTZeroCurve *p);

/**
***************************************************************************
** Destructor for IrxTZeroCurve
***************************************************************************
*/
void irxZeroCurveDestroy(IrxTZeroCurve *p);

/**
***************************************************************************
** Constructor for IrxTMarketConv
***************************************************************************
*/
IrxTMarketConv* irxMarketConvMake(
/** Date interval for the fixed leg */
IrxTDateInterval fixedIvl,
/** Date interval for the floating leg */
IrxTDateInterval floatIvl,
/** Date interval for currency basis swaps (into USD) */
IrxTDateInterval cbsIvl,
/** Day count convention for the fixed leg */
IrxTDayCountConv fixedDcc,
/** Day count convention for the floating leg */
IrxTDayCountConv floatDcc,
/** Day count convention for currency basis swaps */
IrxTDayCountConv cbsDcc,
/** Day count convention for the money market payment */
IrxTDayCountConv mmDcc,
/** Bad day convention for payments */
IrxTBadDayConv  paymentBdc,
/** Bad day convention for accruals */
IrxTBadDayConv  accrualBdc,
/** Bad day convention for reset of floating leg */
IrxTBadDayConv  resetBdc,
/** Bad day convention for money market payments.
    
    Note that if this results in a payment adjustment to on or before the
    start date, then it will automatically switch to use following.
    (Prime example - the 1D rate when the start date is the last business day
    of the month, but not the last calendar day of the month). */
IrxTBadDayConv  mmBdc,
/** Number of days to spot */
int             daysToSpot
);

/**
***************************************************************************
** Memory allocator for IrxTMarketConv
***************************************************************************
*/
IrxTMarketConv* irxMarketConvMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTMarketConv
***************************************************************************
*/
IrxTMarketConv* irxMarketConvCopy(IrxTMarketConv const* src);

/**
***************************************************************************
** Destructor for IrxTMarketConv
***************************************************************************
*/
void irxMarketConvFree(IrxTMarketConv *p);

/**
***************************************************************************
** Constructor for IrxTRateType
***************************************************************************
*/
IrxTRateType* irxRateTypeMake(
/** Type of the rate - e.g. simple, continously compounded, compounded N periods
    per year etc. */
IrxTRateTypeType type,
/** Frequency - only used for rates with N periods per year. This is N. */
int             frequency
);

/**
***************************************************************************
** Memory allocator for IrxTRateType
***************************************************************************
*/
IrxTRateType* irxRateTypeMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTRateType
***************************************************************************
*/
IrxTRateType* irxRateTypeCopy(IrxTRateType const* src);

/**
***************************************************************************
** Destructor for IrxTRateType
***************************************************************************
*/
void irxRateTypeFree(IrxTRateType *p);

/**
***************************************************************************
** Constructor for IrxTSwapZeroCurve
***************************************************************************
*/
IrxTSwapZeroCurve* irxSwapZeroCurveMake(
/** Today - normally before baseDate of the curves. */
IrxTDate        today,
/** Market conventions */
IrxTMarketConv const* marketConv,
/** Holiday */
IrxTCalendar const* calendar,
/** Bootstrapped zero curve for discounting */
IrxTZeroCurve const* discountCurve,
/** Bootstrapped zero curve for estimation.
    If undefined, then we use the discount curve. */
IrxTZeroCurve const* indexCurve
);

/**
***************************************************************************
** Memory allocator for IrxTSwapZeroCurve
***************************************************************************
*/
IrxTSwapZeroCurve* irxSwapZeroCurveMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTSwapZeroCurve
***************************************************************************
*/
IrxTSwapZeroCurve* irxSwapZeroCurveCopy(IrxTSwapZeroCurve const* src);

/**
***************************************************************************
** Destructor for IrxTSwapZeroCurve
***************************************************************************
*/
void irxSwapZeroCurveFree(IrxTSwapZeroCurve *p);

/**
***************************************************************************
** Constructor for IrxTSwap
***************************************************************************
*/
IrxTSwap* irxSwapMake(
/** Market conventions - defines things such as DCC, BDC, date intervals */
IrxTMarketConv const* marketConv,
/** Coupon rate for fixed payments */
double          couponRate,
/** Spread above index curve for floating payments */
double          spread,
/** Start date of the swap */
IrxTDate        startDate,
/** Maturity date of the swap */
IrxTDate        maturityDate,
/** Roll date for cash flow dates - typically the startDate */
IrxTDate        rollDate,
/** Location of the stub - typically at the front, but might be short or
    long. */
IrxTStubLocation stubLocation,
/** Payment for the first coupon. Whether we get a full coupon (or not),
    and if we get a full coupon whether we also get accrued interest. */
IrxTStubPayment fixedStubPayment,
/** Payment for the first coupon. Whether we get a full coupon (or not),
    and if we get a full coupon whether we also get accrued interest.
    Note that for a floating payment a BOND stub is only possible when
    there is no requirement to estimate the first coupon payment. */
IrxTStubPayment floatStubPayment
);

/**
***************************************************************************
** Memory allocator for IrxTSwap
***************************************************************************
*/
IrxTSwap* irxSwapMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTSwap
***************************************************************************
*/
IrxTSwap* irxSwapCopy(IrxTSwap const* src);

/**
***************************************************************************
** Destructor for IrxTSwap
***************************************************************************
*/
void irxSwapFree(IrxTSwap *p);

/**
***************************************************************************
** Constructor for IrxTCouponPayment
***************************************************************************
*/
IrxTCouponPayment* irxCouponPaymentMake(
/** Payment is multiplied by the notional */
double          notional,
/** The coupon rate. Note that when we use a floating payment, this number
    is meaningless. */
double          couponRate,
/** The spread. This is added to the coupon rate to compute the payment. */
double          spread,
/** Index curve number. Value of 0 indicates that this is a fixed payment.
    Otherwise the index curve number is the curve that we use to compute
    the coupon. */
int             indexCurveNumber,
/** Rate is calculated as SIMPLE using this day count convention */
IrxTDayCountConv dcc,
/** Interest accrues from accrueStartDate */
IrxTDate        accrueStartDate,
/** Interest accrues to accrueEndDate */
IrxTDate        accrueEndDate,
/** Rate is observed effective at rateStartDate */
IrxTDate        rateStartDate,
/** Rate matures at rateEndDate */
IrxTDate        rateEndDate
);

/**
***************************************************************************
** Memory allocator for IrxTCouponPayment
***************************************************************************
*/
IrxTCouponPayment* irxCouponPaymentMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTCouponPayment
***************************************************************************
*/
IrxTCouponPayment* irxCouponPaymentCopy(IrxTCouponPayment const* src);

/**
***************************************************************************
** Destructor for IrxTCouponPayment
***************************************************************************
*/
void irxCouponPaymentFree(IrxTCouponPayment *p);

/**
***************************************************************************
** Constructor for IrxTSwapPayments
***************************************************************************
*/
IrxTSwapPayments* irxSwapPaymentsMake(
int             numFlows,
/** Array of size numFlows.
    Payment dates */
IrxTDate const* dates,
/** Array of size numFlows.
    Amounts includes notional payments only. */
double const*   amounts,
/** Array of size numFlows.
    Includes fixed and floating coupon payments and accrued interest. */
IrxTCouponPayment const* couponPayments
);

/**
***************************************************************************
** Memory allocator for IrxTSwapPayments
***************************************************************************
*/
IrxTSwapPayments* irxSwapPaymentsMakeEmpty(
int             numFlows
);

/**
***************************************************************************
** Copy constructor for IrxTSwapPayments
***************************************************************************
*/
IrxTSwapPayments* irxSwapPaymentsCopy(IrxTSwapPayments const* src);

/**
***************************************************************************
** Destructor for IrxTSwapPayments
***************************************************************************
*/
void irxSwapPaymentsFree(IrxTSwapPayments *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _IRX_IRXFLOW_H */
