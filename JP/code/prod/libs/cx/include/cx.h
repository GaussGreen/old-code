/*
***************************************************************************
** HEADER FILE: cx.h
***************************************************************************
*/

#ifndef _CX_CX_H
#define _CX_CX_H


#include <cxutils/include/cxutils.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** Defines the type of credit curve.
    We have two bootstrapping methods for credit curves. The flow method
    is the standard for flow products (flat forwards, continuous integration
    of survival probability). The exotic method is designed for compatibility
    with the credit exotic pricing methods. */
typedef enum
{
/** Used for bootstrapping curves for use within more exotic products - none
    of these are in this library. */
    CX_CURVE_TYPE_EXOTIC,                    /* X */
/** Used for bootstrapping curves for use within credit flows products. */
    CX_CURVE_TYPE_FLOW                       /* F */
} CxTCreditCurveType;

typedef enum
{
/** Linear interpolation. */
    CX_RECOVERY_LINEAR_INTERP,               /* Linear */
/** Step - uses right hand value. */
    CX_RECOVERY_RIGHT_INTERP,                /* Right */
/** Step - uses left hand value. */
    CX_RECOVERY_LEFT_INTERP                  /* Left */
} CxTRecoveryInterp;

/** Protection payment convention for CDS in case of default. */
typedef enum
{
/** Protection payment is due at default (after a delay for arranging the
    details). */
    CX_PROT_PAY_DEF,                         /* Default */
/** Protection payment is only due at the maturity of the protection leg. */
    CX_PROT_PAY_MAT                          /* Maturity */
} CxTProtPayConv;

/** Accrual payment convention for CDS in case of default. */
typedef enum
{
/** No accrual in case of default. */
    CX_ACCRUAL_PAY_NONE,                     /* None */
/** Interest since the last accrual date is due in case of default. */
    CX_ACCRUAL_PAY_ALL                       /* All */
} CxTAccrualPayConv;

/** Credit curve - contains clean spreads and is used for computing
    survival probabilities.
   
    Note that the survival probability = 1 at the
    end of the base date of the curve, and that subsequently survival
    probability should be measured at the end of the day. */
typedef struct _CxTCreditCurve
{
    /** Defines how the curve was built, and is used to determine how the curve
        is used to compute survival probabilities. */
    CxTCreditCurveType type;
    /** We use an ALIB TCurve for representing the survival curve. This gives
        maximum ease of use with the credit exotic pricing modules which need
        us to use TCurve. */
    TCurve*         tc;
    /** Date interval used in integration (where relevant). */
    TDateInterval*  timestep;
} CxTCreditCurve;

typedef struct _CxTIrMarketData
{
    char*           ccy;
    TDate           valueDate;
    TCurve*         discCurve;
    TCurve*         liborCurve;
} CxTIrMarketData;

/** Defines term structure of recovery. For some names, there are sound
    business reasons for believing that the recovery rate for a default
    in the short term will be different than the recovery rate for a
    default in the long term. Typically this is due to the belief that
    a later default will follow a period of asset capitalization, which
    would reduce the recovery rate. */
typedef struct _CxTRecoveryCurve
{
    /** Defines array size for dates, recoveryRates. */
    int             numItems;
    /** Array of size numItems.
        Dates for which the recovery rate is defined changes. */
    TDate*          dates;
    /** Array of size numItems.
        Recovery rates at the corresponding dates. */
    double*         recoveryRates;
    /** Interpolation type for recovery. There are three types:-
        
        LEFT which indicates a step interpolation and the value to the left
        (i.e. on or before) the given date is used.
        
        RIGHT which indicates a step interpolation and the value to the right
        (i.e on or after) the given date is used.
        
        LINEAR which indicates linear interpolation between dates. */
    CxTRecoveryInterp interpType;
    /** Do we allow extrapolation before the first date? In this case, the value
        at the first date is used. */
    TBoolean        extrapBefore;
    /** Do we allow extrapolation after the last date? In this case, the value
        at the last date is used. */
    TBoolean        extrapAfter;
} CxTRecoveryCurve;

typedef struct _CxTCreditMarketData
{
    char*           id;
    TDate           defaultDate;
    CxTCreditCurve* cleanSpreadCurve;
    CxTRecoveryCurve* recoveryCurve;
} CxTCreditMarketData;

/** Contingent leg (a.k.a. protection leg). Defines notional amount and
    protection start and end dates. */
typedef struct _CxTContingentLeg
{
    /** Start date of protection. You are protected from the end of this date. */
    TDate           startDate;
    /** Defines array size for dates, notionals.
        Number of dates corresponding to potential changes in the
        notional amounts. For a vanilla CDS, you would have nbDates=1. */
    int             nbDates;
    /** Array of size nbDates.
        Dates where potentially the notional changes. You are protected
        for the corresponding notional from the end of the previous
        date in the list (or startDate if you are at the start of the list)
        to the end of this date. */
    TDate*          dates;
    /** Array of size nbDates.
        Notionals corresponding to the dates. */
    double*         notionals;
    CxTProtPayConv  payType;
    long            payDelay;
} CxTContingentLeg;

/** Fee leg with fixed payments. */
typedef struct _CxTFeeLeg
{
    /** Defines array size for accStartDates, accEndDates, payDates, notionals, couponRates. */
    int             nbDates;
    /** Array of size nbDates.
        Start date for calculating accrued interest. */
    TDate*          accStartDates;
    /** Array of size nbDates.
        End date for calculating accrued interest. */
    TDate*          accEndDates;
    /** Array of size nbDates.
        Payment date for each fee payment. */
    TDate*          payDates;
    /** Array of size nbDates.
        Notional for each fee payment. */
    double*         notionals;
    /** Array of size nbDates.
        Coupon rate for each fee payment. */
    double*         couponRates;
    /** Day count convention for computing fee payments and accruals in
        case of default. */
    CxTDayCountConv dcc;
    /** Determines how accruals are handled in case of default. */
    CxTAccrualPayConv accrualPayConv;
    /** Denotes whether observation of defaults is at the start of the day
        or the end of the day for the accrual start and end dates. */
    TBoolean        obsStartOfDay;
} CxTFeeLeg;

/** Describes the parameters for a CDS. */
typedef struct _CxTCdsType
{
    /** The interval between coupon payments. */
    TDateInterval*  couponInterval;
    /** The day count convention used for calculating coupon payments and accrual */
    CxTDayCountConv paymentDcc;
    /** The location and type of stub if the CDS start date and end date are not
        exactly on cycle. */
    CxTStubType     stubType;
    /** Whether there is an extra day of accrued interest. In all cases protection
        goes from the beginning of the startDate to the end of the endDate. In
        principle this means that there should be an extra day of accrued interest
        as well which is paid in the final coupon period. Should be TRUE. */
    TBoolean        extraDay;
    /** Should we pay accrued interest on default. Usually TRUE. */
    TBoolean        payAccOnDefault;
    /** Bad day convention used for adjusting intermediate payment and accrual
        dates, and for adjusting final payment date. */
    CxTBadDayConv   badDayConv;
} CxTCdsType;

typedef struct _CxTCds
{
    CxTCdsType*     type;
    TDate           startDate;
    TDate           endDate;
    double          notional;
} CxTCds;

/** Transformation of market CDS rates into an alternative set of dates.
   
    A typical situation is that we have true benchmark instruments but
    existing systems in the bank are only capable of using instruments with
    a particular set of dates (usually defined at calendar intervals from
    the start date).
   
    This structure is returned by the transformation function and contains
    a set of dates and CDS spreads which when used to bootstrap a curve
    are guaranteed to reprice the true benchmark instruments.
   
    Optionally in addition we return a tweak matrix which enables you to
    convert    positions computed in terms of the fictional benchmark
    instruments into  positions in terms of the true benchmarks. */
typedef struct _CxTCdsCurveDateAdj
{
    /** Defines array size for dates, spreads. */
    long            nbDate;
    /** Array of size nbDate.
        Maturity dates of the fictional benchmark instruments. */
    TDate*          dates;
    /** Array of size nbDate.
        CDS spreads corresponding to dates. */
    double*         spreads;
    /** Tweak matrix which enables you to convert between positions in terms of the
        fictional instruments into positions in terms of the true benchmark
        instruments. */
    TMatrix2D*      tweakMatrix;
} CxTCdsCurveDateAdj;

/** Outputs for calculating CDS carry */
typedef struct _CxTCdsCarryOutputs
{
    /** couponsPV is the PV of the physical coupons paid between
        today and the horizon date discounted to value date
        and supposing no default. */
    double          couponsPV;
    /** horizonAccrual is the PV of the accrual received
        on the following coupon date. */
    double          horizonAccrual;
    /** carry is always the sum of couponsPV and horizonAccrual. */
    double          carry;
} CxTCdsCarryOutputs;

/** Conventions for bootstrapping a CDS curve. We will use this mainly in the
    CdsSlideBreakeven function as a way to make the function call less heavy.
    For the moment, all the fields are compulsory. */
typedef struct _CxTCdsConventions
{
    /** Should accrued interest be paid on default. Usually set to TRUE. */
    TBoolean        payAccOnDefault;
    /** Interval between coupon payments. Can be undefined when 3M is
        assumed. */
    TDateInterval*  couponInterval;
    CxTDayCountConv paymentDCC;
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    CxTStubType     stubType;
    CxTCreditCurveType curveType;
    TDateInterval*  timestep;
    TDateInterval*  smoothInterval;
    /** Should protection include the start date. */
    TBoolean        protectStart;
    TBoolean        isPriceClean;
    /** Delay in payments after default (in days). */
    long            delay;
    /** Bad day convention for adjusting coupon payment dates. */
    CxTBadDayConv   badDayConv;
} CxTCdsConventions;

/**
***************************************************************************
** Converts CreditCurveType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxCreditCurveTypeToString (CxTCreditCurveType value);

/**
***************************************************************************
** Converts CreditCurveType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxCreditCurveTypeFromString (const char* str, CxTCreditCurveType *val);

/**
***************************************************************************
** Converts RecoveryInterp to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxRecoveryInterpToString (CxTRecoveryInterp value);

/**
***************************************************************************
** Converts RecoveryInterp from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxRecoveryInterpFromString (const char* str, CxTRecoveryInterp *val);

/**
***************************************************************************
** Converts ProtPayConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxProtPayConvToString (CxTProtPayConv value);

/**
***************************************************************************
** Converts ProtPayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxProtPayConvFromString (const char* str, CxTProtPayConv *val);

/**
***************************************************************************
** Converts AccrualPayConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxAccrualPayConvToString (CxTAccrualPayConv value);

/**
***************************************************************************
** Converts AccrualPayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxAccrualPayConvFromString (const char* str, CxTAccrualPayConv *val);

/**
***************************************************************************
** Constructor for CxTCreditCurve
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveMake(
/** Defines how the curve was built, and is used to determine how the curve
    is used to compute survival probabilities. */
CxTCreditCurveType type,
/** We use an ALIB TCurve for representing the survival curve. This gives
    maximum ease of use with the credit exotic pricing modules which need
    us to use TCurve. */
TCurve*         tc,
/** Date interval used in integration (where relevant). */
TDateInterval*  timestep
);

/**
***************************************************************************
** Memory allocator for CxTCreditCurve
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCreditCurve
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveCopy(CxTCreditCurve* src);

/**
***************************************************************************
** Destructor for CxTCreditCurve
***************************************************************************
*/
void CxCreditCurveFree(CxTCreditCurve *p);

/**
***************************************************************************
** Constructor for CxTIrMarketData
***************************************************************************
*/
CxTIrMarketData* CxIrMarketDataMake(
char*           ccy,
TDate           valueDate,
TCurve*         discCurve,
TCurve*         liborCurve
);

/**
***************************************************************************
** Memory allocator for CxTIrMarketData
***************************************************************************
*/
CxTIrMarketData* CxIrMarketDataMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTIrMarketData
***************************************************************************
*/
CxTIrMarketData* CxIrMarketDataCopy(CxTIrMarketData* src);

/**
***************************************************************************
** Destructor for CxTIrMarketData
***************************************************************************
*/
void CxIrMarketDataFree(CxTIrMarketData *p);

/**
***************************************************************************
** Constructor for CxTRecoveryCurve
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveMake(
/** Defines array size for dates, recoveryRates. */
int             numItems,
/** Array of size numItems.
    Dates for which the recovery rate is defined changes. */
TDate*          dates,
/** Array of size numItems.
    Recovery rates at the corresponding dates. */
double*         recoveryRates,
/** Interpolation type for recovery. There are three types:-
    
    LEFT which indicates a step interpolation and the value to the left
    (i.e. on or before) the given date is used.
    
    RIGHT which indicates a step interpolation and the value to the right
    (i.e on or after) the given date is used.
    
    LINEAR which indicates linear interpolation between dates. */
CxTRecoveryInterp interpType,
/** Do we allow extrapolation before the first date? In this case, the value
    at the first date is used. */
TBoolean        extrapBefore,
/** Do we allow extrapolation after the last date? In this case, the value
    at the last date is used. */
TBoolean        extrapAfter
);

/**
***************************************************************************
** Memory allocator for CxTRecoveryCurve
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveMakeEmpty(
/** Defines array size for dates, recoveryRates. */
int             numItems
);

/**
***************************************************************************
** Copy constructor for CxTRecoveryCurve
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveCopy(CxTRecoveryCurve* src);

/**
***************************************************************************
** Destructor for CxTRecoveryCurve
***************************************************************************
*/
void CxRecoveryCurveFree(CxTRecoveryCurve *p);

/**
***************************************************************************
** Constructor for CxTCreditMarketData
***************************************************************************
*/
CxTCreditMarketData* CxCreditMarketDataMake(
char*           id,
TDate           defaultDate,
CxTCreditCurve* cleanSpreadCurve,
CxTRecoveryCurve* recoveryCurve
);

/**
***************************************************************************
** Memory allocator for CxTCreditMarketData
***************************************************************************
*/
CxTCreditMarketData* CxCreditMarketDataMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCreditMarketData
***************************************************************************
*/
CxTCreditMarketData* CxCreditMarketDataCopy(CxTCreditMarketData* src);

/**
***************************************************************************
** Destructor for CxTCreditMarketData
***************************************************************************
*/
void CxCreditMarketDataFree(CxTCreditMarketData *p);

/**
***************************************************************************
** Constructor for CxTContingentLeg
***************************************************************************
*/
CxTContingentLeg* CxContingentLegMake(
/** Start date of protection. You are protected from the end of this date. */
TDate           startDate,
/** Defines array size for dates, notionals.
    Number of dates corresponding to potential changes in the
    notional amounts. For a vanilla CDS, you would have nbDates=1. */
int             nbDates,
/** Array of size nbDates.
    Dates where potentially the notional changes. You are protected
    for the corresponding notional from the end of the previous
    date in the list (or startDate if you are at the start of the list)
    to the end of this date. */
TDate*          dates,
/** Array of size nbDates.
    Notionals corresponding to the dates. */
double*         notionals,
CxTProtPayConv  payType,
long            payDelay
);

/**
***************************************************************************
** Memory allocator for CxTContingentLeg
***************************************************************************
*/
CxTContingentLeg* CxContingentLegMakeEmpty(
/** Defines array size for dates, notionals.
    Number of dates corresponding to potential changes in the
    notional amounts. For a vanilla CDS, you would have nbDates=1. */
int             nbDates
);

/**
***************************************************************************
** Copy constructor for CxTContingentLeg
***************************************************************************
*/
CxTContingentLeg* CxContingentLegCopy(CxTContingentLeg* src);

/**
***************************************************************************
** Destructor for CxTContingentLeg
***************************************************************************
*/
void CxContingentLegFree(CxTContingentLeg *p);

/**
***************************************************************************
** Constructor for CxTFeeLeg
***************************************************************************
*/
CxTFeeLeg* CxFeeLegMake(
/** Defines array size for accStartDates, accEndDates, payDates, notionals, couponRates. */
int             nbDates,
/** Array of size nbDates.
    Start date for calculating accrued interest. */
TDate*          accStartDates,
/** Array of size nbDates.
    End date for calculating accrued interest. */
TDate*          accEndDates,
/** Array of size nbDates.
    Payment date for each fee payment. */
TDate*          payDates,
/** Array of size nbDates.
    Notional for each fee payment. */
double*         notionals,
/** Array of size nbDates.
    Coupon rate for each fee payment. */
double*         couponRates,
/** Day count convention for computing fee payments and accruals in
    case of default. */
CxTDayCountConv dcc,
/** Determines how accruals are handled in case of default. */
CxTAccrualPayConv accrualPayConv,
/** Denotes whether observation of defaults is at the start of the day
    or the end of the day for the accrual start and end dates. */
TBoolean        obsStartOfDay
);

/**
***************************************************************************
** Memory allocator for CxTFeeLeg
***************************************************************************
*/
CxTFeeLeg* CxFeeLegMakeEmpty(
/** Defines array size for accStartDates, accEndDates, payDates, notionals, couponRates. */
int             nbDates
);

/**
***************************************************************************
** Copy constructor for CxTFeeLeg
***************************************************************************
*/
CxTFeeLeg* CxFeeLegCopy(CxTFeeLeg* src);

/**
***************************************************************************
** Destructor for CxTFeeLeg
***************************************************************************
*/
void CxFeeLegFree(CxTFeeLeg *p);

/**
***************************************************************************
** Constructor for CxTCdsType
***************************************************************************
*/
CxTCdsType* CxCdsTypeMake(
/** The interval between coupon payments. */
TDateInterval*  couponInterval,
/** The day count convention used for calculating coupon payments and accrual */
CxTDayCountConv paymentDcc,
/** The location and type of stub if the CDS start date and end date are not
    exactly on cycle. */
CxTStubType     stubType,
/** Whether there is an extra day of accrued interest. In all cases protection
    goes from the beginning of the startDate to the end of the endDate. In
    principle this means that there should be an extra day of accrued interest
    as well which is paid in the final coupon period. Should be TRUE. */
TBoolean        extraDay,
/** Should we pay accrued interest on default. Usually TRUE. */
TBoolean        payAccOnDefault,
/** Bad day convention used for adjusting intermediate payment and accrual
    dates, and for adjusting final payment date. */
CxTBadDayConv   badDayConv
);

/**
***************************************************************************
** Memory allocator for CxTCdsType
***************************************************************************
*/
CxTCdsType* CxCdsTypeMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCdsType
***************************************************************************
*/
CxTCdsType* CxCdsTypeCopy(CxTCdsType* src);

/**
***************************************************************************
** Destructor for CxTCdsType
***************************************************************************
*/
void CxCdsTypeFree(CxTCdsType *p);

/**
***************************************************************************
** Constructor for CxTCds
***************************************************************************
*/
CxTCds* CxCdsMake(
CxTCdsType*     type,
TDate           startDate,
TDate           endDate,
double          notional
);

/**
***************************************************************************
** Memory allocator for CxTCds
***************************************************************************
*/
CxTCds* CxCdsMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCds
***************************************************************************
*/
CxTCds* CxCdsCopy(CxTCds* src);

/**
***************************************************************************
** Destructor for CxTCds
***************************************************************************
*/
void CxCdsFree(CxTCds *p);

/**
***************************************************************************
** Constructor for CxTCdsCurveDateAdj
***************************************************************************
*/
CxTCdsCurveDateAdj* CxCdsCurveDateAdjMake(
/** Defines array size for dates, spreads. */
long            nbDate,
/** Array of size nbDate.
    Maturity dates of the fictional benchmark instruments. */
TDate*          dates,
/** Array of size nbDate.
    CDS spreads corresponding to dates. */
double*         spreads,
/** Tweak matrix which enables you to convert between positions in terms of the
    fictional instruments into positions in terms of the true benchmark
    instruments. */
TMatrix2D*      tweakMatrix
);

/**
***************************************************************************
** Memory allocator for CxTCdsCurveDateAdj
***************************************************************************
*/
CxTCdsCurveDateAdj* CxCdsCurveDateAdjMakeEmpty(
/** Defines array size for dates, spreads. */
long            nbDate
);

/**
***************************************************************************
** Copy constructor for CxTCdsCurveDateAdj
***************************************************************************
*/
CxTCdsCurveDateAdj* CxCdsCurveDateAdjCopy(CxTCdsCurveDateAdj* src);

/**
***************************************************************************
** Destructor for CxTCdsCurveDateAdj
***************************************************************************
*/
void CxCdsCurveDateAdjFree(CxTCdsCurveDateAdj *p);

/**
***************************************************************************
** Constructor for CxTCdsCarryOutputs
***************************************************************************
*/
CxTCdsCarryOutputs* CxCdsCarryOutputsMake(
/** couponsPV is the PV of the physical coupons paid between
    today and the horizon date discounted to value date
    and supposing no default. */
double          couponsPV,
/** horizonAccrual is the PV of the accrual received
    on the following coupon date. */
double          horizonAccrual,
/** carry is always the sum of couponsPV and horizonAccrual. */
double          carry
);

/**
***************************************************************************
** Memory allocator for CxTCdsCarryOutputs
***************************************************************************
*/
CxTCdsCarryOutputs* CxCdsCarryOutputsMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCdsCarryOutputs
***************************************************************************
*/
CxTCdsCarryOutputs* CxCdsCarryOutputsCopy(CxTCdsCarryOutputs* src);

/**
***************************************************************************
** Destructor for CxTCdsCarryOutputs
***************************************************************************
*/
void CxCdsCarryOutputsFree(CxTCdsCarryOutputs *p);

/**
***************************************************************************
** Constructor for CxTCdsConventions
***************************************************************************
*/
CxTCdsConventions* CxCdsConventionsMake(
/** Should accrued interest be paid on default. Usually set to TRUE. */
TBoolean        payAccOnDefault,
/** Interval between coupon payments. Can be undefined when 3M is
    assumed. */
TDateInterval*  couponInterval,
CxTDayCountConv paymentDCC,
/** If the startDate and endDate are not on cycle, then this parameter
    determines location of coupon dates. */
CxTStubType     stubType,
CxTCreditCurveType curveType,
TDateInterval*  timestep,
TDateInterval*  smoothInterval,
/** Should protection include the start date. */
TBoolean        protectStart,
TBoolean        isPriceClean,
/** Delay in payments after default (in days). */
long            delay,
/** Bad day convention for adjusting coupon payment dates. */
CxTBadDayConv   badDayConv
);

/**
***************************************************************************
** Memory allocator for CxTCdsConventions
***************************************************************************
*/
CxTCdsConventions* CxCdsConventionsMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCdsConventions
***************************************************************************
*/
CxTCdsConventions* CxCdsConventionsCopy(CxTCdsConventions* src);

/**
***************************************************************************
** Destructor for CxTCdsConventions
***************************************************************************
*/
void CxCdsConventionsFree(CxTCdsConventions *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _CX_CX_H */
