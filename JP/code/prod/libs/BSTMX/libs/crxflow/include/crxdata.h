/*
***************************************************************************
** HEADER FILE: crxdata.h
**
** Defines data structures for crxflow library.
***************************************************************************
*/

#ifndef _CRX_CRXDATA_H
#define _CRX_CRXDATA_H


#include <cxutils/include/cxutils.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** Determines how a CDS option pays off at expiry. */
typedef enum
{
/** Uses the market spread at expiry to compute the settlement value of the CDS. */
    CRX_CDS_OPTION_PAYOFF_MARKET,            /* Market */
/** Uses the strike spread at expiry to compute the settlement value of the CDS.
    This has some surprising implications for the calculation of the forward
    spread used within the model. */
    CRX_CDS_OPTION_PAYOFF_STRIKE             /* Strike */
} CrxTCdsOptionPayoff;

/** Whether a CDS option is struck on price or spread. The default is spread
    (for historical purposes). */
typedef enum
{
/** The option is struck on spread - quoted as a decimal, e.g. 0.03 is 300
    basis points. */
    CRX_CDS_OPTION_STRIKE_TYPE_SPREAD,       /* Spread */
/** The option is struck n price - quoted as a decimal with a par of 1. */
    CRX_CDS_OPTION_STRIKE_TYPE_PRICE         /* Price */
} CrxTCdsOptionStrikeType;

/** The distribution of the forward is assumed to be lognormal at option
    expiry. The Q-distribution is ignored. */
#define CRX_DIST_TYPE_LOGNORMAL (long)'L'              /* Lognormal */
/** The distribution of the forward is assumed to be normal at option expiry.
    The Q-distribution is ignored. */
#define CRX_DIST_TYPE_NORMAL (long)'N'              /* Normal */
/** The distribution of the forward is assumed to follow the Q-distribution
    which contains upto 6 values of q which apply at different points in the
    distribution. */
#define CRX_DIST_TYPE_Q      (long)'Q'              /* Q */

/** Call option. Note that for a CDS option a call is the right to buy risk,
    which translates to the right to sell protection at a particular spread. */
#define CRX_OPTION_TYPE_CALL (long)'C'              /* Call */
/** Put option. Note that for a CDS option a put is the right to sell risk,
    which translates to the right to buy protection at a particular spread. */
#define CRX_OPTION_TYPE_PUT  (long)'P'              /* Put */
/** Straddle - a call and a put at the same strike. Mostly useful for
    calculating the implied volatility of a straddle, since of course you
    can just add the prices of the put and the call. */
#define CRX_OPTION_TYPE_STRADDLE (long)'S'              /* Straddle */

#define CRX_ACT_ACT          (long)1                /* ACT/ACT */
#define CRX_ACT_365F         (long)2                /* ACT/365F */
#define CRX_ACT_360          (long)3                /* ACT/360 */
#define CRX_B30_360          (long)4                /* 30/360 */
#define CRX_B30E_360         (long)5                /* 30E/360 */

/** Defines multi-Q distribution. */
typedef struct _CrxTQDist
{
    /** Defines array size for Qs.
        Number of Qs. */
    int             nQs;
    /** Array of size nQs.
        Q-factors. A value of Q=0 indicates no skew - i.e. lognormal.
        A value of Q=1 indicates a normal distribution.
        Typical range of values would be [-1,2] but we have seen values outside
        this range. There should be precisely 6 Q-values.
        
        If you want to use just 2Qs, then you are at liberty to set the values
        on the left equal to each other, and the values on the right likewise. */
    double*         Qs;
    /** Defines array size for Ds.
        Number of Ds. It should be that nDs = nQs-1 */
    int             nDs;
    /** Array of size nDs.
        Delta-values. These indicate the transition points between the Q-values.
        Typically the middle D will be 0.5. Values should be increasing and in
        the range 0-1 exclusive. There should be one less D-value that Q-value,
        i.e. precisely 5 values. */
    double*         Ds;
} CrxTQDist;

/** Instrument details for a CDS option. This describes a European option on a
    CDS (index or single name). A call option is the right to buy risk (sell
    protection). A put option is the right to sell risk (buy protection). */
typedef struct _CrxTCdsOption
{
    /** Issue date of the underlying CDS. */
    TDate           issueDate;
    /** Maturity date of the underlying CDS. */
    TDate           maturityDate;
    /** Fee interval of the underlying CDS. Usual value is 3M. */
    TDateInterval*  feeInterval;
    /** Day count convention of the underlying CDS. Usual value is `ACT/360'. */
    long            dcc;
    /** Does the underlying CDS pay accrued interest on default. Usually TRUE. */
    TBoolean        payAccOnDefault;
    /** Does the underlying CDS represent an index instead of a single name. */
    TBoolean        isIndex;
    /** Loss so far - only relevant for index CDS. */
    double          lossSoFar;
    /** Option type - i.e. call or put. */
    long            optionType;
    /** Exercise date. */
    TDate           exerciseDate;
    /** Payment date - when the exercise payment is made. Must be on or after the
        exercise date. */
    TDate           paymentDate;
    /** Flag to denote that whether the option knock-out on default.
        Only relevant for single name options. */
    TBoolean        koOnDefault;
    /** Strike spread of the option. Measured in decimals, e.g. 100bp
        is represented as 0.01 */
    double          strike;
    /** Option payoff method - only relevant for index options. */
    CrxTCdsOptionPayoff optionPayoff;
    /** Coupon of the underlying CDS - only relevant for index options
        which pay-off using the strike method. */
    double          coupon;
    /** Strike type - i.e. price or spread */
    CrxTCdsOptionStrikeType strikeType;
} CrxTCdsOption;

/** Outputs from the CDS option calculators.
   
    Note that depending on whether the CDS option is price struck or spread
    struck you will get different values populated here. For spread strike,
    the forward spread and annuity is returned. For price strike, the forward
    price is returned. The vol is either a spread vol or price vol depending
    upon the type of strike. */
typedef struct _CrxTCdsOptionCalc
{
    /** Price of the option in money terms. */
    double          price;
    /** At the money volatility - either implied or as input. */
    double          vol;
    /** Forward spread used in the option formula. */
    double          fwdSpread;
    /** Forward spread without taking into account index forward adjustments. */
    double          fwdSpreadUnadj;
    /** Risky annuity measured at the value date of the option for the forward. */
    double          annuity;
    /** Time to expiry (exerciseDate-today in ACT/365F) */
    double          timeToExpiry;
    /** Forward price at expiry */
    double          fwdPrice;
} CrxTCdsOptionCalc;

/** State of the Q-distribution optimization.
   
    The optimizer is designed to return the current state of the optimization
    between function calls. Then you can provide the output of one call to
    the optimizer as an input to the next call of the optimizer. This mechanism
    is provided to enable you to get feedback between steps of the optimization,
    since the optimization is potentially quite a slow procedure.
   
    Before calling the optimizer for the first time, then you need to
    initialize the state. This is done by providing just the first two fields
    of the structure (qdist, vol). */
typedef struct _CrxTCdsOptionQOptimization
{
    /** Q-distribution */
    CrxTQDist*      qdist;
    /** At the money volatility */
    double          vol;
    /** Initial directions for search (on initialization). Final directions for
        search (on exit). */
    TMatrix2D*      direction;
    /** Number of iterations so far. */
    int             iter;
    /** Root mean square difference for the optimization. The difference is
        computed for strike volatilities. */
    double          value;
    /** Most recent absolute improvement in the algorithm. */
    double          vdiff;
} CrxTCdsOptionQOptimization;

/** Q-optimization model parameters.
   
    These are used to control the optimizer routine which finds the best fit
    for the multi-Q parameters given market data. */
typedef struct _CrxTCdsOptionQOptModel
{
    /** Denotes whether we optimize the ATM volatility. Sometimes a better fit
        can be obtained by optimizing the ATM volatility, but in general we
        would expect the ATM volatility to be defined and not changed. */
    TBoolean        optimizeVol;
    /** Defines array size for qsNoOpt.
        The number of Q-values we should not optimize. */
    int             numQsNoOpt;
    /** Array of size numQsNoOpt.
        Which Q-values should we not optimize.
        
        This parameter can take two formats.
        
        The first is just a single number. In this case the relevant Q-parameter is
        not optimized but forced equal to its initial value.
        
        The second is of the format n=m where n and m are integers. In this case
        the n-parameter is set equal to the m-parameter. */
    char**          qsNoOpt;
    /** The maximum number of iterations for a single call to the optimizer.
        
        Note that the iteration is quite slow - this number should not be more than
        about 5 during an interactive session. */
    int             maxIter;
    /** The relative tolerance. If the minimum value changes by less than the
        tolerance multiplied by the previous minimum value, then the optimization
        will stop. */
    double          tolerance;
    /** What is the maximum time for this optimization. Note that a value of 0
        implies no upper limit is placed on the optimizer call. */
    double          maxTime;
} CrxTCdsOptionQOptModel;

/** Contains the current bond price and the settlement date for this
    price. */
typedef struct _CrxTBondPrice
{
    /** Settlement date for the trade */
    TDate           settleDate;
    /** Clean price for the trade */
    double          cleanPrice;
} CrxTBondPrice;

/** Curve describing bond repo rates for a bond.
   
    Note that conventionally, bond repo rates are used in a specific manner
    for computing the forward price of the bond. The repo rate is not only
    used for the funding costs of the bond but also for the re-investment of
    any intermediate coupons. */
typedef struct _CrxTBondRepoCurve
{
    TDate           spotSettleDate;
    /** Day count convention for the rates. This will typically be a money-market
        convention such as `ACT/360' or `ACT/365F' */
    long            dcc;
    /** Defines array size for dates, rates. */
    int             numDates;
    /** Array of size numDates.
        Dates representing the forward settlement dates. */
    TDate*          dates;
    /** Array of size numDates.
        Repo rates for the corresponding forward settlement date. Fundamentally
        these are simple rates. Use 0.04 to represent 4%. */
    double*         rates;
} CrxTBondRepoCurve;

/** Price volatility curve for a bond. */
typedef struct _CrxTBondPriceVolCurve
{
    /** Start date of volatility */
    TDate           today;
    /** Defines array size for dates, vols. */
    int             numDates;
    /** Array of size numDates.
        Dates representing the option expiry dates. */
    TDate*          dates;
    /** Array of size numDates.
        Volatilities corresponding to the option expiry dates. Use 0.04 to
        represent 4%. These are volatilities of the forward clean price. */
    double*         vols;
} CrxTBondPriceVolCurve;

/** Spread volatility curve for a bond. */
typedef struct _CrxTBondSpreadVolCurve
{
    /** Start date of volatility */
    TDate           today;
    /** Defines array size for dates, vols. */
    int             numDates;
    /** Array of size numDates.
        Dates representing the option expiry dates. */
    TDate*          dates;
    /** Array of size numDates.
        Volatilities corresponding to the option expiry dates. Use 0.4 to
        represent 40%. These are volatilities of the forward bond spread.
        Typically we measure the spread vis-a-vis a reference bond of
        similar maturity. */
    double*         vols;
} CrxTBondSpreadVolCurve;

/** Defines a bond option struck on the price of the bond. */
typedef struct _CrxTBondPriceOption
{
    TBond*          bond;
    /** Call or put. */
    long            optionType;
    TDate           exerciseDate;
    /** Payment date should be on or after the exercise date. */
    TDate           paymentDate;
    /** Clean price at which the bond option is struck. */
    double          strikePrice;
} CrxTBondPriceOption;

/** Outputs from the bond price option calculation. */
typedef struct _CrxTBondPriceOptionCalc
{
    /** Option price discounted to spot. */
    double          optionPrice;
    /** The interpolated repo rate. */
    double          repoRate;
    /** The calculated forward price. */
    double          fwdPrice;
    /** The forward premium before discounting. */
    double          fwdPremium;
} CrxTBondPriceOptionCalc;

/** Defines a bond option struck at the difference in spread between two bonds.
   
    One of the bonds will tend to be a risky bond where the bond yield is a
    combination of interest rate, credit and other components. The other bond,
    known as the reference bond, will be risk-free and its yield will only
    be dependent on interest rates.
   
    Hence the spread option is designed to eliminate the interest rate element. */
typedef struct _CrxTBondSpreadOption
{
    /** Typically the risky bond. */
    TBond*          bond;
    /** Typically the risk-free reference bond. */
    TBond*          refBond;
    /** Call or put. */
    long            optionType;
    TDate           exerciseDate;
    /** Payment date should be on or after the exercise date. */
    TDate           paymentDate;
    /** Strike spread as a decimal (e.g. use 0.02 for 200bp). */
    double          strikeSpread;
} CrxTBondSpreadOption;

/** Outputs from the bond spread option calculation. */
typedef struct _CrxTBondSpreadOptionCalc
{
    /** Option price discounted to spot. */
    double          optionPrice;
    /** Repo rate for bond. */
    double          repoRate;
    /** Repo rate for refBond. */
    double          refRepoRate;
    /** Forward price (clean) for bond. */
    double          fwdPrice;
    /** Forward price (clean) for refBond. */
    double          refFwdPrice;
    /** Forward spread corresponding to the forward prices. */
    double          fwdSpread;
    double          dpdy;
    /** The forward spread premium before conversion to price and discounting. */
    double          spreadPremium;
} CrxTBondSpreadOptionCalc;

/**
***************************************************************************
** Converts CdsOptionPayoff to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxCdsOptionPayoffToString (CrxTCdsOptionPayoff value);

/**
***************************************************************************
** Converts CdsOptionPayoff from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxCdsOptionPayoffFromString (const char* str, CrxTCdsOptionPayoff *val);

/**
***************************************************************************
** Converts CdsOptionStrikeType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxCdsOptionStrikeTypeToString (CrxTCdsOptionStrikeType value);

/**
***************************************************************************
** Converts CdsOptionStrikeType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxCdsOptionStrikeTypeFromString (const char* str, CrxTCdsOptionStrikeType *val);

/**
***************************************************************************
** Converts DistType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxDistTypeToString (long value);

/**
***************************************************************************
** Converts DistType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxDistTypeFromString (const char* str, long *val);

/**
***************************************************************************
** Converts OptionType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxOptionTypeToString (long value);

/**
***************************************************************************
** Converts OptionType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxOptionTypeFromString (const char* str, long *val);

/**
***************************************************************************
** Converts DayCountConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxDayCountConvToString (long value);

/**
***************************************************************************
** Converts DayCountConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxDayCountConvFromString (const char* str, long *val);

/**
***************************************************************************
** Constructor for CrxTQDist
***************************************************************************
*/
CrxTQDist* CrxQDistMake(
/** Defines array size for Qs.
    Number of Qs. */
int             nQs,
/** Array of size nQs.
    Q-factors. A value of Q=0 indicates no skew - i.e. lognormal.
    A value of Q=1 indicates a normal distribution.
    Typical range of values would be [-1,2] but we have seen values outside
    this range. There should be precisely 6 Q-values.
    
    If you want to use just 2Qs, then you are at liberty to set the values
    on the left equal to each other, and the values on the right likewise. */
double*         Qs,
/** Defines array size for Ds.
    Number of Ds. It should be that nDs = nQs-1 */
int             nDs,
/** Array of size nDs.
    Delta-values. These indicate the transition points between the Q-values.
    Typically the middle D will be 0.5. Values should be increasing and in
    the range 0-1 exclusive. There should be one less D-value that Q-value,
    i.e. precisely 5 values. */
double*         Ds
);

/**
***************************************************************************
** Memory allocator for CrxTQDist
***************************************************************************
*/
CrxTQDist* CrxQDistMakeEmpty(
/** Defines array size for Qs.
    Number of Qs. */
int             nQs,
/** Defines array size for Ds.
    Number of Ds. It should be that nDs = nQs-1 */
int             nDs
);

/**
***************************************************************************
** Copy constructor for CrxTQDist
***************************************************************************
*/
CrxTQDist* CrxQDistCopy(CrxTQDist* src);

/**
***************************************************************************
** Destructor for CrxTQDist
***************************************************************************
*/
void CrxQDistFree(CrxTQDist *p);

/**
***************************************************************************
** Constructor for CrxTCdsOption
***************************************************************************
*/
CrxTCdsOption* CrxCdsOptionMake(
/** Issue date of the underlying CDS. */
TDate           issueDate,
/** Maturity date of the underlying CDS. */
TDate           maturityDate,
/** Fee interval of the underlying CDS. Usual value is 3M. */
TDateInterval*  feeInterval,
/** Day count convention of the underlying CDS. Usual value is `ACT/360'. */
long            dcc,
/** Does the underlying CDS pay accrued interest on default. Usually TRUE. */
TBoolean        payAccOnDefault,
/** Does the underlying CDS represent an index instead of a single name. */
TBoolean        isIndex,
/** Loss so far - only relevant for index CDS. */
double          lossSoFar,
/** Option type - i.e. call or put. */
long            optionType,
/** Exercise date. */
TDate           exerciseDate,
/** Payment date - when the exercise payment is made. Must be on or after the
    exercise date. */
TDate           paymentDate,
/** Flag to denote that whether the option knock-out on default.
    Only relevant for single name options. */
TBoolean        koOnDefault,
/** Strike spread of the option. Measured in decimals, e.g. 100bp
    is represented as 0.01 */
double          strike,
/** Option payoff method - only relevant for index options. */
CrxTCdsOptionPayoff optionPayoff,
/** Coupon of the underlying CDS - only relevant for index options
    which pay-off using the strike method. */
double          coupon,
/** Strike type - i.e. price or spread */
CrxTCdsOptionStrikeType strikeType
);

/**
***************************************************************************
** Memory allocator for CrxTCdsOption
***************************************************************************
*/
CrxTCdsOption* CrxCdsOptionMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTCdsOption
***************************************************************************
*/
CrxTCdsOption* CrxCdsOptionCopy(CrxTCdsOption* src);

/**
***************************************************************************
** Destructor for CrxTCdsOption
***************************************************************************
*/
void CrxCdsOptionFree(CrxTCdsOption *p);

/**
***************************************************************************
** Constructor for CrxTCdsOptionCalc
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalcMake(
/** Price of the option in money terms. */
double          price,
/** At the money volatility - either implied or as input. */
double          vol,
/** Forward spread used in the option formula. */
double          fwdSpread,
/** Forward spread without taking into account index forward adjustments. */
double          fwdSpreadUnadj,
/** Risky annuity measured at the value date of the option for the forward. */
double          annuity,
/** Time to expiry (exerciseDate-today in ACT/365F) */
double          timeToExpiry,
/** Forward price at expiry */
double          fwdPrice
);

/**
***************************************************************************
** Memory allocator for CrxTCdsOptionCalc
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalcMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTCdsOptionCalc
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalcCopy(CrxTCdsOptionCalc* src);

/**
***************************************************************************
** Destructor for CrxTCdsOptionCalc
***************************************************************************
*/
void CrxCdsOptionCalcFree(CrxTCdsOptionCalc *p);

/**
***************************************************************************
** Constructor for CrxTCdsOptionQOptimization
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationMake(
/** Q-distribution */
CrxTQDist*      qdist,
/** At the money volatility */
double          vol,
/** Initial directions for search (on initialization). Final directions for
    search (on exit). */
TMatrix2D*      direction,
/** Number of iterations so far. */
int             iter,
/** Root mean square difference for the optimization. The difference is
    computed for strike volatilities. */
double          value,
/** Most recent absolute improvement in the algorithm. */
double          vdiff
);

/**
***************************************************************************
** Memory allocator for CrxTCdsOptionQOptimization
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTCdsOptionQOptimization
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationCopy(CrxTCdsOptionQOptimization* src);

/**
***************************************************************************
** Destructor for CrxTCdsOptionQOptimization
***************************************************************************
*/
void CrxCdsOptionQOptimizationFree(CrxTCdsOptionQOptimization *p);

/**
***************************************************************************
** Constructor for CrxTCdsOptionQOptModel
***************************************************************************
*/
CrxTCdsOptionQOptModel* CrxCdsOptionQOptModelMake(
/** Denotes whether we optimize the ATM volatility. Sometimes a better fit
    can be obtained by optimizing the ATM volatility, but in general we
    would expect the ATM volatility to be defined and not changed. */
TBoolean        optimizeVol,
/** Defines array size for qsNoOpt.
    The number of Q-values we should not optimize. */
int             numQsNoOpt,
/** Array of size numQsNoOpt.
    Which Q-values should we not optimize.
    
    This parameter can take two formats.
    
    The first is just a single number. In this case the relevant Q-parameter is
    not optimized but forced equal to its initial value.
    
    The second is of the format n=m where n and m are integers. In this case
    the n-parameter is set equal to the m-parameter. */
char**          qsNoOpt,
/** The maximum number of iterations for a single call to the optimizer.
    
    Note that the iteration is quite slow - this number should not be more than
    about 5 during an interactive session. */
int             maxIter,
/** The relative tolerance. If the minimum value changes by less than the
    tolerance multiplied by the previous minimum value, then the optimization
    will stop. */
double          tolerance,
/** What is the maximum time for this optimization. Note that a value of 0
    implies no upper limit is placed on the optimizer call. */
double          maxTime
);

/**
***************************************************************************
** Memory allocator for CrxTCdsOptionQOptModel
***************************************************************************
*/
CrxTCdsOptionQOptModel* CrxCdsOptionQOptModelMakeEmpty(
/** Defines array size for qsNoOpt.
    The number of Q-values we should not optimize. */
int             numQsNoOpt
);

/**
***************************************************************************
** Copy constructor for CrxTCdsOptionQOptModel
***************************************************************************
*/
CrxTCdsOptionQOptModel* CrxCdsOptionQOptModelCopy(CrxTCdsOptionQOptModel* src);

/**
***************************************************************************
** Destructor for CrxTCdsOptionQOptModel
***************************************************************************
*/
void CrxCdsOptionQOptModelFree(CrxTCdsOptionQOptModel *p);

/**
***************************************************************************
** Constructor for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceMake(
/** Settlement date for the trade */
TDate           settleDate,
/** Clean price for the trade */
double          cleanPrice
);

/**
***************************************************************************
** Memory allocator for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceCopy(CrxTBondPrice* src);

/**
***************************************************************************
** Destructor for CrxTBondPrice
***************************************************************************
*/
void CrxBondPriceFree(CrxTBondPrice *p);

/**
***************************************************************************
** Constructor for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveMake(
TDate           spotSettleDate,
/** Day count convention for the rates. This will typically be a money-market
    convention such as `ACT/360' or `ACT/365F' */
long            dcc,
/** Defines array size for dates, rates. */
int             numDates,
/** Array of size numDates.
    Dates representing the forward settlement dates. */
TDate*          dates,
/** Array of size numDates.
    Repo rates for the corresponding forward settlement date. Fundamentally
    these are simple rates. Use 0.04 to represent 4%. */
double*         rates
);

/**
***************************************************************************
** Memory allocator for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveMakeEmpty(
/** Defines array size for dates, rates. */
int             numDates
);

/**
***************************************************************************
** Copy constructor for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveCopy(CrxTBondRepoCurve* src);

/**
***************************************************************************
** Destructor for CrxTBondRepoCurve
***************************************************************************
*/
void CrxBondRepoCurveFree(CrxTBondRepoCurve *p);

/**
***************************************************************************
** Constructor for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveMake(
/** Start date of volatility */
TDate           today,
/** Defines array size for dates, vols. */
int             numDates,
/** Array of size numDates.
    Dates representing the option expiry dates. */
TDate*          dates,
/** Array of size numDates.
    Volatilities corresponding to the option expiry dates. Use 0.04 to
    represent 4%. These are volatilities of the forward clean price. */
double*         vols
);

/**
***************************************************************************
** Memory allocator for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveMakeEmpty(
/** Defines array size for dates, vols. */
int             numDates
);

/**
***************************************************************************
** Copy constructor for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveCopy(CrxTBondPriceVolCurve* src);

/**
***************************************************************************
** Destructor for CrxTBondPriceVolCurve
***************************************************************************
*/
void CrxBondPriceVolCurveFree(CrxTBondPriceVolCurve *p);

/**
***************************************************************************
** Constructor for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveMake(
/** Start date of volatility */
TDate           today,
/** Defines array size for dates, vols. */
int             numDates,
/** Array of size numDates.
    Dates representing the option expiry dates. */
TDate*          dates,
/** Array of size numDates.
    Volatilities corresponding to the option expiry dates. Use 0.4 to
    represent 40%. These are volatilities of the forward bond spread.
    Typically we measure the spread vis-a-vis a reference bond of
    similar maturity. */
double*         vols
);

/**
***************************************************************************
** Memory allocator for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveMakeEmpty(
/** Defines array size for dates, vols. */
int             numDates
);

/**
***************************************************************************
** Copy constructor for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveCopy(CrxTBondSpreadVolCurve* src);

/**
***************************************************************************
** Destructor for CrxTBondSpreadVolCurve
***************************************************************************
*/
void CrxBondSpreadVolCurveFree(CrxTBondSpreadVolCurve *p);

/**
***************************************************************************
** Constructor for CrxTBondPriceOption
***************************************************************************
*/
CrxTBondPriceOption* CrxBondPriceOptionMake(
TBond*          bond,
/** Call or put. */
long            optionType,
TDate           exerciseDate,
/** Payment date should be on or after the exercise date. */
TDate           paymentDate,
/** Clean price at which the bond option is struck. */
double          strikePrice
);

/**
***************************************************************************
** Memory allocator for CrxTBondPriceOption
***************************************************************************
*/
CrxTBondPriceOption* CrxBondPriceOptionMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTBondPriceOption
***************************************************************************
*/
CrxTBondPriceOption* CrxBondPriceOptionCopy(CrxTBondPriceOption* src);

/**
***************************************************************************
** Destructor for CrxTBondPriceOption
***************************************************************************
*/
void CrxBondPriceOptionFree(CrxTBondPriceOption *p);

/**
***************************************************************************
** Constructor for CrxTBondPriceOptionCalc
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalcMake(
/** Option price discounted to spot. */
double          optionPrice,
/** The interpolated repo rate. */
double          repoRate,
/** The calculated forward price. */
double          fwdPrice,
/** The forward premium before discounting. */
double          fwdPremium
);

/**
***************************************************************************
** Memory allocator for CrxTBondPriceOptionCalc
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalcMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTBondPriceOptionCalc
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalcCopy(CrxTBondPriceOptionCalc* src);

/**
***************************************************************************
** Destructor for CrxTBondPriceOptionCalc
***************************************************************************
*/
void CrxBondPriceOptionCalcFree(CrxTBondPriceOptionCalc *p);

/**
***************************************************************************
** Constructor for CrxTBondSpreadOption
***************************************************************************
*/
CrxTBondSpreadOption* CrxBondSpreadOptionMake(
/** Typically the risky bond. */
TBond*          bond,
/** Typically the risk-free reference bond. */
TBond*          refBond,
/** Call or put. */
long            optionType,
TDate           exerciseDate,
/** Payment date should be on or after the exercise date. */
TDate           paymentDate,
/** Strike spread as a decimal (e.g. use 0.02 for 200bp). */
double          strikeSpread
);

/**
***************************************************************************
** Memory allocator for CrxTBondSpreadOption
***************************************************************************
*/
CrxTBondSpreadOption* CrxBondSpreadOptionMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTBondSpreadOption
***************************************************************************
*/
CrxTBondSpreadOption* CrxBondSpreadOptionCopy(CrxTBondSpreadOption* src);

/**
***************************************************************************
** Destructor for CrxTBondSpreadOption
***************************************************************************
*/
void CrxBondSpreadOptionFree(CrxTBondSpreadOption *p);

/**
***************************************************************************
** Constructor for CrxTBondSpreadOptionCalc
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalcMake(
/** Option price discounted to spot. */
double          optionPrice,
/** Repo rate for bond. */
double          repoRate,
/** Repo rate for refBond. */
double          refRepoRate,
/** Forward price (clean) for bond. */
double          fwdPrice,
/** Forward price (clean) for refBond. */
double          refFwdPrice,
/** Forward spread corresponding to the forward prices. */
double          fwdSpread,
double          dpdy,
/** The forward spread premium before conversion to price and discounting. */
double          spreadPremium
);

/**
***************************************************************************
** Memory allocator for CrxTBondSpreadOptionCalc
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalcMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CrxTBondSpreadOptionCalc
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalcCopy(CrxTBondSpreadOptionCalc* src);

/**
***************************************************************************
** Destructor for CrxTBondSpreadOptionCalc
***************************************************************************
*/
void CrxBondSpreadOptionCalcFree(CrxTBondSpreadOptionCalc *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _CRX_CRXDATA_H */
