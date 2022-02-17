#ifndef _drtypes_h
#define _drtypes_h

#ifdef  __cplusplus
extern "C" {
#endif

typedef long drdate;


/** Coupon payment stub types. Describe how coupon payments are made if they 
    occure on non coupon payment dates.

    NOTE: The order of enum types is important for the simple translation
          to DR and ALib types to work.
*/
typedef enum
{
    NONE = 0,   /**< Receive no coupon when selling. Receive full coupon on 
                     payment date when buying. This is dirty price of a bond. */
    BOND,       /**< Receive partial coupon with accrual period from accrual
                     start to current date when selling. Pay partial coupon 
                     with accrual period from accrual start to current date,
                     receive full coupon on payment date. Coupon has to be 
                     known in advance.  This is clean price of a bond.        */
    SIMPLE,     /**< Receive partial coupon with accrual period from current
                     date to accrual end when buying.                         */
    PAR         /**< Receive par rate from current date to payment date and 
                     accrued over the remaining period.  This will ensure that
                     the floating leg is priced at par.                       */

} KStubType;
 
/** String aray for quick translation from KStubType to name
*/
extern char const* KStubTypeStrings[];


/** Stub location. Partial coupon payment results from non integral number 
    for coupon periods per tenor. While determining coupon payment dates we
    can start from the maturity date going backward. This will result in 
    shorter first coupon period or the 'front' stub. Alternatively, we start 
    from the inception date and move forward to maturity. This will result in 
    shorter last coupon period or the 'back' stub.  "Long" stub refers to
    a short stub plus one regular coupon interval.

    NOTE: The order of enum types is important for the simple translation
          to DR and ALib types to work.
*/
typedef enum
{
    SHORT_FRONT = 0,
    SHORT_BACK,
    LONG_FRONT,
    LONG_BACK

} KStubLocation;

/** String aray for quick translation from KStubLocation to name
*/
extern char const* KStubLocationStrings[];


/** Day count convention

    NOTE: The order of enum types is important for the simple translation
          to DR and ALib types to work.
*/
typedef enum
{
    DCC_30_360 = 0,
    DCC_ACT_360,
    DCC_ACT_365,
    DCC_ACT_ACT

} KDcc;

/** String aray for quick translation from KDcc to name
*/
extern char const* KDccStrings[];

typedef enum
{
    OPT_CALL,
    OPT_PUT
} KOptType;

/** String array for quick translation from KOptType to name
 */
extern char const* KOptTypeStrings[];


typedef enum {
    OPT_PRICE,
    OPT_DELTA,
    OPT_GAMMA,
    OPT_VEGA,
	OPT_EXERCISE_PROBABILITY /**<This is N(d2) for a call and N(-d2) for a put */
} KOptResult;

KOptType
Char2OptType(char optionType);

/** Money market basis

    NOTE: The order of enum types is important for the simple translation
          to DR and ALib types to work.
*/
typedef enum
{
    MMB_ACT = 0,
    MMB_360,
    MMB_365

} KMMBasis;

/** String aray for quick translation from KDcc to name
*/
extern char const* KMMBasisStrings[];



/** Probability distribution type

    NOTE: The order of enum types is important for the simple translation
          to DR and ALib types to work.
*/
typedef enum
{
    NORMAL = 0,
    LOGNORMAL

} KProbDistType;

/** String aray for quick translation from KProbDistType to name
*/
extern char const* KProbDistTypeStrings[];



/** Frequency type. Set of well defined frequencies. For other intervals
    one must use DateInterval type.

    NOTE: The order of enum types is important for the simple translation
          to DR and ALib types to work.
*/
typedef enum
{
    ANNUAL = 0,
    SEMI_ANNUAL,
    QUARTERLY,
    IMM,
    MONTHLY,
    WEEKLY,
    DAILY

} KFrequency;


/** String aray for quick translation from KFrequency to name
*/
extern char const* KFrequencyStrings[];



/** Knockout type.
*/
typedef enum
{
    KNOCK_IN,
    KNOCK_OUT
} KKnockType;

/** Knockout range type.
*/
typedef enum
{
    KNOCK_INSIDE,
    KNOCK_OUTSIDE
} KKnockRangeType;



#ifdef  __cplusplus
}
#endif


#endif
