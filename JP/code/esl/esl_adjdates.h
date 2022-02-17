/****************************************************************************/
/*      Header for ESL adjustable and flexible dates.                       */
/****************************************************************************/
/*      ESL_ADJDATES_H                                                      */
/****************************************************************************/

#ifndef ADJ_ADJDATES_H
#define ADJ_ADJDATES_H

#include "esl_macros.h"
#include "esl_types.h"


/*********************************************************************
 *
 *  The general idea here is to define component structures that
 *  can be used in a generic date adjustment layer.
 *  These structures are adapter structures that point to the 
 *  relevant fields in the underlying deal structure(s), so adjusting 
 *  the component structures automatically adjusts the underlying deal
 *  structure(s)
 *
 *  NOTE:  This only includes logic that is self contained within ESL
 *         and library specific adjustment functionality (for example
 *         using Fwd_FX) must be part of the specific library code
 *********************************************************************/


/*****************
 * ENUMS
 *****************/

/* Reset dates are defined in Kapital as offsets to physical dates (either
   accrual start, accrual end or payment dates
   Kapital supplies us with which physical date the offset is applied from,
   and from this information we can decide if we are moving the reset to
   in advance or arrears.
   Some advanced streams also support absolute dates or date schedules that
   are independent of the accrual/payment schedules, and hence it is not likely
   to be an obvious decision between moving it to in advance or arrears.  No 
   effort is made to adjust these complicated types of reset schedules and 
   the reset type in this case is undefined.
*/
typedef enum 
{
    InAdvance,
    InArrears,
    UnDefined

} ADJ_RESET_TYPE;

/* defines coupon type with regard to past resets - whether it supports a single
   entry past reset or an array of multiple past reset dates and rates */
typedef enum
{
    Single,
    Multiple

} ADJ_PAST_FIXING_TYPE;


/**********************************
 * COUPON STRUCTURES AND FUNCTIONS
 **********************************/

/* past reset may be interest rate or FX or equity */
typedef struct
{
    long   *ResetDate;
    long   *ResetEffDate;
    double *ResetRate;

} ADJ_PAST_RESET;


/* generic coupon adapter structure */
typedef struct
{
 
    /* coupon data */
    int  *NbCoupons;
    long *AccStartDates;
    long *AccEndDates;
    long *PaymentDates;
    long *ResetDates;
    long *ResetEffDates;

    ADJ_RESET_TYPE ResetType;

    /* reset data */
    int  IndexMaturity;
    char IndexDayCount;
    char IndexFrequency;
 

    ADJ_PAST_FIXING_TYPE PastFixingType;
    int                  *NbPastFixings;
    ADJ_PAST_RESET       PastResets;  /* handles both single and multiple fixing cases */
    
    /* adjustment factors */
    double YieldRatio[MAXNBDATE];  /* = Y(deal reset date)/Y(model reset date) */
    double ConstantReset[MAXNBDATE]; /* used when constant future fixing required */

} ADJ_COUPON_DATA;


/* calculate and apply the date adjustment for flexible reset dates */
int ADJ_AdjustCouponDates(
    ADJ_COUPON_DATA*        coupon_data,            /* (I/O) coupon structure     */
    long                    valueDate,              /* (I) value date             */
    T_CURVE const*    zeroCurve,              /* (I) zero curve for reset   */
    int                     includeValueDateEvents, /* (I) true/false */
    const char*             errorMsg);              /* (I) err msg help string    */


/**********************************
 * OPTION STRUCTURES AND FUNCTIONS
 **********************************/

/* generic option structure */
typedef struct
{
    int  *NbExercise;
    long *ExerDates;
    long *NotifDates;
    long *NotifEffDates;

    double *Strikes;

} ADJ_OPTION_DATA;


/* apply date adjustment for flexible notification dates */
int ADJ_AdjustOptionDates(
    ADJ_OPTION_DATA  *option_data,    /* (I/O) adapter structure           */
    long             valueDate);      /* (I) value date                    */


/********************
 * UTILITY FUNCTIONS
 ********************/

/* same as library version except it doesn't print out error
   message to DR_Error - simply returns SUCCES or FAILURE */
int ADJ_DrDatesInSet(
    int        NbDates1,      /* (I) Number of dates to check   */
    long      *Dates1,        /* (I) Dates to check             */
    int        NbDates2,      /* (I) Number of org dates        */ 
    long      *Dates2);       /* (I) Org dates                  */



#endif

