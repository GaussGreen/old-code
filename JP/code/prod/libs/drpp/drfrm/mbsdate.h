/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsdate.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Robert Lenk
 *			   Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Code version    :  1.20
 *	Extracted	:  6/20/97 at 12:42:02
 *	Last Updated	:  6/20/97 at 12:41:55
 ***************************************************************************
 *      Cashflow/accrual/reset date utility functions
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef _MBSDATE_H
#define _MBSDATE_H
extern "C" {
#include "cdate.h"
}
/************************************************************
 *  PUBLIC FUNCS
 ************************************************************/

/************************************************************
 * GetAccrualStartDateOfCurrPeriod
 * (Monthly accrual periods only)
 * Determines the accrual start date of the accrual period
 * in which given date lies
 *   Special case: if the accrual start day-of-month
 * is beyond last valid day in month (e.g., 31 in February),
 * function uses last valid day in the month
 ************************************************************/
int EXPORT GetAccrualStartDateOfCurrPeriod
   (TDate     currDate,           /* (I) any date */
    long      accrualStartDom,    /* (I) day-of-month (1-31) on which
                                   * accrual begins */
    TDate    *accrualStartDate);  /* (O) start date of current accrual pd. */

/************************************************************
 *  GetAccrualStartDateFromCfDate
 * (Monthly accrual periods only)
 *  Determines the accrual start date for the indicated cashflow
 ************************************************************/
int EXPORT GetAccrualStartDateFromCfDate
   (TDate    cashflowDate,        /* (I) cashflow date */
    long     payDelayDays,        /* (I) # days between start of next
                                   * accrual period and cashflow (e.g., 19) */
    long     accrualStartDom,     /* (I) day-of-month (1-31) on which
                                   * accrual begins */
    TDate   *accrualStartDate);   /* (O) start date of accrual pd for CF */

/************************************************************
 *  GetCfDateFromRegAccrStDate
 *  (Monthly accrual periods only)
 *  Determines the cashflow date for the specified reg. accrual start date
 ************************************************************/
int EXPORT GetCfDateFromRegAccrStDate
   (TDate    accrualStartDate,    /* (I) start date of regular accrual pd */
    long     payDelayDays,        /* (I) # days between start of next
                                   * accrual period and cashflow (e.g., 19) */
    TDate   *cashflowDate);       /* (O) cashflow date */

/************************************************************
 * GetNextEffResetDate
 * Computes next eff reset date from curr reset date;
 * can accept 0 (i.e., init cpn period) as curr reset date
 * OK to use same var for input & output date
 ************************************************************/
int EXPORT GetNextEffResetDate
   (TDate    effResetDate,        /* (I) curr eff reset date */
    TDate    firstEffResetDate,   /* (I) first effective reset date of 
                                   * instrument; can use 0 if monthly */
    long     monthsBetweenResets, /* (I) # months between resets (>= 1)
                                   * accrual begins */
    TDate   *nextEffResetDate);   /* (O) next effective reset date */

/************************************************************
 *  GetEffResetDateFromAccrualStartDate
 *  (Monthly accrual periods only)
 *  Determines the effective reset date associated with the given
 *  accrual start date; if effective reset of this accrual period
 *  comes BEFORE first effective reset of instrument, then
 *  return 0 (denoting accrual on inital coupon)
 *    NB: can handle resets which occur monthly, or at
 *  intervals of an integral number of months (semi, annual, etc)
 ************************************************************/
int EXPORT GetEffResetDateFromAccrualStartDate
   (TDate    accrualStartDate,    /* (I) start date of accrual pd */
    long     effResetLkbkMons,    /* (I) effective reset (typ. 2 bus. days
                                   * after reset lookup) occurs this
                                   * # months before start of accrual period */
    long     effResetDom,         /* (I) day-of-month (1-31) of effective
                                   * reset (e.g., 1 for 1st of month) */
    TDate    firstEffResetDate,   /* (I) first effective reset date of 
                                   * instrument; use 0 none or monthly */
    long     monthsBetweenResets, /* (I) # months between resets (>= 1) */
    TDate   *effResetDate);       /* (O) effective reset date */

/************************************************************
 *  GetEffResetDateFromCfDate
 *  (Monthly accrual periods only)
 *  Determines the reset date associated with the given
 *  cashflow date; if reset of this cashflow date
 *  comes BEFORE first reset of instrument, then
 *  return 0 (denoting cashflow determined by inital coupon)
 *    NB: can handle resets which occur monthly, or at
 *  intervals of an integral number of months (semi, annual, etc)
 ************************************************************/
int EXPORT GetEffResetDateFromCfDate
   (TDate    cashflowDate,        /* (I) cashflow date */
    long     payDelayDays,        /* (I) # days between start of next
                                   * accrual period and cashflow (e.g., 19) */
    long     accrualStartDom,     /* (I) day-of-month (1-31) on which
                                   * accrual begins */
    long     effResetLkbkMons,    /* (I) effective reset (typ. 2 bus. days
                                   * after reset lookup) occurs this
                                   * # months before start of accrual period */
    long     effResetDom,         /* (I) day-of-month (1-31) of effective
                                   * reset (e.g., 1 for 1st of month) */
    TDate    firstEffResetDate,   /* (I) first effective reset date of 
                                   * instrument; use 0 none or monthly */
    long     monthsBetweenResets, /* (I) # months between resets (>= 1) */
    TDate   *effResetDate);       /* (O) effective reset date */

/************************************************************
 *  GetAccrualStartDateFromEffResetDate
 *  (Monthly accrual periods only)
 *  Determines the (first) accrual start date
 *  corresponding to the given reset date 
 ************************************************************/
int EXPORT GetAccrualStartDateFromEffResetDate
   (TDate    effResetDate,        /* (I) date of effective reset */
    long     effResetLkbkMons,    /* (I) effective reset occurs this # 
                                   * months before start of accrual period */
    long     accrualStartDom,     /* (I) day-of-month (1-31) of start
                                   * of accrual period */
    TDate    effOrigDate,         /* (I) orig. date of MBS */
    TDate   *accrualStartDate);   /* (O) accrual start date */


/************************************************************
 *  GetCfDateFromEffResetDate
 *  (Monthly accrual periods only)
 *  Determines the (first) cashflow date
 *  corresponding to the given reset date 
 *  Note: can accept reset date of 0 (denotes
 *  initial cpn period of MBS)
 ************************************************************/
int EXPORT GetCfDateFromEffResetDate
   (TDate    effResetDate,        /* (I) date of effective reset */
    long     effResetLkbkMons,    /* (I) effective reset occurs this # 
                                   * months before start of accrual period */
    long     accrualStartDom,     /* (I) day-of-month (1-31) of start
                                   * of accrual period */
    TDate    effOrigDate,         /* (I) orig. date of MBS */
    long     payDelayDays,        /* (I) pay delay, in days */
    TDate   *cashflowDate);       /* (O) cashflow date */


/************************************************************
 * FindMatchingDate
 * Attempts to find given date in date list;
 * returns index (in date list) of this date,
 * or -1 if not found
 ************************************************************/
int EXPORT FindMatchingDate
   (TDate       currDate,       /* (I) date to find in list */
    TDateList  *dateList,       /* (I) date list to search */
    long       *dateIdx);       /* (O) index of this date in datelist */

/************************************************************
 * FindClosestDOM
 * Finds date closest to initial date that has desired day-of-mon
 * NB: OK to use same date for input & output
 ************************************************************/
int EXPORT FindClosestDOM
   (TDate     date,             /* (I) initial guess for date */
    long      dom,              /* (I) day-of-month desired */
    TDate    *closestDate);     /* (O) closest date w/desired d-o-m */

/************************************************************
 * FindCloseDate
 * Finds loose-tolerance matchine date in date list;
 * if none found, returns input date
 ************************************************************/
int EXPORT FindCloseDate
   (TDate      date,            /* (I) initial guess for date */
    TDateList *dateList,        /* (I) list of dates */
    long       maxDayDiff,      /* (I) Max # days difference between dates
                                 * that's still considered a match */
    TBoolean  *foundMatch,      /* (O) */
    TDate     *matchingDate);   /* (O) */

/************************************************************
 * CheckDates
 * Checks that dates in list are ascending, w/o repeats
 ************************************************************/
int EXPORT CheckDates
   (long      numDates,         /* (I) # dates in array */
    TDate    *dates,            /* (I) array of dates to check */
    char     *labelStr);        /* (I) name of this array, for error msg */



#endif
