/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsdate.c
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Robert Lenk
 *			   Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Code version    :  1.20
 *	Extracted	:  6/20/97 at 12:42:02
 *	Last Updated	:  6/20/97 at 12:41:54
 ***************************************************************************
 *      Cashflow/accrual/reset date utility functions
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#include <math.h>               /* exp */
extern "C" {
#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "cdate.h"              /* DayCountFraction */
#include "ldate.h"              /* DayCountFraction */
#include "datelist.h"           /* GtoNewDateList */
#include "cerror.h"             /* GtoErrMsg */
#include "cmemory.h"            /* MALLOC */
#include "convert.h"            /* GtoFormatDate */
#include "yearfrac.h"           /* GtoFormatDayCountConv */
}

#include "mbsconst.h"

#include "mbsdate.h"            /* Prototype consistency */



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
    TDate    *accrualStartDate)   /* (O) start date of current accrual pd. */
{
    F_INIT("GetAccrualStartDateOfCurrPeriod");
    long    currDom;

    /* Check inputs */
    XONTRUE( accrualStartDom < 1 || accrualStartDom > 31,
            "Invalid accrualStartDom" );

    /* First, determine the accrual start date in curr month */
    *accrualStartDate = currDate;
    XONFAIL( SetDOM(accrualStartDom,TRUE,accrualStartDate) );

    /* If currDate before this date, then accrual starts 
     * in previous month */
    XONFAIL( GetDOM(currDate, &currDom) );
    if( currDom < accrualStartDom )
    {
        XONFAIL( NxtMth(*accrualStartDate, -1, accrualStartDate) );
    }        

    F_END;
}


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
    TDate   *accrualStartDate)    /* (O) start date of accrual pd for CF */
{
    F_INIT("GetAccrualStartDateFromCfDate");

    /* move back to start of the accrual period AFTER the one we want */
    *accrualStartDate = cashflowDate - payDelayDays;
    /* then move back to prev accrual start date (1mon before),
     * which is the one we want */
    XONFAIL( NxtMth(*accrualStartDate, -1, accrualStartDate) );
    /* then (in case pay delay days is approx) move this to closest
     * date with desired day-of-month */
    XONFAIL(FindClosestDOM(*accrualStartDate,
                           accrualStartDom,
                           accrualStartDate) );

    F_END;
}


/************************************************************
 *  GetCfDateFromRegAccrStDate
 *  (Monthly accrual periods only)
 *  Determines the cashflow date for the specified reg. accrual start date
 ************************************************************/
int EXPORT GetCfDateFromRegAccrStDate
   (TDate    accrualStartDate,    /* (I) start date of regular accrual pd */
    long     payDelayDays,        /* (I) # days between start of next
                                   * accrual period and cashflow (e.g., 19) */
    TDate   *cashflowDate)        /* (O) cashflow date */
{
    F_INIT("GetCfDateFromRegAccrStDate");

    /* (monthly accrual) move fwd to start of next accrual period */
    XONFAIL( NxtMth(accrualStartDate, 1, cashflowDate) );
    /* then move ahead by delay days to reach cashflow date */
    *cashflowDate += payDelayDays;

    F_END;
}


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
    TDate   *nextEffResetDate)    /* (O) next effective reset date */
{
    F_INIT("GetNextEffResetDate");
    TDate    tmpDate;

    /* if curr date is 0 (init cpn period), use "first" date */
    if( effResetDate IS 0 )
    {
        tmpDate = firstEffResetDate;
    }
    else
    {
        XONFAIL( NxtMth(effResetDate,
                        monthsBetweenResets,
                        &tmpDate) );
    }
    *nextEffResetDate = tmpDate;

    F_END;
}



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
    TDate   *effResetDate)        /* (O) effective reset date */
{
    F_INIT("GetEffResetDateFromAccrualStartDate");
    TDate         computedEffResetDate;
    TDate         nextEffResetDate;

    /* Check inputs */
    XONTRUE( effResetLkbkMons < 0, "Cannot have negative lookback period" );
    XONTRUE( monthsBetweenResets < 1,
            "Invalid monthsBetweenResets" );

    /* If interval between resets > 1 month,
     * MUST specify a first reset date, to define the offset
     * of the reset cycle
     */
    XONTRUE(monthsBetweenResets > 1 &&
            firstEffResetDate IS 0,
            "For resets occurring > monthly, MUST supply first reset date" );

    /* compute reset date assumuming monthly resets */
    XONFAIL( NxtMth(accrualStartDate,
                    -effResetLkbkMons, 
                    &computedEffResetDate) );
    XONFAIL( SetDOM(effResetDom, TRUE, &computedEffResetDate) );

    /* If computed reset comes before first reset, 
     * return 0 to denote accrual on initial coupon
     */
    if( computedEffResetDate < firstEffResetDate )
    {
        *effResetDate = 0;
    }
    
    /* (else computed reset on/after first reset) */
    else {

        /* For monthly resets, or if computed reset
         * date happens to be first reset, we're done */
        if(monthsBetweenResets IS 1 ||
           computedEffResetDate IS firstEffResetDate )
        {
            *effResetDate = computedEffResetDate;
        }
        /* else (resets occur at longer intervals than monthly
         *       AND computed reset is after first reset):
         * -->search forward to find latest reset
         * on/before our computed reset date 
         */
        else
        {
            nextEffResetDate = firstEffResetDate;
            while( nextEffResetDate <= computedEffResetDate )
            {
                *effResetDate = nextEffResetDate;
                XONFAIL( NxtMth(nextEffResetDate, 
                                monthsBetweenResets, 
                                &nextEffResetDate) );
                /* always ensure day-of-month is correct */
                XONFAIL( SetDOM(effResetDom, TRUE, &nextEffResetDate) );
            }
        }

    }    /* (else computed reset on/after first reset) */

    F_END;
}


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
    TDate   *effResetDate)        /* (O) effective reset date */
{
    F_INIT("GetEffResetDateFromCfDate");
    TDate  accrualStartDate;

    /* get accrual start date for this CF */
    XONFAIL( GetAccrualStartDateFromCfDate(cashflowDate,
                                           payDelayDays,
                                           accrualStartDom,
                                           &accrualStartDate) );
    /* then get eff. reset date of this cashflow */
    XONFAIL( GetEffResetDateFromAccrualStartDate(accrualStartDate,
                                                 effResetLkbkMons,
                                                 effResetDom,
                                                 firstEffResetDate,
                                                 monthsBetweenResets,
                                                 effResetDate) );

    F_END;
}


/************************************************************
 *  GetAccrualStartDateFromEffResetDate
 *  (Monthly accrual periods only)
 *  Determines the (first) accrual start date
 *  corresponding to the given reset date 
 *  Note: can accept reset date of 0 (denotes
 *  initial cpn period of MBS)
 ************************************************************/
int EXPORT GetAccrualStartDateFromEffResetDate
   (TDate    effResetDate,        /* (I) date of effective reset */
    long     effResetLkbkMons,    /* (I) effective reset occurs this # 
                                   * months before start of accrual period */
    long     accrualStartDom,     /* (I) day-of-month (1-31) of start
                                   * of accrual period */
    TDate    effOrigDate,         /* (I) orig. date of MBS */
    TDate   *accrualStartDate)    /* (O) accrual start date */
{
    F_INIT("GetAccrualStartDateFromEffResetDate");

    /* Check inputs */
    XONTRUE( effResetLkbkMons < 0, "Cannot have negative lookback period" );
    XONTRUE( accrualStartDom < 1 || accrualStartDom > 31, 
            "Invalid accrual start day-of-month" );
    
    /* special case: if reset date = 0, accrual start is
     * the origination date of the bond */
    if( effResetDate IS 0 )
    {
        *accrualStartDate = effOrigDate;
        XONFAIL( SetDOM(1,TRUE,accrualStartDate) );
    }
    else
    {
        /* Move (fwd) from eff. reset date to month of accrual start */
        XONFAIL( NxtMth(effResetDate,effResetLkbkMons,accrualStartDate) );
        /* But accrual starts on a specific day-of-month */
        XONFAIL( SetDOM(accrualStartDom,TRUE,accrualStartDate) );
    }

    F_END;
}


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
    TDate   *cashflowDate)        /* (O) cashflow date */
{
    F_INIT("GetCfDateFromEffResetDate");
    TDate    regAccStDate;

    /* first, get first accrual start date */
    XONFAIL(GetAccrualStartDateFromEffResetDate
            (effResetDate,
             effResetLkbkMons,
             accrualStartDom,
             effOrigDate,
             &regAccStDate) );

    /* then get cf date for this regular acc date */
    XONFAIL(GetCfDateFromRegAccrStDate
            (regAccStDate,
             payDelayDays,
             cashflowDate) );

    F_END;
}



/************************************************************
 * FindMatchingDate
 * Attempts to find given date in date list;
 * returns index (in date list) of this date,
 * or -1 if not found
 ************************************************************/
int EXPORT FindMatchingDate
   (TDate       currDate,       /* (I) date to find in list */
    TDateList  *dateList,       /* (I) date list to search */
    long       *dateIdx)        /* (O) index of this date in datelist */
{
    F_INIT("FindMatchingDate");
    long    idx;

    /* reset output */
    *dateIdx = -1;

    XONTRUE(dateList IS NULL, "Failed to provide datelist" );

    for(idx=0; idx<dateList->fNumItems; idx++)
    {
        if( currDate IS dateList->fArray[idx] )
        {
            *dateIdx = idx;
            break;
        }
    }

    F_END;
}



/************************************************************
 * FindClosestDOM
 * Finds date closest to initial date that has desired day-of-mon
 * NB: OK to use same date for input & output
 ************************************************************/
int EXPORT FindClosestDOM
   (TDate     date,             /* (I) initial guess for date */
    long      dom,              /* (I) day-of-month desired */
    TDate    *closestDate)      /* (O) closest date w/desired d-o-m */
{
    F_INIT("FindClosestDOM");
    TDate   priorMonDate;
    TDate   currMonDate;
    TDate   nextMonDate;
    long    dPrior;
    long    dCurr;
    long    dNext;
    long    dMin;

    /* dates w/this day-of-month on curr & surrounding months */
    XONFAIL( NxtMth(date,-1,&priorMonDate) );
    XONFAIL( SetDOM(dom,TRUE,&priorMonDate) );
    currMonDate = date;
    XONFAIL( SetDOM(dom,TRUE,&currMonDate) );
    XONFAIL( NxtMth(date,1,&nextMonDate) );
    XONFAIL( SetDOM(dom,TRUE,&nextMonDate) );

    /* # actual days diff */
    dPrior = abs((int)(date-priorMonDate));
    dCurr = abs((int)(date-currMonDate));
    dNext = abs((int)(date-nextMonDate));
    
    /* find closest */
    dMin = dPrior;
    if( dCurr < dMin )
    {
        dMin = dCurr;
    }
    if( dNext < dMin )
    {
        dMin = dNext;
    }

    /* pick closest, preferably curr month if a tie */
    if(dCurr IS dMin)
    {
        *closestDate = currMonDate;
    }
    else if(dPrior IS dMin)
    {
        *closestDate = priorMonDate;
    }
    else
    {
        *closestDate = nextMonDate;
    }

    F_END;
}



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
    TDate     *matchingDate)    /* (O) */
{
    F_INIT("FindCloseDate");
    long    idx;

    *matchingDate = date;
    *foundMatch = FALSE;

    XONTRUE(maxDayDiff <= 0,
            "Cannot have max # days diff be < 1");

    for(idx=0; idx<dateList->fNumItems; idx++)
    {
        if( abs((int)(date - dateList->fArray[idx])) <= maxDayDiff )
        {
            *matchingDate = dateList->fArray[idx];
            *foundMatch = TRUE;
            break;
        }
    }

    F_END;
}



/************************************************************
 * CheckDates
 * Checks that dates in list are ascending, w/o repeats
 ************************************************************/
int EXPORT CheckDates
   (long      numDates,         /* (I) # dates in array */
    TDate    *dates,            /* (I) array of dates to check */
    char     *labelStr)         /* (I) name of this array, for error msg */
{
    F_INIT("CheckDates");
    long   iDt;

    if(numDates > 0)
    {
        XONTRUE(dates IS NULL,
                "Failed to supply any dates" );
        for(iDt=0; iDt<numDates-1; iDt++)
        {
            if(dates[iDt] IS dates[iDt+1])
            {
                sprintf(outmesg,
                        "Error: repeated date found in %s: %s\n",
                        labelStr,
                        GtoFormatDate(dates[iDt]));
                XONTRUE(TRUE,outmesg);
            }
            else if(dates[iDt] > dates[iDt+1])
            {
                sprintf(outmesg,
                        "Error: dates # %ld,%ld have wrong order in %s:\n%s, %s\n",
                        iDt, iDt+1,
                        labelStr,
                        GtoFormatDate(dates[iDt]),
                        GtoFormatDate(dates[iDt+1]));
                XONTRUE(TRUE,outmesg);
            }
        }
    }

    F_END;
}

