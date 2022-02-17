/****************************************************************************/
/*      Do all conversions necessary to use the dev-kit                     */
/*      conversion, ...                                                     */
/****************************************************************************/
/*      DEVKITCONV.c                                                        */
/****************************************************************************/

/*
$Header:
*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "esl_devkitconv.h"
#include "esl_date.h"
#include "esl_error.h"


/**Make these long because they're compared against TDates, which
   are also long. Since this is a time-critical part of the routine,
   we want to make it as fast as possible.
 */
static long  leapCumDays[] = {
    -1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};

static long  cumDays[] = {
    -1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364};

static long leapDays[] = {
    0, 31,  29,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31};
    /* JAN  FEB  MAR  APR  MAY  JUN  JUL  AUG  SEP  OCT  NOV  DEC */

static long  days[] = {
    0, 31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31};
    /* JAN  FEB  MAR  APR  MAY  JUN  JUL  AUG  SEP  OCT  NOV  DEC */

/*****  TDate2DRDate  ***************************************************/
/**
 *      Converts an alib TDate to a DR Date (YYYYMMDD)
 */
int TDate2DRDate(TDATE    date    /** (I) Alib date format */ 
                ,long    *DrDate  /** (O) YYYYMMDD         */
		)
{
    int    status = FAILURE;
    long   day;
    long   month;
    long   year = GTO_TDATE_BASE_YEAR;
    long   fourYearBlocks;
    long   count;
    long  *cumDaysp;

    if (DrDate == NULL)
    {
        DR_Error("TDate2DRDate: output date ptr is NULL.");
        goto RETURN;
    }

    if (date < 0)
    {
        DR_Error("TDate2DRDate: TDate is < 0.");
        goto RETURN;
    }

    /* Get year
     */
    while (date >= DAYS_IN_400_YEARS)
    {
        date -= DAYS_IN_400_YEARS;
        year += 400;
    }

    /* Go through this loop at most 3 times so that Dec 31 in the
     * year 2000, 2400, etc doesn't get moved to the year 2001, 2401.
     */
    for (count = 3; date >= DAYS_IN_100_YEARS && count > 0; count--)
    {
        date -= DAYS_IN_100_YEARS;
        year += 100;
    }

    /* Dont have to make sure we go through at most 24 times since
     * the last 4 years (of 100) has *less* (or the same number of)
     * days than the other groups of 4 years.
     */
    if (date >= DAYS_IN_4_YEARS)
    {
        fourYearBlocks = date/DAYS_IN_4_YEARS;
        date -= MULTIPLY_BY_1461(fourYearBlocks);
        year += (int)fourYearBlocks << 2;   /* Multiply by 4 */
    }

    /* Go through this loop at most 3 times so that Dec 31 in a leap
     * year does not get moved to the next year.
     */
    for (count = 3; date >= DAYS_IN_1_YEAR && count > 0; count--)
    {
        date -= DAYS_IN_1_YEAR;
        year += 1;
    }

    /* Get month and date
     */

    /* date/32 is a good lower bound for month. */
    month = (date >> 5) + 1;

    if (Isleap(year))
        cumDaysp = leapCumDays + month;
    else
        cumDaysp = cumDays + month;

    /* There is an extra increment and decrement of cumDaysp here, but
       it's necessary in order to set month correctly. */
    for ( ; date > *cumDaysp; month++)
       cumDaysp++;
    day = date - *(--cumDaysp);

    *DrDate = year * 10000 + month * 100 + day;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("TDate2DRDate: failed.");
    }
    return (status);

} /* TDate2DRDate */


/*****  DRDate2TDate  ***************************************************/
/**
 *      Converts a DR Date (YYYYMMDD) to an alib TDate
 */
int  DRDate2TDate(long     DrDate  /** (I) YYYYMMDD  */
                 ,TDATE   *odate   /** (O) Alib date */
		)
{
    int      status = FAILURE;

    TDATE    date = 0;               /* Assign to odate at end */
    long     fourYearBlocks;
    long     DrYear, DrMonth, DrDay;
    long     year, month, day;
    int      isLeapYear;

    if (odate == NULL)
    {
        DR_Error("DRDate2TDate: output date ptr is NULL.");
        goto RETURN;
    }

    Dsplit(DrDate, &DrMonth, &DrDay, &DrYear);
    year  = DrYear;
    month = DrMonth;
    day   = DrDay;

    year = year - GTO_TDATE_BASE_YEAR;  /* Avoid ptr derefs*/
    isLeapYear = Isleap(DrYear);

    /* Make sure day is in range.
     */
    if (day >= 1 && day <= 28)
       /*EMPTY*/;                      /* Guaranteed to be OK */
                    /* Avoid doing check below */
    else if (day < 1 ||
        (isLeapYear ? day > leapDays[month] : day > (days[month])))
    {
        DR_Error("DRDate2TDate: Invalid DrDate.");
        goto RETURN;
    }

    /* Make sure month and year are in range.
     */
    if (month < 1 || month > GTO_MONTHS_PER_YEAR ||
        DrYear < GTO_TDATE_BASE_YEAR)
    {
        DR_Error("DRDate2TDate: DrDate is out of range.");
        goto RETURN;
    }


    /* Take years into account
     */
    while (year >= 400)
    {
        year -= 400;
        date += DAYS_IN_400_YEARS;
    }

    while (year >= 100)
    {
        year -= 100;
        date += DAYS_IN_100_YEARS;
    }

    if (year >= 4)
    {
        fourYearBlocks = (long)(year>>2);       /* Divide by 4 */
        year -= (int)(fourYearBlocks<<2);       /* Multiply by 4 */
        date += MULTIPLY_BY_1461(fourYearBlocks);
    }

    while (year >= 1)
    {
        year -= 1;
        date += DAYS_IN_1_YEAR;
    }

    if (isLeapYear)
        date += leapCumDays[month-1] + day;
    else
        date += cumDays[month-1] + day;

    *odate = date;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("DRDate2TDate: failed.");
    }
    return (status);

} /* DRDate2TDate */


/*****  IsPtrNull  ******************************************************/
/**
*       Returns: TRUE if ptr is NULL or if both 1) type is valid and 2) 
*                ptr points to the array {1,0} where both elements are 
*                represented in their corresponding type;
*                FALSE otherwise.
*       WARNING: if pointer being tested intentionally holds {1,0}, or type
*                is invalid, then the function will return TRUE
*/
int    IsPtrNull (DEVKIT_TYPE   type    /** (I) Type    */
                 ,void         *ptr     /** (I) Pointer */
		)
{
    if (ptr == NULL) return (TRUE);

    switch (type)
    {
        case DEVKIT_INT:
        {
            int *ptrL = (int *)ptr;
            if ( (ptrL[0] == 1) && (ptrL[1] == 0) ) 
                 return (TRUE);
            else return (FALSE);
        }
        case DEVKIT_LONG:
        {
            long *ptrL = (long *)ptr;            
            if ( (ptrL[0] == 1L) && (ptrL[1] == 0L) ) 
                 return (TRUE);
            else return (FALSE);
        }
        case DEVKIT_DOUBLE:
        {
            double *ptrL = (double *)ptr;
            if ( (IS_EQUAL(ptrL[0],1.)) && (IS_EQUAL(ptrL[1],0.)) ) 
                 return (TRUE);
            else return (FALSE);
        }
        case DEVKIT_CHAR:
        {
            char *ptrL = (char *)ptr;
            if ( (ptrL[0] == 1) && (ptrL[1] == '\0') ) /* platform dependent ??? */
                 return (TRUE);
            else return (FALSE);
        }
        case DEVKIT_TDATE:
        {
            TDATE *ptrL = (TDATE *)ptr;
            if ( (ptrL[0] == (TDATE)1) && (ptrL[1] == (TDATE)0) ) 
                 return (TRUE);
            else return (FALSE);
        }
        default:
        {
            return (FALSE);
        }
    }  /* switch */

}  /* IsPtrNull */

