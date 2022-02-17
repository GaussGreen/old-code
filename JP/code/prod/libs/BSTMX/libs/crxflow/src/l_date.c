/****************************************************************************/
/*      Do all the manipulations on dates: settlement, next business day,   */
/*      conversion, ...                                                     */
/****************************************************************************/
/*      DATE.c                                                              */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <string.h>

#include "l_date.h"

#define HOLIDAY_FILE "HOLIDAY"
#define OKAY 0
#define ERR -1

long noleap[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
long leap[13] =  {0,31,29,31,30,31,30,31,31,30,31,30,31};
long cumdays[12] = {0,31,59,90,120,151,181,212,243,273,304,334};


/*********************************************************************************
 *  Checks for validity of YYYYMMDD formatted dates  
 *  returns 0 if all is well                        
 *  returns 1 if bad year                           
 *  returns 2 if bad month                          
 *  returns 3 if bad day                            
 *
 ********************************************************************************/
int Dateok(long date)
{

    long *daysin;
    long mm,dd,yy;

    Dsplit(date,&mm,&dd,&yy);
    if (yy < 1600 || yy > 3000)
        return (1); /* bad year */
    if (mm < 1 || mm > 12)
        return (2); /* bad month */
    /* see if leap year */
    daysin = Isleap(yy) ? leap : noleap;
    if (dd < 1 || dd > daysin[mm])
        return (3); /* bad month */
    /* all is well if we are here */
    return (0);
}

/*********************************************************************************
 * split date into month, day, year
 *
 ********************************************************************************/

void Dsplit(long   date_i,                /* (I) london date                    */
            long   *mm_o,                 /* (O) month                          */
            long   *dd_o,                 /* (O) day                            */
            long   *yy_o)                 /* (O) year                           */
{
    /*  Returns month, day, year from an integer in                */
    /*  the format YYY(Y)MMDD as in 840701 or 20011115 or 1120104  */

    *yy_o = date_i/10000;
    *mm_o = (date_i - *yy_o*10000)/100;
    *dd_o = date_i - *yy_o*10000 - *mm_o*100;
    return;
}

/*********************************************************************************
 * check whether year is leap year
 *
 ********************************************************************************/
int Isleap(long year)
{
    /*  Returns zero if 4 digit argument is not a leap year  */
    /*  otherwise returns 1.                                 */

    if (year%4 != 0)
        return (0); /* not divisible by 4 */
    if (year%100 != 0)
        return (1); /* divisible by 4 but not 100 */
    if (year%400 != 0)
        return (0); /* divisible by 100 but not 400 */
    return (1);         /* divisible by 400, so is a leap year */
}


/*********************************************************************************
 * return the number of days between date1_i and date2_i
 *
 ********************************************************************************/
long Daysact(long date1_i,                /* (I) start date                     */
             long date2_i)                /* (I) end date                       */
{
    long month1,day1,year1,month2,day2,year2,days,days1,days2;

    /* split dates into components */
    Dsplit(date1_i,&month1,&day1,&year1);
    Dsplit(date2_i,&month2,&day2,&year2);

    /* now convert */
    /* page 19 of SIA Blue Book */
    /* convert first date */
    days1 = day1 + cumdays[month1-1];
    days1 = days1 + (year1/400) - (year1/100) + (year1/4);

    if (month1 < 3)
        days1 = days1 - Isleap(year1);

    /* convert second date */
    days2 = day2 + cumdays[month2-1];
    days2 = days2 + (year2/400) - (year2/100) + (year2/4);

    if (month2 < 3)
        days2 = days2 - Isleap(year2);

    /* now compute difference */
    days = (year2 - year1)*365 + (days2 - days1);

    return (days);
}





/*****  GetDLOffset  ***************************************************/
/*
 *      Given a datelist and a targetDate,
 *      returns the offset of the element nearest to the targetdate
 *      if mode =
 *          CbkEXACT:  returns -999 if no matching targetdate is found
 *
 *          CbkLOWER:  returns the nearest offset that is LOWER or equal;
 *                     or -999 if all dates > targetdate
 *          CbkHIGHER: returns the nearest offset that is HIGHER or equal;
 *                     or -999 if all dates < targetdate
 */

int     GetDLOffset(int         NbDates,
                    long       *DL,
                    long        targetDate,
                    int         mode)
{
    int  offset = -999;
    int  i;

    if ((NbDates <=0) || (DL == NULL)) goto RETURN;

    /* find the largest i s.t. DL[i] <= targetDate
     * the i leaving the loop is always valid [0 .. NbDates-1]
     *
     * this final "largest" i may be:
     * ( >0): then DL[i] <= targetDate
     * (==0): then DL[k] >  targetDate for all 1<=k<=NbDates-1
     *             but we don't know about DL[0]
     */

    for (i=NbDates-1; i>0; i--)
    {
        if (DL[i] <= targetDate) break;
    }

    switch (mode)
    {
    case CbkEXACT:
        if (DL[i] == targetDate) offset = i;
        break;
    case CbkLOWER:
        if (DL[i] <= targetDate) offset = i;
        break;
    case CbkHIGHER:
        if (DL[i] >= targetDate) offset = i;
        else
        if (i<NbDates-1)   offset = i+1;
        break;
    default:
        goto RETURN;
    } /* switch (mode) */


RETURN:

    return (offset);

} /* GetDLOffset */
