/****************************************************************************/
/*      Do all the manipulations on dates: settlement, next business day,   */
/*      conversion, ...                                                     */
/*                                                                          */
/* IMPORTANT: This file contains conditionally compiled code depending upon */
/*            whether one is using IRX dates (i.e. -DESL_NEW_DATE).         */
/*                                                                          */
/*            With the exception of common include files and some common    */
/*            definitions, the code is organized so that the conditionally  */
/*            compiled code is placed near the top of the file, and common  */
/*            code (without any conditional compilation) is placed near the */
/*            end.  PLEASE respect this convention when adding/modifying    */
/*            this file.                                                    */
/*                                                                          */
/*            See esl_date.h for more information.                          */
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
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "esl_date.h"
#include "esl_alloc.h"
#include "esl_util.h"
#include "esl_error.h"
#include "irx/date.h"

#ifdef ESL_NEW_DATE
#   include "limits.h"
#   include "irx/dateutils.h"
#endif

#if defined(WIN32)
#define STRICMP stricmp
#elif defined(SOLARIS) || (defined(sun) && defined(unix))
/* Neither stricmp nor strcasecmp seem to be defined on Solaris. */
#include "irx/strutils.h" 
#define STRICMP irxStrcmpi
#else
#define STRICMP strcasecmp
#endif

#define HOLIDAY_FILE "HOLIDAY"
#define OKAY 0
#define ERR -1

static long JULIAN_EXCEL_OFFSET = 2415019L;

static long noleap[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
static long leap[13] =  {0,31,29,31,30,31,30,31,31,30,31,30,31};

static const char MonthNames[][4] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};


YMDDate YMDDateFromIRDate(IRDate date)
{
#ifdef ESL_NEW_DATE
    IrxTYearMonthDay ymd;

    irxDateToYMD(date, &ymd);

    return ymd.year * 10000 + ymd.month * 100 + ymd.day;
#else
    return date;
#endif
}


IRDate IRDateFromYMDDate(YMDDate date)
{
#ifdef ESL_NEW_DATE
    IRDate            iDate = INVALID_DATE;
    IrxTYearMonthDay ymd;
    
    if (date <= 0)
        return iDate;  

    ymd.year  = date / 10000;
    ymd.month = (date - ymd.year * 10000) / 100;
    ymd.day   = date - ymd.year * 10000 - ymd.month * 100;

    irxYMDToDate(&ymd, &iDate);

    return iDate;
#else
    return date;
#endif
}


IRDate IRDateFromTDate(long tDate)
{
#ifdef ESL_NEW_DATE
    return tDate;
#else
    IrxTYearMonthDay    ymd;

    if (irxDateToYMD(tDate, &ymd) == FAILURE)
    {
        DR_Error ("Could not convert TDate %ld to the internal date representation (YYYYMMDD).", tDate);
        return INVALID_DATE;
    }

    return ymd.year * 10000L + ymd.month * 100L + ymd.day;
#endif
}


/** Returns month, day, year from an integer in                
*   the format YYY(Y)MMDD as in 840701 or 20011115 or 1120104  */
void Dsplit(IRDate date_i, long *mm_o, long *dd_o, long *yy_o)
{
#ifdef ESL_NEW_DATE
    IrxTYearMonthDay ymd;
    
    irxDateToYMD(date_i, &ymd);

    *yy_o = ymd.year;
    *mm_o = ymd.month;
    *dd_o = ymd.day;
#else
    *yy_o = date_i/10000;
    *mm_o = (date_i - *yy_o*10000)/100;
    *dd_o = date_i - *yy_o*10000 - *mm_o*100;
#endif
}


/** Returns zero if 4 digit argument is not a leap year  *
 *  otherwise returns 1.                                 */
int Isleap(long year)
{
#ifdef ESL_NEW_DATE
    return IRX_IS_LEAP(year);
#else
    if (year%4 != 0)
        return 0; /* not divisible by 4 */
    if (year%100 != 0)
        return 1; /* divisible by 4 but not 100 */
    if (year%400 != 0)
        return 0; /* divisible by 100 but not 400 */
    return 1;         /* divisible by 400, so is a leap year */
#endif
}


IRDate   ThirdWed(long mth, long year)
{
#ifdef ESL_NEW_DATE
    return irxDateNthWeekDay(year, mth, IRX_WEDNESDAY, 3);
#else
    IRDate  date;

    date = Datepack(mth, 1L, year);

    while (Dayofwk(date) != 3)
        date++;

    /* Now we are at the 1st Wednesday */
    /* and must increment by two weeks */
    date = date + 14L;

    return date;
#endif
}


/** Packs mm_i,dd_i,yy_i into YYYYMMDD format  * 
 *  calls Y2toy4 on yy_i before packing        */
IRDate Datepack(long mm_i, long dd_i, long yy_i)
{
    IRDate date;

#ifdef ESL_NEW_DATE
    if (IRX_INVALID_DATE(date = irxDate(Y2toy4(yy_i), mm_i, dd_i)))
        return INVALID_DATE;
#else
    long y4;

    y4 = Y2toy4(yy_i);
    date = y4*10000 + mm_i*100 + dd_i;
#endif

    return date;
}


/**  Calculates day of week (0-6) of given date  */
long Dayofwk(IRDate date)
{
#ifdef ESL_NEW_DATE
    return IRX_DAY_OF_WEEK(date);
#else
    long mm,dd,yy;
    long *daysin;

    Dsplit(date,&mm,&dd,&yy);
    daysin = Isleap(yy) ? leap : noleap;
    while (mm > 1)
        dd+=daysin[--mm];
    if (yy > 0) {
        --yy;
        dd+=yy;
        dd+= yy/4 - yy/100 + yy/400; /* adjust for leap years */
    }
    return dd%7; /* 0 to 6*/
#endif
}


/*****************************************************************************/
/**                                                                          * 
 *  FUNCTION   Days360                                                       * 
 *                                                                           * 
 *     Counts days between two dates following the 30/360 covention as       * 
 *     specified in  the ISDA rule book. The method below agrees  with       * 
 *     the Analytics Library and with the routines implemented by  the       * 
 *     STIRT group (thanks to Neill Penney for useful discussions!)          * 
 *                                                                           * 
 *                                                                           * 
 *****************************************************************************/
long Days360(IRDate date1_i, IRDate date2_i)
{
    long days;

#ifdef ESL_NEW_DATE
    if (irxDayCountDays(date1_i, date2_i, IRX_B30_360, &days) != SUCCESS)
        days = LONG_MAX;
#else
    long month1,day1,year1,month2,day2,year2;

    /* Trouble shooting for 2/28 coupon */
    if (date1_i == date2_i)
        return 0L;

    /* split dates into components */
    Dsplit(date1_i,&month1,&day1,&year1);
    Dsplit(date2_i,&month2,&day2,&year2);

    if (day1 == 31)
        day1= 30;

    if ((day2 == 31) && (day1 == 30))
        day2 = 30;

    days = (year2-year1)*360 + (month2-month1)*30 + day2-day1;
#endif

    return days;
}


long Daysact(IRDate date1_i, IRDate date2_i)
{
    long days;

#ifdef ESL_NEW_DATE
    if (irxDayCountDays(date1_i, date2_i, IRX_ACT_365F, &days) != SUCCESS)
        days = LONG_MAX;
#else
    static long cumdays[12] = {0,31,59,90,120,151,181,212,243,273,304,334};

    long month1,day1,year1,month2,day2,year2,days1,days2;

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
#endif

    return days;
}


/** Returns a date in YYYYMMDD format based on moving        * 
 *  forward or backwards a number of calendar days.          */
IRDate Nxtday(IRDate date, long days)
{
#ifdef ESL_NEW_DATE
    IrxTDate          ret;
    IrxTDateInterval  interval;

    interval.prd     = days;
    interval.prd_typ = IRX_PRD_TYPE_DAY;
    interval.eom     = 0;

    irxDateAddInterval(date, interval, &ret);

    return ret;
#else
    long nxtm,nxtd,nxty,edays;
    long *daysin;

    Dsplit(date,&nxtm,&nxtd,&nxty);
    edays = nxtd + days ; /* positive or negative */
    daysin = Isleap(nxty) ? leap : noleap;
    if (days >= 0) { /* move forward */
        while (edays > daysin[nxtm]) {
            edays -= daysin[nxtm];
            nxtm++;
            /* see if into next year */
            if (nxtm > 12) {
                nxtm = 1;
                nxty++;
                daysin = Isleap(nxty) ? leap : noleap;
            }
        }
        /* nxtm, edays, nxty contain the proper date */
    }
    else { /* move backwards */
        while (edays < 1) {
            nxtm--;
            if (nxtm < 1) {
                nxty--;
                nxtm = 12;
                daysin = Isleap(nxty) ? leap : noleap;
            }
            edays +=daysin[nxtm];
        }
        /* nxtm, edays, nxty contain the proper date */
    }

    /* return result */
    return Datepack(nxtm,edays,nxty);
#endif
}


/** Returns a date based on moving forward or backwards (mths) from (date)
 *  The eom==0 is not supported for IRX dates as in the old code it was a
 *  'spill-over' flag. For eom==1 we get swap convention behavior 30+1M=30,
 *  if you need bond convention call irxDateAddMonths directly
 */
IRDate Nxtmth(IRDate date, long mths, long eom)
{
#ifdef ESL_NEW_DATE
    if (eom == 0)
    {
        DR_Error("Nxtmth - unsupported eom flag\n");
        return INVALID_DATE;
    }
    return irxDateAddMonths(date, mths, 0);
#else
    /** Returns a date in YYYYMMDD format based on moving        *
     *  forward or backwards (mths) from (date).                 *
     *  if eom > 0 the returned does not spill over to the next  *
     *  month. That is moving 1m from the 31rst of march gives   *
     *  the 30th of april. If eom=0 it does spill over: 1m from  *
     *  the 31rst of march is the 1rst of may.	             *
     *  a positive entry for mths => move forward                *
     *  a negative entry for mths => move backward.              */

    IRDate nxtdte;
    long *daysin,*daysout;
    long mm,dd,yy;
    long sign,nxtmm,nxtdd,nxtyy;

    Dsplit(date,&mm,&dd,&yy);
    if (mths > 0)
        sign = 1;
    else
        sign = -1;
    /* move forward by mths */
    nxtmm = mm + mths;
    nxtyy = yy;
    while (nxtmm > 12 || nxtmm < 1){
        nxtmm -= 12*sign;
        nxtyy+=sign;
    }
    /* now nxtyy and nxtmm have appropriate values */
        /* see if have to adjust days for end of month */

    daysin = Isleap(yy) ? leap : noleap;
    daysout = Isleap(nxtyy) ? leap : noleap;
    nxtdd = dd;

    if (eom) {
        if (nxtdd > daysout[nxtmm])
            nxtdd = daysout[nxtmm];
    }
    else {
        if (nxtdd > daysout[nxtmm]) {
            nxtdd = nxtdd - daysout[nxtmm];
            ++nxtmm;
            /*  Year will never be incremented because  */
            /*  December has 31 days.                   */
        }
    }

    nxtdte = nxtyy*10000 + nxtmm*100 + nxtdd;

    return nxtdte;
#endif
}


/**************************************************************************/
/**                                                                       * 
 * FUNCTION    DrDayCountFraction                                         * 
 *                                                                        * 
 * Calculates the daycount fraction between two dates according to one of * 
 * three methods: Actual/365Fixed, Actual/360 and 30/360.                 * 
 *                                                                        * 
 * RETURNS:  Success or Failure                                           * 
 *                                                                        */
/**************************************************************************/

int  DrDayCountFraction(IRDate Date1,   /**< (I) Start date             */
                        IRDate Date2,   /**< (I) End date               */
                        char     Conv,    /**< (I) Convention (A,0,3,5)   */
                        double  *frac)    /**< (O) Corresponding fraction */
{
#ifdef ESL_NEW_DATE
    int status = FAILURE;

    switch (Conv)
    {
        case 'A':
/*
            status = irxDayCountFraction(Date1, Date2, IRX_ACT_ACT, frac);
*/
            DR_Error("DrDayCountFraction: Convention Act/Act not supported!");
            return FAILURE;

        case '5':
            status = irxDayCountFraction(Date1, Date2, IRX_ACT_365F, frac);
            break;
        case '0':
            status = irxDayCountFraction(Date1, Date2, IRX_ACT_360, frac);
            break;  
        case '3':
            status = irxDayCountFraction(Date1, Date2, IRX_B30_360, frac);
            break;  
        default:
            DR_Error("DrDayCountFraction: "
                    "Convention (%c) for day count fraction not supported!", Conv);
            return FAILURE;
    }

    return status;
#else
    long    Ndays;

    switch (Conv)
    {
        case 'A':
            DR_Error("DrDayCountFraction: Convention Act/Act not supported!\n");
            return FAILURE;

        case '5':
            Ndays = Daysact(Date1,
                            Date2);
            *frac = (double)Ndays/365.0;
            break;
        case '0':
            Ndays = Daysact(Date1,
                            Date2);
            *frac = (double)Ndays/360.0;
            break;
        case '3':
            Ndays = Days360(Date1,
                            Date2);
            *frac = (double)Ndays/360.0;
            break;
        default:
            DR_Error("DrDayCountFraction: "
                    "Convention (%c) for day count fraction not supported!\n", Conv);
            return FAILURE;
    }

    return SUCCESS;
#endif
}


double  DrDcf(IRDate from, IRDate to, ESL_DCC dcc)
{
    double fraction = FLT_MAX;

#ifdef ESL_NEW_DATE
    switch (dcc)
    {
        case ESL_DCC_ACT_365:
            irxDayCountFraction(from, to, IRX_ACT_365F, &fraction);
        case ESL_DCC_ACT_ACT:
            irxDayCountFraction(from, to, IRX_ACT_ACT, &fraction);
            break;
        case ESL_DCC_ACT_360:
            irxDayCountFraction(from, to, IRX_ACT_360, &fraction);
            break;  
        case ESL_DCC_30_360:
            irxDayCountFraction(from, to, IRX_B30_360, &fraction);
            break;  
    }
#else
    assert(dcc == ESL_DCC_30_360 || dcc == ESL_DCC_ACT_360 || dcc == ESL_DCC_ACT_365);

    switch (dcc)
    {
        case ESL_DCC_ACT_365:
        case ESL_DCC_ACT_ACT: /* can not handle this - assert will capture it */
            fraction = (double)Daysact(from, to)/365.0;
            break;
        case ESL_DCC_ACT_360:
            fraction = (double)Daysact(from, to)/360.0;
            break;
        case ESL_DCC_30_360:
            fraction = (double)Days360(from, to)/360.0;
            break;
    }
#endif
    return fraction;
}


/****************************************************************************/
/*                                                                          */
/*   FUNCTION    DrDateFwdAny                                               */
/*                                                                          */
/*   Starting from a given date, this function calculates the date obtained */
/*   by incrementing a number of periods.                                   */
/*                                                                          */
int     DrDateFwdAny (  IRDate  StartDate, /* (I) Start date          */
                        int     NbPers,    /* (I) Number of periods   */
                        char    Freq,      /* (I) Frequency           */
                        char    FwdOrBwd,  /* (I) Forward or backward */
                        IRDate  *DateOut)  /* (O) Calculated date     */
{
    int     status = FAILURE;
    int     Factor;    
#ifdef ESL_NEW_DATE
    IrxTDateInterval  interval;
#else
    char    IntvalPerType;     /* Period type corresponding to frequency   */
    long    IntvalNbPers;      /* Nb of periods of type above corr to freq */
#endif

    if (FwdOrBwd == 'F')
        Factor = 1;
    else if (FwdOrBwd == 'B')
        Factor = -1;
    else
    {
        DR_Error ("DrDateFwdAny: Fourth argument must be F(fwd) or B(bwd)!");
        goto RETURN;
    }

#ifdef ESL_NEW_DATE
    switch (Freq)
    {
        case 'A':
            interval.prd_typ = IRX_PRD_TYPE_ANNUAL;
            break;
        case 'S': 
            interval.prd_typ = IRX_PRD_TYPE_SEMI_ANNUAL;
            break;
        case 'Q': 
            interval.prd_typ = IRX_PRD_TYPE_QUARTER;
            break;
        case 'I': /* NbPers IMM date */
            break;
        case 'M':
            interval.prd_typ = IRX_PRD_TYPE_MONTH;
            break;
        case 'W':
            interval.prd_typ = IRX_PRD_TYPE_WEEK;
            break;
        case 'D':
            interval.prd_typ = IRX_PRD_TYPE_DAY;
            break;
        default:
            DR_Error ("DrDateFwdAny: Unrecognised frequency: %c!", Freq);
            goto RETURN;
    }  /* switch */

    if (Freq == 'I')
         *DateOut =  Nxtimm(StartDate, (long)(Factor*NbPers));
    else
    {
        interval.prd = Factor * NbPers;
        interval.eom = 0;
        irxDateAddInterval(StartDate, interval, DateOut);
    }
#else
    switch (Freq)
    {
        case 'A': IntvalNbPers   = 12 * NbPers;
                  IntvalPerType  = 'M';
                  break;
        case 'S': IntvalNbPers   = 6 * NbPers;
                  IntvalPerType  = 'M';
                  break;
        case 'Q': IntvalNbPers   = 3 * NbPers;
                  IntvalPerType  = 'M';
                  break;
        case 'I': IntvalNbPers   = NbPers;
                  IntvalPerType  = 'I';
                  break;
        case 'M': IntvalNbPers   = NbPers;
                  IntvalPerType  = 'M';
                  break;
        case 'W': IntvalNbPers   = 7 * NbPers;
                  IntvalPerType  = 'D';
                  break;
        case 'D': IntvalNbPers   = NbPers;
                  IntvalPerType  = 'D';
                  break;
        default:
                  DR_Error ("DrDateFwdAny: Unrecognised frequency: %c!", Freq);
                  goto RETURN;
    }  /* switch */

    if (IntvalPerType == 'M')
        *DateOut =  Nxtmth(StartDate, (long)(Factor*IntvalNbPers), 1L);
    else if (IntvalPerType == 'I')
        *DateOut =  Nxtimm(StartDate, (long)(Factor*IntvalNbPers));
    else
        *DateOut =  Nxtday(StartDate, (long)(Factor*IntvalNbPers));
#endif

    status = SUCCESS;

  RETURN:
    if (status == FAILURE)
        DR_Error("DrDateFwdAny: Failed.");

    return status;
}


/** Creates an IRDate from a date string in MM/DD/YY format.     */ 
IRDate eval_date(char *datest)
{
    long year, month, day;

    year = (datest[6] - '0') * 10 + (datest[7] - '0');
    year = (year > 73) ? 1900 + year : 2000 + year;
    month = (datest[0] - '0') * 10 + (datest[1] - '0');
    day = (datest[3] - '0') * 10 + (datest[4] - '0');

#ifdef ESL_NEW_DATE
    return irxDate(year, month, day);
#else
    return year * 10000 + month * 100 + day;
#endif
}


/**  Create an IRDate from a date string in "DD-MON-YYYY" (e.g. "05-DEC-2006") format */
IRDate eval_date2(char *datest)
{
    long year, month, day;
    char string[4];

    year  = (datest[7] - '0') * 1000 + (datest[8] - '0') * 100 
            + (datest[9] - '0') * 10 + (datest[10] - '0');
    day = (datest[0] - '0') * 10 + (datest[1] - '0');
        
    string[0] = datest[3];
    string[1] = datest[4];
    string[2] = datest[5];
    string[3] = '\0';

    for (month = 0; month < 12; ++month)
        if (STRICMP(string, MonthNames[month]) == 0)
            break;
    ++month;
                
#ifdef ESL_NEW_DATE
    return irxDate(year, month, day);
#else
    if (month > 12)
        month = 0;
    return year * 10000 + month * 100 + day;
#endif
}


/************************************************************************************
**
** IMPORTANT:  ALL CODE BELOW THIS POINT IS SUPPOSED TO WORK WITH BOTH OLD (YYYYMMDD)
**             DATES (YMDDates) AND NEW IRX DATES.
**
************************************************************************************/


/**  Returns TRUE if the date is an IMM date    */
int Isimm(IRDate    date)
{
        int   status = FALSE;

        long  day;
        long  mth;
        long  year;

        Dsplit(date, &mth, &day, &year);

        if ((mth==3) || (mth==6) || (mth==9) || (mth==12))
        {
            if (date == ThirdWed(mth, year))
            {
                status = TRUE;
            }
        }

        return (status);
}


/** This routine checks whether the date passed in is a holiday  *
 *  or not according to the HOLIDAYS file. It returns (0) if     *
 *  it is not and (-1) if it is.  Dates in the holiday file are  *
 *  expected to be in YYYYMMDD format.                           */  
int IsHoliday(IRDate date_i)
{
    FILE *fp1;
    int xflag = ERR;
    int rcode = ERR;
    IRDate holiday;
    char line1[82];

    if ((fp1=fopen(HOLIDAY_FILE,"r")) == NULL)
    {
        fprintf(stdout,"Can't Open file %s\n",HOLIDAY_FILE);
        return(rcode);
    }
    while (xflag != OKAY)
    {
        if (fgets(line1,80,fp1) == NULL)
	{
            xflag = OKAY;
            rcode = OKAY;
        }
        if (xflag != OKAY)
	{
            line1[8]='\0';
            holiday = IRDateFromYMDDate(atol(line1));
            if (holiday > date_i)
	    {
                    rcode=OKAY;
                    xflag=OKAY;
            }
            else if (holiday == date_i)
	    {
                rcode = ERR;;
                xflag = OKAY;
            }
        }
    }
    fclose(fp1);
    return(rcode);
}

int Date_CheckAndReport(
		IRDate	date) /**< (I) Date */
{
	int rc;     /* Result code */

	/* Check it */
	rc = Dateok(date);

	/* Report it */
	switch (rc)
	{
	case 0:		/* OK */
		break;
	case 1:		/* Bad year */
		DR_Error((char*) "Date_CheckAndReport: Date %ld has a bad year.\n", date);
				break;
	case 2:		/* Bad month */
		DR_Error((char*) "Date_CheckAndReport: Date %ld has a bad month.\n", date);
				break;
	case 3:		/* Bad day */
		DR_Error((char*) "Date_CheckAndReport: Date %ld has a bad day.\n", date);
				break;
	default:	/* Unknown */
		DR_Error((char*) "Date_CheckAndReport: %ld is a bad date.\n", date);
				break;
	}

	/* Done */
	return rc;
}


/** Converts a 2 digit year (84) to a 4 digit year (1984);       *
 *  assumes that 2 digit years less than 50 are after year 2000  */
long Y2toy4(long year_i)
{
    long y4_o,y2;

    y2 = year_i;
    if (y2 >= 100)
        y4_o = y2;
    else {
    if (y2 <= 50)
        y2 = y2 + 100;
        y4_o = y2 + 1900;
    }
    return (y4_o);
}


/** Converts a 4 digit year to a 2 digit year;  *
 *  years after 99 are 00, 01, etc.             */
long Y4toy2(long year_i)
{
    long y2_o;

    y2_o = year_i%100;
    return (y2_o);
}

/**  Convert input date to "MM/DD/YY"  */
void Y2date_str(IRDate date, char *string)
{
    long mm,dd,yy;

    Dsplit(date,&mm,&dd,&yy);
    yy = Y4toy2(yy);
    sprintf(string,"%02ld/%02ld/%02ld",mm,dd,yy);

    return;
}


long Months360(IRDate date1_i, IRDate date2_i)
{
    long month1,day1,year1,month2,day2,year2,months;

    /* Trouble shooting for 2/28 coupon */
    if (date1_i == date2_i)
    {
        return(0L);
    }

    /* split dates into components */
    Dsplit(date1_i,&month1,&day1,&year1);
    Dsplit(date2_i,&month2,&day2,&year2);

    months = (year2-year1)*12 + (month2-month1);

    return (months);
}


/** a positive entry for nbPeriods => move forward     *
 *  a negative entry for nbPeriods => move backward.   */
IRDate  Nxtimm(IRDate  date, long nbPeriods)
{
    IRDate   immdate;
    long      day,    mth,    year;
    long           immmth, immyear;

    int    Fwd;

    int    i;


    if (nbPeriods >= 0)
    {
        Fwd = TRUE;
    }
    else
    {
        nbPeriods = -nbPeriods;
        Fwd = FALSE;
    }


    immdate = date;

    for (i=0; i<nbPeriods; i++)
    {
        Dsplit(date, &mth, &day, &year);

        immmth  = (((mth-1)/3) + 1) * 3;
        immyear = year;

        immdate = ThirdWed(immmth, immyear);

        if (Fwd)
        {
            if (immdate <= date)
            {
                immmth = immmth + 3;
                if (immmth == 15)
                {
                    immmth = 3;
                    immyear ++;
                }

                immdate = ThirdWed(immmth, immyear);
            }
        }
        else
        {
            if (immdate >= date)
            {
                immmth = immmth - 3;
                if (immmth == 0)
                {
                    immmth = 12;
                    immyear --;
                }

                immdate = ThirdWed(immmth, immyear);
            }
        }


        date = immdate;

    } /* For i */


    return(immdate);

}


/** This routine returns a new date after advancing or going  *
 *  backward a number of weekdays (Monday - Friday).          */
IRDate Nxtwkday(IRDate date, long advance)
{
    IRDate work_date;
    long dow;

    work_date = date;
    if (advance == 0)
        return(date);
    else if (advance < 0) {  /* from date, go backward */
        while (advance != 0) {
            dow = Dayofwk(work_date);
            advance++;
            if (dow == 1)  /* Monday */
                work_date = Nxtday(work_date,(long)-3);
            else if (dow == 0)  /* Sunday */
                work_date = Nxtday(work_date,(long)-2);
            else  /* Tuesday - Saturday */
                work_date = Nxtday(work_date,(long)-1);
        }
    }
    else {  /* from date, go forward */
        while (advance != 0) {
            dow = Dayofwk(work_date);
            advance--;
            if (dow == 5)  /* Friday */
                work_date = Nxtday(work_date,(long)3);
            else if (dow == 6)  /* Saturday */
                work_date = Nxtday(work_date,(long)2);
            else  /* Sunday - Thursday */
                work_date = Nxtday(work_date,(long)1);
        }
    }
    return(work_date);
}

/** This routine returns a new date after advancing or going  *
 *  backward a number of business days (Monday - Friday &     *
 *  holidays). It checks for holidays in HOLIDAY file.        *
 *  The function returns FAILURE if it is not possible to     *
 *  call the IsHoliday function successfully.                 */
IRDate Nxtbusday(IRDate date, long advance)
{
    IRDate work_date;

    work_date = date;

    if (advance == 0)
        return(date);
    else if (advance < 0) {  /* from date, go backward */
        while (advance != 0) {
            advance++;
            work_date = Nxtwkday(work_date,(long)-1);
            while(IsHoliday(work_date))
                work_date=Nxtwkday(work_date,(long)-1);
        }
    }
    else {  /* from date, go forward */
        while (advance != 0) {
            advance--;
            work_date = Nxtwkday(work_date,(long)1);
            while(IsHoliday(work_date))
                work_date=Nxtwkday(work_date,(long)1);
        }
    }
    return(work_date);
}


/*****************************************************************************/
/**                                                                          *
 *   FUNCTION    DrNewEventListFromFreqWithInterpType                        *
 *                                                                           *
 * Better interface to DrNewEventListFromFreq where an interpolation         *
 * method should be specified (instead of clients negating NbInpDates).      *
 *                                                                           */
 EVENT_LIST * DrNewEventListFromFreqWithInterpType(
	ESL_INTERP_TYPE interpType,  /**< (I) Interpolation type to use*/
        int         NbInpDates,  /**< (I) Nb of dates input directly   */
        IRDate       *InpDates,    /**< (I) Dates given directly by user */
        char        Freq,        /**< (I) Frequency of event           */
        char        Stub,        /**< (I) Stub location Front or Back  */
        char        DatesIn,     /**< (I) Y=input dates must be in list*/
        double     *Curve0,      /**< (I) Set of values for event      */
        double     *Curve1,      /**< (I) Set of values for event      */
        double     *Curve2,      /**< (I) Set of values for event      */
        double     *Curve3,      /**< (I) Set of values for event      */
        double     *Curve4)      /**< (I) Set of values for event      */
{
    if (interpType == LINEAR_INTERP || interpType == STAIRCASE_INTERP)
      return DrNewEventListFromFreq(
                (interpType == STAIRCASE_INTERP ? -1 : 1) * NbInpDates,
                InpDates,
                Freq,
                Stub,
                DatesIn,
                Curve0,
                Curve1,
                Curve2,
                Curve3,
                Curve4);
    else
    {
        DR_Error("Interpolation type must be either linear or staircase!");
        return NULL;
    }
}


/*****************************************************************************/
/**                                                                          *
 *   FUNCTION    DrNewEventListFromFreq                                      *
 *                                                                           *
 *   Added 8/04 -- allow staircase intepolation (use NEGATIVE NbInpDates)    *
 *                                                                           */
 EVENT_LIST * DrNewEventListFromFreq
               (int    NbInpDates, /**< (I) Nb of dates input directly   */
                IRDate   *InpDates,  /**< (I) Dates given directly by user */
                char    Freq,      /**< (I) Frequency of event           */
                char    Stub,      /**< (I) Stub location Front or Back  */
                char    DatesIn,   /**< (I) Y=input dates must be in list*/
                double *Curve0,    /**< (I) Set of values for event      */
                double *Curve1,    /**< (I) Set of values for event      */
                double *Curve2,    /**< (I) Set of values for event      */
                double *Curve3,    /**< (I) Set of values for event      */
                double *Curve4)    /**< (I) Set of values for event      */
{
    int status = FAILURE; /* Until proven SUCCESSFUL */

    EVENT_LIST *EventListLocal = NULL;

    int         NbCurves;
    int         CurrentIsMatched;

    int         i;       /* Convenience variables */
    int         AuxIdx;
    int         useStaircase = FALSE;


    /* 8/04: enable staircase interpolation without changing prototype
     * and affecting current products */
    if (NbInpDates < 0)
    {
        NbInpDates   = - NbInpDates;
        useStaircase = TRUE;
    }


    /* Evaluate number of user curves */
    if (Curve0 == NULL)
    {
        NbCurves = 0;
    }
    else
    {
        if (Curve1 == NULL)
        {
            NbCurves = 1;
        }
        else
        {
            if (Curve2 == NULL)
            {
                NbCurves = 2;
            }
            else
            {
                if (Curve3 == NULL)
                {
                    NbCurves = 3;
                }
                else
                {
                    if (Curve4 == NULL)
                    {
                        NbCurves = 4;
                    }
                    else
                    {
                        NbCurves = 5;
                    }
                }
            }
        }
    }


    /* Allocate and initialise structure to be returned */
    EventListLocal = (EVENT_LIST *)malloc((size_t)sizeof(EVENT_LIST));
    if (EventListLocal == NULL)
    {
        DR_Error ("DrNewEventListFromFreq: could not allocate memory "
                    "for event list!");
        goto FREE_MEM_AND_RETURN;
    }
    EventListLocal->Dates    = NULL;
    for (i=0; i<MAXNBEVCURVES; i++)
    {
        EventListLocal->Curve[i] = NULL;
    }
    EventListLocal->NbCurves = NbCurves;



    /* Three steps required:                                            */
    /*     1 - Generate dates as per frequency                          */
    /*     2 - Check that user dates are a subset of the above date set */
    /*     3 - Interpolate user supplied curves                         */


    /* STEP 1 - Date generation */
    if ((Freq == 'I') || /* All these  three cases */
        (Freq == 'N') || /* imply only input dates */
        (Freq == 'E'))   /* should appear in list  */
    {
        EventListLocal->NbEntries = NbInpDates;
        EventListLocal->Dates = (IRDate *) DR_Array (IDATE, 0, NbInpDates-1);
        if (EventListLocal->Dates == NULL)
        {
            DR_Error ("DrNewEventListFromFreq: unable to allocate memory for "
                        "dates in event list!");
            goto FREE_MEM_AND_RETURN;
        }
        for (i=0; i<NbInpDates; i++)
        {
            EventListLocal->Dates[i] = InpDates[i];
        }
    }
    else
    {
        if (DateListFromFreq(InpDates[0],
                             InpDates[NbInpDates-1],
                             Freq,
                             Stub,
                             &(EventListLocal->NbEntries), /* Output */
                             &(EventListLocal->Dates))     /* Output */
                                             == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }


    /* STEP 2 - Check that user dates are subset (1st Date is OK) */
    if ((DatesIn == 'Y') && (Freq != 'I') && (Freq != 'E'))
    {
        for (i=1, AuxIdx=1; i<NbInpDates; i++)
        {
            CurrentIsMatched = FALSE;
            while (AuxIdx < EventListLocal->NbEntries)
            {
                if (InpDates[i] == EventListLocal->Dates[AuxIdx])
                {
                    CurrentIsMatched = TRUE;
                    break;
                }
                AuxIdx++;
            }
            if (!CurrentIsMatched)
            {
                DR_Error ("DrNewEventListFromFreq: could not match input date "
                            "with given frequency!");
                goto FREE_MEM_AND_RETURN;
            }
        }
    }


    /* STEP 3 - Allocate memory  and interpolate  */
    /*          high and low barriers and rebate  */
    for (i=0; i<NbCurves; i++)
    {
        EventListLocal->Curve[i] =
                (double *) DR_Array (DOUBLE, 0, EventListLocal->NbEntries-1);
        if (EventListLocal->Curve[i] == NULL)
        {
            DR_Error ("DrNewEventListFromFreq: unable to allocate memory for "
                        "curves in event list!");
            goto FREE_MEM_AND_RETURN;
        }
    }
    /* Interpolation is only necessary if freq != I,E */
    if ((Freq == 'I') || (Freq == 'E'))
    {
        IRDate  CurrDate;
        for (i=0; i<EventListLocal->NbEntries; i++)
        {
            CurrDate = EventListLocal->Dates[i];
            switch(NbCurves)
            {
                case 5:
                    EventListLocal->Curve[4][i] = Curve4[i];
                case 4:
                    EventListLocal->Curve[3][i] = Curve3[i];
                case 3:
                    EventListLocal->Curve[2][i] = Curve2[i];
                case 2:
                    EventListLocal->Curve[1][i] = Curve1[i];
                case 1:
                    EventListLocal->Curve[0][i] = Curve0[i];
            }
        }
    }
    else
    {
        /* When  frequency is  != from I/E,  then a  minimum */
        /* of two dates will have been provided  by the user */
        int   IdxLeft  = 0;
        int   IdxRight = 1;
        IRDate  CurrDate;
        for (i=0; i<EventListLocal->NbEntries; i++)
        {
            CurrDate = EventListLocal->Dates[i];

            /* decide intepolation type: staircase or linear */
            if (useStaircase)
            {
                if (CurrDate >= InpDates[IdxRight])
                {
                    IdxLeft++;
                    IdxRight++;
                }

                switch(NbCurves)
                {
                    case 5:
                        EventListLocal->Curve[4][i] = Curve4[IdxLeft];
                    case 4:
                        EventListLocal->Curve[3][i] = Curve3[IdxLeft];
                    case 3:
                        EventListLocal->Curve[2][i] = Curve2[IdxLeft];
                    case 2:
                        EventListLocal->Curve[1][i] = Curve1[IdxLeft];
                    case 1:
                        EventListLocal->Curve[0][i] = Curve0[IdxLeft];
                }
            }
            else /* linear */
            {
                if (CurrDate > InpDates[IdxRight])
                {
                    IdxLeft++;
                    IdxRight++;
                }

                switch(NbCurves)
                {
                    case 5:
                        dlinterp (CurrDate,
                                  &(EventListLocal->Curve[4][i]),
                                  InpDates[IdxLeft],
                                  InpDates[IdxRight],
                                  Curve4[IdxLeft],
                                  Curve4[IdxRight]);
                    case 4:
                        dlinterp (CurrDate,
                                  &(EventListLocal->Curve[3][i]),
                                  InpDates[IdxLeft],
                                  InpDates[IdxRight],
                                  Curve3[IdxLeft],
                                  Curve3[IdxRight]);
                    case 3:
                        dlinterp (CurrDate,
                                  &(EventListLocal->Curve[2][i]),
                                  InpDates[IdxLeft],
                                  InpDates[IdxRight],
                                  Curve2[IdxLeft],
                                  Curve2[IdxRight]);
                    case 2:
                        dlinterp (CurrDate,
                                  &(EventListLocal->Curve[1][i]),
                                  InpDates[IdxLeft],
                                  InpDates[IdxRight],
                                  Curve1[IdxLeft],
                                  Curve1[IdxRight]);
                    case 1:
                        dlinterp (CurrDate,
                                  &(EventListLocal->Curve[0][i]),
                                  InpDates[IdxLeft],
                                  InpDates[IdxRight],
                                  Curve0[IdxLeft],
                                  Curve0[IdxRight]);
                }
            }
        }
    }


    status = SUCCESS;

    FREE_MEM_AND_RETURN:


    if (status == FAILURE)
    {
        if (EventListLocal != NULL)
        {
            DrFreeEventList(EventListLocal);
        }
        DR_Error("DrNewEventListFromFreq: Failed.");
        return(NULL);
    }
    else
    {
        return(EventListLocal);
    }

}   /* End of DrNewEventListFromFreq() */


/****************************************************************************/
/*                                                                          *
 *   FUNCTION    DateListFromFreq                                           *
 *                                                                          *
 *   Generates a date list including the start and end dates given and in   *
 *   accordance with the frequency and stub convention passed in.           *
 *                                                                          */
 int DateListFromFreq
              (IRDate       StartDate,     /**< (I) Start of date list       */
               IRDate       EndDate,       /**< (I) End of date list         */
               char       Freq,          /**< (I) Frequency of list        */
               char       StubConv,      /**< (I) Stub location Front/back */
               int       *NbDates,       /**< (O) Number of dates in list  */
               IRDate     **DateList)      /**< (O) List of dates asc order  */
{
    int status = FAILURE;

    char       IntvalPerType;   /* Period type corresponding to frequency   */
    long       IntvalNbPers;    /* Nb of periods of type above corr to freq */

    IRDate       AuxDate;                 /* Variables for local use          */
    int        Factor;
    int        IndexStart;

    IRDate      *DateListLocal = NULL;    /* Local copy of output list        */
    int        NbDatesLocal=0;          /* Local copy of output value       */

    int        StubsNotAllowed;  /* If caller requires that start and end   */
                                 /* dates define an exact number of pers    */


    int        i;



    /* Simple input checks*/
    if (StartDate > EndDate)
    {
        DR_Error ("DateListFromFreq: Start date must be before end date!");
       	goto RETURN;
    }
    if ((StubConv != 'F') && (StubConv != 'B') && (StubConv != 'N'))
    {
        DR_Error ("DateListFromFreq: Stub convention must be F,B or N!");
       	goto RETURN;
    }
    if (StubConv == 'N')
    {
        StubsNotAllowed = TRUE;
        StubConv = 'B';
    }
    else
    {
        StubsNotAllowed = FALSE;
    }


    switch(Freq)
    {
        case 'A': IntvalNbPers   = 12L;
                  IntvalPerType  = 'M';
                  break;
        case 'S': IntvalNbPers   = 6L;
                  IntvalPerType  = 'M';
                  break;
        case 'Q': IntvalNbPers   = 3L;
                  IntvalPerType  = 'M';
                  break;
        case 'I': IntvalNbPers   = 1L;
                  IntvalPerType  = 'I';
                  break;
        case 'M': IntvalNbPers   = 1L;
                  IntvalPerType  = 'M';
                  break;
        case 'W': IntvalNbPers   = 7L;
                  IntvalPerType  = 'D';
                  break;
        case 'D': IntvalNbPers   = 1L;
                  IntvalPerType  = 'D';
                  break;
        default:
                  DR_Error ("DateListFromFreq: Unrecognised frequency (%c) !", Freq);
       	       	  goto RETURN;
    }


    /* Count dates */
    if (StubConv == 'B')
    {
        i = 0;
        do
        {
            i++;
            if (IntvalPerType == 'M')
            {
                AuxDate = Nxtmth(StartDate,
                                 (long)(i*IntvalNbPers),
                                 1L);
            }
            else if (IntvalPerType == 'I')
            {
                AuxDate = Nxtimm(StartDate,
                                 (long)(i*IntvalNbPers));
            }
            else
            {
                AuxDate = Nxtday(StartDate,
                                 (long)(i*IntvalNbPers));
            }


        } while (AuxDate < EndDate);
        NbDatesLocal = i+1;

        if (StubsNotAllowed)
        {
            if (AuxDate != EndDate)
            {
                DR_Error("DateListFromFreq: Start and end date must define "
                    "an integer number of periods, if stubs are not allowed!\n");
                goto RETURN;
            }
        }

    }
    else
    {
        i = 0;
        do
        {
            i++;
            if (IntvalPerType == 'M')
            {
                AuxDate = Nxtmth(EndDate,
                                 (long)(-i*IntvalNbPers),
                                 1L);
            }
            else if (IntvalPerType == 'I')
            {
                AuxDate = Nxtimm(EndDate,
                                 (long)(-i*IntvalNbPers));
            }
            else
            {
                AuxDate = Nxtday(EndDate,
                                 (long)(-i*IntvalNbPers));
            }


        } while (AuxDate > StartDate);
        NbDatesLocal = i+1;
    }


    /* Allocate memory for date list */
    DateListLocal = (IRDate *) DR_Array (IDATE, 0, NbDatesLocal-1);
    if (DateListLocal == NULL)
    {
        DR_Error ("DateListFromFreq: Unable to allocate memory for datelist!");
       	goto RETURN;
    }


    /* Finally generate list */
    DateListLocal[0]              = StartDate;
    DateListLocal[NbDatesLocal-1] = EndDate;
    if (StubConv == 'B')
    {
        AuxDate = StartDate;
        Factor  = 1;
        IndexStart = 0;
    }
    else
    {
        AuxDate = EndDate;
        Factor = -1;
        IndexStart = NbDatesLocal-1;
    }
    for (i=1; i<NbDatesLocal-1; i++)
    {
        if (IntvalPerType == 'M')
        {
            DateListLocal[IndexStart+i*Factor] =
                            Nxtmth(AuxDate,
                                   (long)(Factor*i*IntvalNbPers),
                                   1L);
        }
        else if (IntvalPerType == 'I')
        {
            DateListLocal[IndexStart+i*Factor] =
                            Nxtimm(AuxDate,
                                   (long)(Factor*i*IntvalNbPers));
        }
        else
        {
            DateListLocal[IndexStart+i*Factor] =
                            Nxtday(AuxDate,
                                   (long)(Factor*i*IntvalNbPers));
        }
    }



    /* All is well if we got here, so pass values back to caller */
    *NbDates  = NbDatesLocal;
    *DateList = DateListLocal;

    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (DateListLocal, IDATE, 0, NbDatesLocal-1);
        DR_Error("DateListFromFreq: Failed.");
    }

    return(status);

}   /* End of DateListFromFreq() */


/****************************************************************************/
/**                                                                         *
 *   FUNCTION    DrDatesInSchedule                                          *
 *                                                                          *
 *   Returns SUCCESS if all dates passed in are contained in the schedule   *
 *   specified by a start date, an end date, a frequency and a stub conv.   *
 *                                                                          *
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       *
 *                                                                          */
 int DrDatesInSchedule
              (int        NbDatesIn,     /**< (I) Number of dates to check   */
               IRDate      *DatesIn,       /**< (I) User dates                 */
               IRDate       StartDate,     /**< (I) Start of date list         */
               IRDate       EndDate,       /**< (I) End of date list           */
               char       Freq,          /**< (I) Frequency of list          */
               char       StubConv)      /**< (I)                            */
{
    int status = FAILURE;

    int        NbDatesLocal;
    IRDate      *DateListLocal = NULL;    /* For local use only */

    int        i, j;

    /* Pre-checks */
    if (NbDatesIn == 0)
        return(SUCCESS); /* Nothing to check */

    if (DatesIn[0] < StartDate)
    {
        DR_Error("DrDatesInSchedule: first supplied date precedes start "
                    "of date list!\n");
        return(FAILURE);
    }
    if (DatesIn[NbDatesIn-1] > EndDate)
    {
        DR_Error("DrDatesInSchedule: last supplied date is beyond end "
                    "of date list!\n");
        return(FAILURE);
    }

    /* Build date list to compare against */
    if (DateListFromFreq(StartDate,
                         EndDate,
                         Freq,
                         StubConv,
                         &NbDatesLocal,
                         &DateListLocal) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Perform the comparison */
    i=0;
    for (j=0; j<NbDatesLocal; j++)
    {
        if (DateListLocal[j] == DatesIn[i])
            i++;

        if (i>NbDatesIn-1)
            break;
    }

    if (i < NbDatesIn)
    {
        DR_Error("DrDatesInSchedule: could not match all dates!\n");
        goto FREE_MEM_AND_RETURN;
    }


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (DateListLocal, IDATE, 0, NbDatesLocal-1);

    return(status);

}   /* End of DrDatesInSchedule() */



/****************************************************************************/
/**                                                                         *
 *   FUNCTION    DrSameDateSchedules                                        *
 *                                                                          *
 *   Returns SUCCESS if dates passed in are the same as dates generated by  *
 *   specifing a start date, an end date, a frequency, a stub conv,         *
 *   arrears reset and a value date, which is used to cut the list before it*
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       *
 *                                                                          */
 int DrSameDateSchedules
              (int       *NbDatesIn,     /**< (I/0) Number of dates to check */
               IRDate      *DatesIn,       /**< (I) User dates                 */
               IRDate       StartDate,     /**< (I) Start of date list         */
               IRDate       EndDate,       /**< (I) End of date list           */
               IRDate       ValueDate,     /**< (I) Value date                 */
               char       Freq,          /**< (I) Frequency of list          */
               char       StubConv,      /**< (I)                            */
               char       Arrears)       /**< (I) 'Y' if reset-in-arrears    */
{
    int status = FAILURE;

    int        NbDatesLocal;
    int        NbResetsLocal;
    IRDate      *DateListLocal = NULL;
    IRDate      *ResetListLocal = NULL;
    int        i;

    /* Build date list to compare against */
    if (DateListFromFreq(StartDate,
                         EndDate,
                         Freq,
                         StubConv,
                         &NbDatesLocal,
                         &DateListLocal) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Create reset date list */
    if(Arrears == 'Y')
    {
        ResetListLocal = DateListLocal + 1;
    }
    else
    {
        ResetListLocal = DateListLocal;
    }
    NbResetsLocal = NbDatesLocal - 1;

    /* Check if given dates are reset dates */
    if (DrDatesInSet(*NbDatesIn,
                     DatesIn,
                     NbResetsLocal,
                     ResetListLocal) == FAILURE)
    {
        DR_Error("DrSameDateSchedules: Input dates are not reset dates!");
        goto FREE_MEM_AND_RETURN;
    }

    /* Perform the comparison */
    i = 0;
    while((i < NbResetsLocal) && (ResetListLocal[i] < ValueDate))
    {
        if (DatesIn[i] != ResetListLocal[i])
        {
            DR_Error ("DrSameDateSchedules: Fixing date %ld is missing!", ResetListLocal[i]);
            goto FREE_MEM_AND_RETURN;
        }
        i++;
    }

    *NbDatesIn = i;

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (DateListLocal, IDATE, 0, NbDatesLocal-1);

    return(status);

}   /* End of DrSameDateSchedules() */


/****************************************************************************/
/**                                                                         *
 *   FUNCTION    DrSameDateSets                                             *
 *                                                                          *
 *   Returns SUCCESS if dates passed in are the same as dates generated by  *
 *   cuting a givne date list up to a value date.                           *
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       *
 *                                                                          */
int DrSameDateSets
              (int       *NbDates1,      /**< (I/0) Number of dates to check */
               IRDate      *Dates1,        /**< (I) User dates                 */
               int        NbDates2,      /**< (I) Number of dates to cmp     */
               IRDate      *Dates2,        /**< (I) Cpm dates                  */
               IRDate       ValueDate)     /**< (I) Value date                 */
{
    int     i;
    int     status = FAILURE;

    /* Check if dates1 are a subset of dates2 */
    if (DrDatesInSet(*NbDates1,
                     Dates1,
                     NbDates2,
                     Dates2) == FAILURE)
    {
        DR_Error("DrSameDateSets: Input dates are not a subset!");
        goto FREE_MEM_AND_RETURN;
    }

    /* Perform the comparison */
    for (i = 0; i < NbDates2; i++)
    {
        if (Dates2[i] < ValueDate)
        {
            if (Dates2[i] != Dates1[i])
            {
                DR_Error ("DrSameDateSets: Fixing date %ld is missing!", Dates2[i]);
                goto FREE_MEM_AND_RETURN;
            }
        }
        else break;
    }

    *NbDates1 = i;

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    return(status);

}   /* End of DrSameDateSets() */


/****************************************************************************/
/**                                                                         *
 *   FUNCTION    DrSameDateSetsFlows                                        *
 *                                                                          *
 *   Returns SUCCESS if dates passed in are the same as dates generated by  *
 *   cuting a givne date list up to a value date.                           *
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       *
 *                                                                          */
int DrSameDateSetsFlows
              (int       *NbDates1,      /**< (I/0) Number of dates to check */
               IRDate      *Dates1,        /**< (I) User dates                 */
               int        NbDates2,      /**< (I) Number of dates to cmp     */
               IRDate      *Dates2,        /**< (I) Cpm dates                  */
               IRDate       ValueDate)     /**< (I) Value date                 */
{
    int     i;
    int     status = FAILURE;

    /* Check if dates1 are a subset of dates2 */
    if (DrDatesInSet(*NbDates1,
                     Dates1,
                     NbDates2,
                     Dates2) == FAILURE)
    {
        DR_Error("DrSameDateSets: Input dates are not a subset!");
        goto FREE_MEM_AND_RETURN;
    }

    /* Perform the comparison */
    for (i = 0; i < NbDates2; i++)
    {
        if (Dates2[i] <= ValueDate)
        {
            if (Dates2[i] != Dates1[i])
            {
                DR_Error ("DrSameDateSets: Fixing date %ld is missing!", Dates2[i]);
                goto FREE_MEM_AND_RETURN;
            }
        }
        else break;
    }

    *NbDates1 = i;

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    return(status);

}   /* End of DrSameDateSetsFlows() */



/****************************************************************************/
/**                                                                         *
 *   FUNCTION    DrDatesIn2FreqSchedule                                     *
 *                                                                          *
 *   Compares the passed-in date list to a date list generated in the       *
 *   current (with respect to the value date) pmt period according to some  *
 *   obs frequency. The current pmt period is found from the start date,    *
 *   end date, pmt frequency, stub convention and a value date              *
 *                                                                          */
int DrDatesIn2FreqSchedule
              (int        NbDatesIn,     /**< (I) Number of dates to check   */
               IRDate      *DatesIn,       /**< (I) User dates                 */
               IRDate       StartDate,     /**< (I) Start of date list         */
               IRDate       EndDate,       /**< (I) End of date list           */
               IRDate       ValueDate,     /**< (I) Value date                 */
               char       PmtFreq,       /**< (I) Frequency of pmts          */
               char       StubConv,      /**< (I)                            */
               char       ObsFreq,       /**< (I) Frequency of obs           */
               int       *NbDatesOut,    /**< (O) Number of found dates      */
               IRDate      *DatesOut)      /**< (O) Found dates                */
{
    int        status = FAILURE;
    int        NbPmtDates = 0; /* always takes values >= 2 but initialized */
    int        NbObsDates = 0; /* here to avoid compiler warning           */
    int        startPmtIdx;
    int        endPmtIdx;
    int        startObsIdx;
    int        endObsIdx;
    int        i, j;

    IRDate       *PmtDateList = NULL;
    IRDate       *ObsDateList = NULL;
    IRDate       startPmtDate;

    /* Input date checks */
    if (EndDate < ValueDate)
    {
        DR_Error("DrDatesIn2FreqSchedule: value date above final maturity!");
        goto FREE_MEM_AND_RETURN;
    }

    /* Build date list to compare against */
    if (DateListFromFreq(StartDate,
                         EndDate,
                         PmtFreq,
                         StubConv,
                         &NbPmtDates,
                         &PmtDateList) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Create observation schedule */
    if (DateListFromFreq(StartDate,
                         EndDate,
                         ObsFreq,
                         'B',
                         &NbObsDates,
                         &ObsDateList) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Check if given dates are observation dates */
    if (DrDatesInSet(NbDatesIn,
                     DatesIn,
                     NbObsDates,
                     ObsDateList) == FAILURE)
    {
        DR_Error("DrDatesIn2FreqSchedule: Input dates are not observation dates!");
        goto FREE_MEM_AND_RETURN;
    }


    /* Find start and end of current period */
    endPmtIdx = 0;
    while ((ValueDate > PmtDateList[endPmtIdx]) &&
           (endPmtIdx < NbPmtDates - 1))
        endPmtIdx++;

    startPmtIdx = endPmtIdx - 1;

    /* Check if the deal is still fwd starting */
    if (endPmtIdx == 0)
    {
        *NbDatesOut = 0;
        status = SUCCESS;

        goto FREE_MEM_AND_RETURN;
    }

    /* Find start and end of observation list in current period */
    startPmtDate = PmtDateList[startPmtIdx];
    startObsIdx = 0;
    while ((startPmtDate > ObsDateList[startObsIdx]) &&
           (startObsIdx < NbObsDates - 1))
        startObsIdx++;

    endObsIdx = 0;
    while ((ValueDate > ObsDateList[endObsIdx]) &&
           (endObsIdx < NbObsDates - 1))
        endObsIdx++;

    /* Compare lists */
    for(i = startObsIdx; i < endObsIdx; i++)
    {
        j = 0;
        while ((ObsDateList[i] > DatesIn[j]) && (j < NbDatesIn - 1))
               j ++;
        if (DatesIn[j] == ObsDateList[i])
        {
            DatesOut[i-startObsIdx] = ObsDateList[i];
        }
        else
        {
            DR_Error ("DrDatesIn2FreqSchedule: Past fixing date %ld is missing!", ObsDateList[i]);
            goto FREE_MEM_AND_RETURN;
        }
    }
    *NbDatesOut = endObsIdx - startObsIdx;


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (ObsDateList, IDATE, 0, NbObsDates-1);
    Free_DR_Array (PmtDateList, IDATE, 0, NbPmtDates-1);

    return(status);

}   /* End of DrDatesIn2FreqSchedule */



/****************************************************************************/
/**                                                                         *
 *   FUNCTION    DrDatesInSet                                               *
 *                                                                          *
 *   Returns SUCCESS if a first set of dates is contained in the second set *
 *   of dates.                                                              *
 *                                                                          *
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       *
 *                                                                          */
int DrDatesInSet
              (int         NbDates1,      /**< (I) Number of dates to check   */
               IRDate const* Dates1,        /**< (I) Dates to check             */
               int         NbDates2,      /**< (I) Number of org dates        */
               IRDate const* Dates2)        /**< (I) Org dates                  */
{

    int     i, j;
    int     status = FAILURE;


    /* Pre-checks */
    if (NbDates1 == 0) return(SUCCESS); /* Nothing to check */

    /* Perform the comparison */
    i = 0;
    for (j = 0; j < NbDates2; j++)
    {
        if (Dates2[j] == Dates1[i]) i++;

        if (i>NbDates1-1) break;
    }

    if (i < NbDates1)
    {
        /*DR_Error ("DrDatesInSet: could not match all dates!");*/
        goto RETURN;
    }


    status = SUCCESS;

    RETURN:

    return(status);

}   /* End of DrDatesInSet() */



/****************************************************************************/
/**                                                                         *
 *   FUNCTION   DrFreeEventList                                             *
 *                                                                          */
void DrFreeEventList(EVENT_LIST  *List)
{
    int i;


    if (List != NULL)
    {
        if (List->Dates != NULL)
        {
            free(List->Dates);
        }
        for (i=0; i<MAXNBEVCURVES; i++)
        {
            if (List->Curve[i] != NULL)
            {
                free(List->Curve[i]);
            }
        }
        free(List);
    }

    return;

}  /* DrFreeEventList */



/*****  AddDateToList  **********************************************/
/**
 *      Allocates memory and add a date to the end of a list of date
 *      Returns SUCCESS or FAILURE
 */

int     AddDateToList
            (int      *NbDates,
             IRDate    **DateList,
             IRDate      NewDate)
{
    int     status = FAILURE;   /* Error status = FAILURE initially */
    IRDate   *DateListLocal = NULL;
    int     Nb;  /* temp var for NbDates */

    if ((NbDates == NULL) || (DateList == NULL)) goto RETURN;

    /* Re-allocate memory */

    Nb = *NbDates + 1;
    if (Nb <= 0) goto RETURN;

    DateListLocal = (IRDate *)DR_REALLOC(*DateList,Nb*sizeof(IRDate));
    if (DateListLocal == NULL)
    {
        DR_Error ("AddDateToList: allocation failure!");
        goto RETURN;
    }

    /* Add the new date */
    DateListLocal[Nb-1] = NewDate;

    *DateList = DateListLocal;
    *NbDates = Nb;

    status = SUCCESS;

RETURN:

    return (status);

} /* AddDateToList */




/*****  SortDateList  **********************************************/
/**
 *      Given a DateList and its size
 *      sort its contents in 1st: date ascending order,
 *                           2nd: value ascending order. (if not NULL)
 *      Returns SUCCESS or FAILURE
 */

int     SortDateList(int    NbDates,
                     IRDate  *DateList,
                     IRDate  *SuppValue)
{
    IRDate   tmpDate, tmpValue;  /* temp storage */
    int    i, j;
    IRDate  *Val = NULL;

#undef  IS_GREATER
#define IS_GREATER(dateA, valA, dateB, valB)          \
                ((dateA > dateB) || ((dateA == dateB) && (valA > valB)))


    if (DateList == NULL) return (FAILURE);

    if (SuppValue == NULL) Val = DateList;  /* trick to sort dates only */
                    else   Val = SuppValue;

    for (j = 1 ; j < NbDates; j++)
    {
        tmpDate = DateList[j];
        tmpValue = Val[j];

        i = j-1;
        while ((i >= 0) &&
               IS_GREATER(DateList[i],Val[i],tmpDate,tmpValue))
        {
            DateList[i+1] = DateList[i];
            Val[i+1]      = Val[i];
            i--;
        }

        DateList[i+1] = tmpDate;
        Val[i+1] = tmpValue;

    }  /* for j */

    return (SUCCESS);

#undef IS_GREATER

} /* SortDateList */


/*********** MergeDateLists ************************************************/
/**                                                                         *
 *      1-merges two date lists,                                            *
 *      2-sorts them                                                        *
 *      3-removes possible duplicates                                       *
 *                                                                          *
 *                                                                          *
 *      Note: that DateList (first input list) is preserved                 *
 *      whereas MergedList (second input list) is modified by the function. *
 *                                                                          *
 *                                                                          *
 ****************************************************************************/

int MergeDateLists(
    int  NbDatesList,    /**< (I)Nb dates in list to be merged with Mergedlist*/
    IRDate *DateList,      /**< (I)List to be merged with MergedList            */
    int  *MergedListSize,/**< (I/O) MergedList size after duplicates removal  */
    IRDate **MergedList)   /**< (I/O) Absorbing list                            */
{
    IRDate    *TempList = NULL;
    int     i,index,TempListSize,status = FAILURE;


    /* Some basic checks first*/
    if (MergedList == NULL ||
        MergedListSize == NULL)
    {
        DR_Error("Invalid pointer inputs to MergeDateLists\n");
        goto RETURN;
    }

    if (NbDatesList < 0 || *(MergedListSize) < 0)
    {
        DR_Error ("Invalid List size inputs to MergedDateLists: "
                    "should be at least 0!\n");
        goto RETURN;
    }

    TempListSize = NbDatesList + (*MergedListSize);
    if (TempListSize == 0)
    {
        status = SUCCESS;   /* nothing to do !*/
        goto RETURN;
    }

    TempList = (IRDate*) DR_REALLOC (*MergedList,TempListSize * sizeof(IRDate));

    if (TempList == NULL)
    {
        goto RETURN;
    }

    *MergedList = TempList;

    for (i = (*MergedListSize); i < TempListSize; i++)
    {
        (*MergedList)[i] = DateList[i - (*MergedListSize)];
    }

    if (SortDateList (TempListSize,(*MergedList),NULL) == FAILURE)
    {
        goto RETURN;
    }

    /* now remove duplicates*/
    index = 1;
    for (i =1 ;i < NbDatesList + (*MergedListSize); i++)
    {
        if ((*MergedList)[i] > (*MergedList)[index - 1])
        {
            (*MergedList)[index] = (*MergedList)[i];
            index++;
        }
        else
            if ((*MergedList)[i] == (*MergedList)[index - 1])
            {
                TempListSize--;
            }
        else
        {
            goto RETURN;
        }
    }/* for i*/

    if (TempListSize < NbDatesList + (*MergedListSize))
    {
        TempList = (IRDate*) DR_REALLOC(*MergedList,TempListSize * sizeof(IRDate));

        if (TempList == NULL)
        {
            goto RETURN;
        }

        *MergedList = TempList;
    }

    *MergedListSize = TempListSize;

    status = SUCCESS;

RETURN:

        return (status);

}/* MergeDateLists */


/*****  GetDLOffset  ***************************************************/
/**
 *      Given a datelist and a targetDate,
 *      returns the offset of the element nearest to the targetdate
 *      if mode =
 *          CbkEXACT:  returns INVALID_DATE if no matching targetdate is found
 *
 *          CbkLOWER:  returns the nearest offset that is LOWER or equal;
 *                     or INVALID_DATE if all dates > targetdate
 *          CbkHIGHER: returns the nearest offset that is HIGHER or equal;
 *                     or INVALID_DATE if all dates < targetdate
 */

int     GetDLOffset(int         NbDates,
                    IRDate const* DL,
                    IRDate        targetDate,
                    int         mode)
{
    int  offset = INVALID_DATE;
    /* int  i = 0; */

    if ((NbDates <=0) || (DL == NULL)) goto RETURN;

    /* find the largest i s.t. DL[i] <= targetDate
     * the i leaving the loop is always valid [0 .. NbDates-1]
     *
     * this final "largest" i may be:
     * ( >0): then DL[i] <= targetDate
     * (==0): then DL[k] >  targetDate for all 1<=k<=NbDates-1
     *             but we don't know about DL[0]
     */

    /*
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
    } 
    */

   
    /* import from srm3 date.c*/
    if (mode == CbkEXACT)
    {
        if (DL[0]         > targetDate) return -999;
        if (DL[NbDates-1] < targetDate) return -999;
    }
    else if (mode == CbkLOWER)
    {
        if (DL[0]         > targetDate) return -999;
        if (DL[NbDates-1] < targetDate) return NbDates-1;
    }
    else if (mode == CbkHIGHER)
    {
        if (DL[0]         > targetDate) return 0;
        if (DL[NbDates-1] < targetDate) return -999;
    }
    else return -999;

    /* trap the case when targetDate = last date */
    /* this is a special case becasue of the downward roundoff in calculating i_mid */
    if (targetDate == DL[NbDates-1]) return NbDates-1;

    {
        int i_begin, i_end, i_mid, i_step;
        i_begin = 0; 
        i_end = NbDates-1; 
        i_mid = (i_begin + i_end)/2; 
        i_step = 1;

        while (i_step != 0) 
        {
            if (DL[i_mid] == targetDate) return i_mid;
            else if (targetDate > DL[i_mid]) 
            {
                i_begin = i_mid;
                i_mid = (i_begin + i_end)/2;

                /* i_step = change in i_mid */
                i_step = i_mid - i_begin;
            } 
            else 
            {
                i_end = i_mid;
                i_mid = (i_begin + i_end)/2;

                /* i_step = change in i_mid */
                i_step = i_end - i_mid;
            }
        }

        switch (mode)
        {
            case CbkEXACT:
                return -999;
            case CbkLOWER:
                return i_begin;
            case CbkHIGHER:
                return i_end;
        }

        return -999;

    }

RETURN:

    return (offset);

} /* GetDLOffset */


/****************************************************************************/
/**                                                                         *
 *   FUNCTION    hasStub                                                    *
 *                                                                          *
 *   Boolean function which returns 0 if the StartDate and EndDate define   *
 *   an integer number (>=0) of periods (i.e. no stub), 1 if there exists   *
 *   a stub.                                                                *
 *   Assumes StartDate > EndDate and Freq = (A,S,Q,M,W)                     *
 *                                                                          */

 int hasStub (IRDate       StartDate,     /**< (I) Start date                   */
              IRDate       EndDate,       /**< (I) End date list                */
              char       Freq)          /**< (I) Frequency                    */
{
    return hasStubLoc(StartDate, EndDate, Freq, 'B');
}


/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    hasStubLoc                                                 * 
 *                                                                          * 
 *   Boolean function which returns 0 if the StartDate and EndDate define   *  
 *   an integer number (>=0) of periods (i.e. no stub), 1 if there exists   * 
 *   a stub. The function accounts for stub location ('F' or 'B')           * 
 *   Assumes StartDate > EndDate and Freq = (A,S,Q,M,W)                     * 
 *                                                                          */

int hasStubLoc (IRDate       StartDate,     /* (I) Start date */
		IRDate       EndDate,       /* (I) End date list */
		char       Freq,          /* (I) Frequency    */
		char       StubLoc)       /* (I) Stub location ('F' or 'B') */
{
  char    IntvalPerType = 'A'; /* Period type corresponding to frequency  */
  long    IntvalNbPers  = 0L;  /* Nb of prds of type above corr to freq   */
  IRDate    AuxDate;             /* Auxilary date                           */
  int     hasStub;             /* Boolean flag for return value           */
	long    i;

	if ((StubLoc != 'F')&&(StubLoc !='B'))
		return 0;

    switch(Freq)
    {
        case 'A': IntvalNbPers   = 12L;
                  IntvalPerType  = 'M';
                  break;
        case 'S': IntvalNbPers   = 6L;
                  IntvalPerType  = 'M';
                  break;
        case 'Q': IntvalNbPers   = 3L;
                  IntvalPerType  = 'M';
                  break;
        case 'M': IntvalNbPers   = 1L;
                  IntvalPerType  = 'M';
                  break;
        case 'W': IntvalNbPers   = 7L;
                  IntvalPerType  = 'D';
                  break;
        default:  return 0;
    }


    if (StubLoc == 'B') 
	{
        i = 0;
		AuxDate = StartDate;
        while (AuxDate < EndDate)
		{
			i++;
            if (IntvalPerType == 'M')
			{
                AuxDate = Nxtmth(StartDate,
                             i*IntvalNbPers,
                             1L);
			}
            else
			{
                AuxDate = Nxtday(StartDate,
                             i*IntvalNbPers);
			}

		}

        hasStub = (AuxDate == EndDate) ? 0 : 1;
	}
	else /* StubLoc = 'F' */
	{
        i = 0;
		AuxDate = EndDate;
        while (AuxDate > StartDate)
		{
			i++;
            if (IntvalPerType == 'M')
			{
                AuxDate = Nxtmth(EndDate,
                             -i*IntvalNbPers,
                             1L);
			}
            else
			{
                AuxDate = Nxtday(EndDate,
                             -i*IntvalNbPers);
			}

		}

        hasStub = (AuxDate == StartDate) ? 0 : 1;
	}

    return(hasStub);

}   /* End of hasStubLoc */


unsigned
ExpandDateSchedule(unsigned sSize, IRDate const* sSched, unsigned tSize, IRDate* tSched, char frequency)
{
    unsigned cnt, ret, i;

    for (i = 0; i < sSize && i < tSize; ++i)
        tSched[i] = sSched[i];

    cnt = i;

    if (strchr("IDWMQSA", frequency) && frequency != 'I')
    {
        IRDate begDate = sSched[0];
        IRDate endDate = sSched[sSize-1];
        IRDate curDate = begDate;

        while (curDate < endDate && cnt < tSize)
        {
            DrDateFwdAny(curDate, 1, frequency, 'F', &curDate);

            if (curDate < endDate) 
                tSched[cnt++] = curDate;
        }
    }

    SortDateList(cnt, tSched, 0); 

    ret = 0;
    for (i = 0; i < cnt; ++i)
    {
        if (i>0 && tSched[i] == tSched[i-1])
            continue;
        tSched[ret++] = tSched[i];
    }

    return ret;
}

/* Converts an IRDate into a string in "DD-MON-YYYY" format. */
void StringFromIRDate(IRDate date, char *string)
{
    long mm, dd, yyyy;
    
    Dsplit(date, &mm, &dd, &yyyy);
    sprintf(string, "%02ld-%s-%04ld", dd, MonthNames[mm-1], yyyy);
}


/** Checks validity of a date                        * 
 *  returns 0 if all is well                         * 
 *  returns 1 if bad year                            * 
 *  returns 2 if bad month                           * 
 *  returns 3 if bad day                             */
int Dateok(IRDate date)
{
    long *daysin;
    long mm,dd,yy;

    Dsplit(date,&mm,&dd,&yy);
    if (yy < 0 ||yy > 3000)
        return 1; /* bad year */
    if (mm < 1 || mm > 12)
        return 2; /* bad month */
    /* see if leap year */
    daysin = Isleap(yy) ? leap : noleap;
    if (dd < 1 || dd > daysin[mm])
        return 3; /* bad month */
    /* all is well if we are here */
    return 0;
}


/* Converts an Excel date as long (1L = 01-Jan-1900) to a date */
int
DrDateFromExcel(long xlDate, IRDate *date)
{
    long    julDate;

    julDate = DrExcel2Julday(xlDate);   /* days offset */

    return DrJulday2Date(julDate, date);
}


/* Converts a date to an Excel date as long */
int
DrDateToExcel(IRDate date, long *xlDate)
{
    long    julDate;    /* Julian date */

    if(DrDate2Julday(date, &julDate) != SUCCESS)
        return (FAILURE);
    else {
        *xlDate = DrJulday2Excel(julDate);  /* days offset */
        return (SUCCESS);
    }
}




/**-------------------------------------------------------------
 * Returns Julian day number corresponding to mm, dd, and yyyy.
 */

int
DrJulday(int mm, int id, int iyyy, long *julDate)
{
    int ja, jy, jm; 
    long    IGREG  = (15+31L*(10+12L*1582));

    /* Convert to Julian day
     */
    if (iyyy == 0) 
    {
        DR_Error("DrJulday: there is no year zero.\n");  /* instead of DR_Error2 defined in masim_s */
        return (FAILURE);
    }

    if (iyyy < 0) ++iyyy;
    if (mm > 2) 
    {
        jy = iyyy;
        jm = mm+1;
    } 
    else 
    {
        jy = iyyy-1;
        jm = mm+13;
    }
    *julDate = (long) (floor(365.25*jy)+floor(30.6001*jm)+id+1720995);
    if (id+31L*(mm+12L*iyyy) >= IGREG) 
    {
        ja = (int)(0.01*jy) ;
        *julDate += 2-ja+(int) (0.25*ja);
    }

    return(SUCCESS) ;
}



/**-------------------------------------------------------------
 * Returns Julian day number corresponding to "DRDate" in format 20020325.
 */
int DrDate2Julday(IRDate drDate, long *julDate)
{
    long    mm, dd, yyyy;

    Dsplit(drDate, &mm, &dd, &yyyy);

    /* Convert to Julian day
     */
    return DrJulday(mm, dd, yyyy, julDate); 

}


/**-------------------------------------------------------------
 * Returns DR date corresponding to a Julian day number "julDate".
 */
int DrJulday2Date(long julDate, IRDate *drDate)
{
    long int    ja, jalpha, jb, jc, jd, je ;
    int     mm, id, iyyy ;
const   long        IGREG  = 2299161 ;

    /* Split Julian day into dd, mm, yyyy
     */
    if (julDate >= IGREG) 
    {
        jalpha = (long int) (((float) (julDate-1867216)-0.25)/36524.25);
        ja = julDate+1+jalpha - (long)(0.25*jalpha);
    } 
    else 
    {
        ja = julDate;
    }

    jb = ja+1524;
    jc = (long int) (6680.0 + ((float) (jb-2439870)-122.1)/365.25);
    jd = (long int) (365L*jc + (0.25*jc));
    je = (long int) ((jb-jd) / 30.6001);
    id = (int) (jb - jd - (int)(30.6001*je));
    mm = (int) (je - 1) ;
    if (mm > 12) mm -= 12;
    iyyy = (int)(jc - 4715);
    if (mm > 2) --(iyyy);
    if (iyyy <= 0) --(iyyy);

    /* Convert to DR date
     */
    *drDate = iyyy*10000 + mm*100 + id;

    *drDate = IRDateFromYMDDate(*drDate);

    return(SUCCESS) ;
}



/*-------------------------------------------------------------
 * Convert an Excel date to to Julian date.
 * The offset between Julian day and Excel day is supplied 
 * by the API through global variable JULIAN_EXCEL_OFFSET.
 */

long DrJulday2Excel(long julDate)
{
    return (julDate - JULIAN_EXCEL_OFFSET);
}



long DrExcel2Julday(long xlDate)
{
    return (xlDate + JULIAN_EXCEL_OFFSET);
}




/************************************************************************************
**
** IMPORTANT:  READ COMMENTS AT TOP OF FILE BEFORE ADDING/MODIFYING CODE TO THIS FILE. 
**
************************************************************************************/
