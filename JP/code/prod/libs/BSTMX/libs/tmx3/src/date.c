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
#include "tmx123head.h"


#define HOLIDAY_FILE "HOLIDAY"
#define OKAY 0
#define ERR -1


long noleap[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
long leap[13] =  {0,31,29,31,30,31,30,31,31,30,31,30,31};
long cumdays[12] = {0,31,59,90,120,151,181,212,243,273,304,334};


void Dsplit(long date_i, long *mm_o, long *dd_o, long *yy_o)
{
    /*  Returns month, day, year from an integer in                */
    /*  the format YYY(Y)MMDD as in 840701 or 20011115 or 1120104  */

    *yy_o = date_i/10000;
    *mm_o = (date_i - *yy_o*10000)/100;
    *dd_o = date_i - *yy_o*10000 - *mm_o*100;
    return;
}

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


int Isimm(long    date)
{

        /*  Returns TRUE if the date is an IMM date    */

        long  day;
        long  mth;
        long  year;

        long  thirdWed;
        int   status = FALSE;

        
        Dsplit(date, &mth, &day, &year);

        if ((mth==3) || (mth==6) || (mth==9) || (mth==12)) 
        {
        
            thirdWed = ThirdWed(mth, year);
            if (thirdWed == date)
            {
                status = TRUE;
            }
        }

        return (status);
            
}


long   ThirdWed(long mth, long year)
{

    long  date;

    date = Datepack(mth, 1L, year);

    while (Dayofwk(date) != 3)
    {
        date++;
    }

    /* Now we are at the 1st Wednesday */
    /* and must increment by two weeks */
    date = date + 14L;

    return(date);

}


int Holiday(long date_i)
{
    /*  This routine checks whether the date passed in is a holiday  */
    /*  or not according to the HOLIDAYS file. It returns (0) if     */
    /*  it is not and (-1) if it is. The format is YYYYMMDD.         */

    FILE *fp1;
    int xflag = ERR;
    int rcode = ERR;
    long holiday;
    char line1[82];

    if ((fp1=fopen(HOLIDAY_FILE,"r")) == NULL) {
        fprintf(stdout,"Can't Open file %s\n",HOLIDAY_FILE);
        return(rcode);
    }
    while (xflag != OKAY) {
        if (fgets(line1,80,fp1) == NULL) {
            xflag = OKAY;
            rcode = OKAY;
        }
        if (xflag != OKAY) {
            line1[8]='\0';
            holiday = atol(line1);
            if (holiday > date_i) {
                    rcode=OKAY;
                    xflag=OKAY;
            }
            else if (holiday == date_i) {
                rcode = ERR;;
                xflag = OKAY;
            }
        }
    }
    fclose(fp1);
    return(rcode);
}


int Dateok(long date)
{
    /*  Checks for validity of YYYYMMDD formatted dates  */
    /*  returns 0 if all is well                         */
    /*  returns 1 if bad year                            */
    /*  returns 2 if bad month                           */
    /*  returns 3 if bad day                             */

    long *daysin;
    long mm,dd,yy;

    Dsplit(date,&mm,&dd,&yy);
    if (yy < 0 ||yy > 3000)
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


int Date_CheckAndReport(
    long                date)               /* (I) Date */
{
    int                 rc;                 /* REsult code */

    /* Check it */
    rc = Dateok(date);

    /* Report it */
    switch (rc)
    {
    case 0:     /* OK */
                break;
    case 1:     /* Bad year */
                DR_Error ("Date_CheckAndReport: Date %ld has a bad year.\n", date);
                break;
    case 2:     /* Bad month */
                DR_Error("Date_CheckAndReport: Date %ld has a bad month.\n", date);
                break;
    case 3:     /* Bad day */
                DR_Error ("Date_CheckAndReport: Date %ld has a bad day.\n", date);
                break;
    default:    /* Unknown */
                DR_Error ("Date_CheckAndReport: %ld is a bad date.\n", date);
                break;
    }

    /* Done */
    return rc;
}


long Y2toy4(long year_i)
{
    /*  Converts a 2 digit year (84) to a 4 digit year (1984);       */
    /*  assumes that 2 digit years less than 50 are after year 2000  */

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


long Datepack(long mm_i, long dd_i, long yy_i)
{
    /*  Packs mm_i,dd_i,yy_i into YYYYMMDD format  */
    /*  calls Y2toy4 on yy_i before packing        */

    long y4,packed_o;

    y4 = Y2toy4(yy_i);
    packed_o = y4*10000 + mm_i*100 + dd_i;
    return (packed_o);
}


long Y4toy2(long year_i)
{
    /*  Converts a 4 digit year to a 2 digit year;  */
    /*  years after 99 are 00, 01, etc.             */

    long y2_o;

    y2_o = year_i%100;
    return (y2_o);
}


void Y2date_str(long date, char *string)
{
    /*  Convert YYYYMMDD to MM/DD/YY  */

    long mm,dd,yy;

    Dsplit(date,&mm,&dd,&yy);
    yy = Y4toy2(yy);
    sprintf(string,"%02ld/%02ld/%02ld",mm,dd,yy);

    return;
}


long eval_date(char *datest)
{
    /*  Convert MM/DD/YY to YYYYMMDD  */

    long date;
    long year, month, day;

    year = (datest[6] - '0') * 10 + (datest[7] - '0');
    year = (year > 73) ? 1900 + year : 2000 + year;
    month = (datest[0] - '0') * 10 + (datest[1] - '0');
    day = (datest[3] - '0') * 10 + (datest[4] - '0');
    date = year * 10000 + month * 100 + day;

    return(date);
}


long eval_date2(char *datest)
{
    /*  Convert DD-MMM-YYYY to YYYYMMDD  */

    long date;
    long year, month, day;
    char string[4];

    year  = (datest[7] - '0') * 1000 + (datest[8] - '0') * 100 
            + (datest[9] - '0') * 10 + (datest[10] - '0');
    day = (datest[0] - '0') * 10 + (datest[1] - '0');
        
    string[0] = datest[3];
    string[1] = datest[4];
    string[2] = datest[5];
    string[3] = '\0';
        
    if (strcmp (string, "Jan") == 0)
        month = 1;
    else if (strcmp (string, "Feb") == 0)
        month = 2;
    else if (strcmp (string, "Mar") == 0)
        month = 3;
    else if (strcmp (string, "Apr") == 0)
        month = 4;
    else if (strcmp (string, "May") == 0)
        month = 5;
    else if (strcmp (string, "Jun") == 0)
        month = 6;
    else if (strcmp (string, "Jul") == 0)
        month = 7;
    else if (strcmp (string, "Aug") == 0)
        month = 8;
    else if (strcmp (string, "Sep") == 0)
        month = 9;
    else if (strcmp (string, "Oct") == 0)
        month = 10;
    else if (strcmp (string, "Nov") == 0)
        month = 11;
    else if (strcmp (string, "Dec") == 0)
        month = 12;
    else         
        month = 0;                          /* Error */
                
    date = year * 10000 + month * 100 + day;

    return(date);
}


long Dayofwk(long date)
{
    /*  Calculates day of week (0-6) of YYYYMMDD formatted date  */

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
    return (dd%7); /* 0 to 6*/
}


/*****************************************************************************/
/*                                                                           */
/*  FUNCTION   Days360                                                       */
/*                                                                           */
/*     Counts days between two dates following the 30/360 covention as       */
/*     specified in  the ISDA rule book. The method below agrees  with       */
/*     the Analytics Library and with the routines implemented by  the       */
/*     STIRT group (thanks to Neill Penney for useful discussions!)          */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
long Days360(long date1_i, long date2_i)
{
    long month1,day1,year1,month2,day2,year2,days;

    /* Trouble shooting for 2/28 coupon */
    if (date1_i == date2_i) 
    {
        return(0L);
    }

    /* split dates into components */
    Dsplit(date1_i,&month1,&day1,&year1);
    Dsplit(date2_i,&month2,&day2,&year2);


    if (day1 == 31)
        day1= 30;

    if ((day2 == 31) && (day1 == 30))
        day2 = 30;
    
    days = (year2-year1)*360 + (month2-month1)*30 + day2-day1;

    return (days);
}


long Months360(long date1_i, long date2_i)
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


long Daysact(long date1_i, long date2_i)
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


long Nxtday(long date, long days)
{
    /*  Returns a date in YYYYMMDD format based on moving        */
    /*  forward or backwards a number of calendar days.          */

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
    return(Datepack(nxtm,edays,nxty));
}


long Nxtmth(long date, long mths, long eom)
{
    /*  Returns a date in YYYYMMDD format based on moving        */
    /*  forward or backwards (mths) from (date).                 */
    /*  if eom > 0 the returned does not spill over to the next  */
    /*  month. That is moving 1m from the 31rst of march gives   */
    /*  the 30th of april. If eom=0 it does spill over: 1m from  */
        /*  the 31rst of march is the 1rst of may.           */
    /*  a positive entry for mths => move forward                */
    /*  a negative entry for mths => move backward.              */

    long *daysin,*daysout;
    long mm,dd,yy;
    long sign,nxtmm,nxtdd,nxtyy,nxtdte;

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
        
    return (nxtdte);
}




long  Nxtimm(long  date, long nbPeriods)
{


    /*  a positive entry for nbPeriods => move forward     */
    /*  a negative entry for nbPeriods => move backward.   */


    long   immdate;
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




/**************************************************************************/
/*                                                                        */
/* FUNCTION    DrDayCountFraction                                         */
/*                                                                        */
/* Calculates the daycount fraction between two dates according to one of */
/* three methods: Actual/365Fixed, Actual/360 and 30/360.                 */
/*                                                                        */
/* RETURNS:  Success or Failure                                           */
/*                                                                        */
/**************************************************************************/

int  DrDayCountFraction(long     Date1,      /*(I) Start date             */
                        long     Date2,      /*(I) End date               */
                        char     Conv,       /*(I) Convention (A,0,3,5)   */
                        double  *Fraction)   /*(O) Corresponding fraction */
{
    long    Ndays;

    switch (Conv)
    {
        case 'A':
            DR_Error("DrDayCountFraction: Convention Act/Act not supported!");
            return(FAILURE);

        case '5':
            Ndays = Daysact(Date1,
                            Date2);
            *Fraction = (double)Ndays/365.0;
            break;
        case '0':
            Ndays = Daysact(Date1,
                            Date2);
            *Fraction = (double)Ndays/360.0;
            break;  
        case '3':
            Ndays = Days360(Date1,
                            Date2);
            *Fraction = (double)Ndays/360.0;
            break;  
        default:
            DR_Error("DrDayCountFraction: "
                    "Convention for day count fraction not supported!");
            return(FAILURE);
    }

    return(SUCCESS);

}




long Nxtwkday(long date, long advance)
{
    /*  This routine returns a new date after advancing or going  */
    /*  backward a number of weekdays (Monday - Friday).          */
    /*  The format for date is YYYYMMDD.                          */

    long dow,work_date;

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


long Nxtbusday(long date, long advance)
{
    /*  This routine returns a new date after advancing or going  */
    /*  backward a number of business days (Monday - Friday &     */
    /*  holidays). It checks for holidays in HOLIDAY file.        */
    /*  The format for date is YYYYMMDD.                          */

    long work_date;

    work_date = date;

    if (advance == 0)
        return(date);
    else if (advance < 0) {  /* from date, go backward */
        while (advance != 0) {
            advance++;
            work_date = Nxtwkday(work_date,(long)-1);
            while(Holiday(work_date))
                work_date=Nxtwkday(work_date,(long)-1);
        }
    }
    else {  /* from date, go forward */
        while (advance != 0) {
            advance--;
            work_date = Nxtwkday(work_date,(long)1);
            while(Holiday(work_date))
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
        INTERP_TYPE interpType,  /**< (I) Interpolation type to use*/
        int         NbInpDates,  /**< (I) Nb of dates input directly   */
        long       *InpDates,    /**< (I) Dates given directly by user */
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
                long   *InpDates,  /**< (I) Dates given directly by user */
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
        EventListLocal->Dates = (long *) DR_Array (LONG, 0, NbInpDates-1);
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
        long  CurrDate;
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
        long  CurrDate;
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
/*                                                                          */
/*   FUNCTION    DrDateFwdAny                                               */
/*                                                                          */
/*   Starting from a given date, this function calculates the date obtained */
/*   by incrementing a number of periods.                                   */
/*                                                                          */
int     DrDateFwdAny (  long    StartDate, /* (I) Start date          */
                        int     NbPers,    /* (I) Number of periods   */
                        char    Freq,      /* (I) Frequency           */
                        char    FwdOrBwd,  /* (I) Forward or backward */
                        long    *DateOut)  /* (O) Calculated date     */
{
    char    IntvalPerType;     /* Period type corresponding to frequency   */
    long    IntvalNbPers;      /* Nb of periods of type above corr to freq */
    int     Factor;    
    int     status = FAILURE;



    if (FwdOrBwd == 'F') 
    {
        Factor = 1;
    }
    else if (FwdOrBwd == 'B')
    {
        Factor = -1;
    }
    else
    {
        DR_Error ("DrDateFwdAny: Fourth argument must be F(fwd) or B(bwd)!");
        goto RETURN;
    }


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
                  DR_Error ("DrDateFwdAny: Unrecognised frequency!");
                  goto RETURN;

    }  /* switch */


     if (IntvalPerType == 'M')
     {
         *DateOut =  Nxtmth(StartDate,
                            (long)(Factor*IntvalNbPers),
                             1L);
     }
     else if (IntvalPerType == 'I')
     {
         *DateOut =  Nxtimm(StartDate,
                            (long)(Factor*IntvalNbPers));
     }
     else
     {
         *DateOut =  Nxtday(StartDate,
                            (long)(Factor*IntvalNbPers));
     }


    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error("DrDateFwdAny: Failed.");
    }

    return(status);

}   /* End of DrDateFwdAny() */



/****************************************************************************/
/*                                                                          */
/*   FUNCTION    DateListFromFreq                                           */
/*                                                                          */
/*   Generates a date list including the start and end dates given and in   */
/*   accordance with the frequency and stub convention passed in.           */
/*                                                                          */
 int DateListFromFreq
              (long       StartDate,     /* (I) Start of date list       */
               long       EndDate,       /* (I) End of date list         */
               char       Freq,          /* (I) Frequency of list        */
               char       StubConv,      /* (I) Stub location Front/back */
               int       *NbDates,       /* (O) Number of dates in list  */
               long     **DateList)      /* (O) List of dates asc order  */
{


    int status = FAILURE;

    char       IntvalPerType;   /* Period type corresponding to frequency   */
    long       IntvalNbPers;    /* Nb of periods of type above corr to freq */
    
    long       AuxDate;                 /* Variables for local use          */
    int        Factor;
    int        IndexStart;

    long      *DateListLocal = NULL;    /* Local copy of output list        */
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
                  DR_Error ("DateListFromFreq: Unrecognised frequency!");
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
                    "an integer number of periods, if stubs are not allowed!");
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
    DateListLocal = (long *) DR_Array (LONG, 0, NbDatesLocal-1);
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
        Free_DR_Array (DateListLocal, LONG, 0, NbDatesLocal-1);
        DR_Error("DateListFromFreq: Failed.");
    }

    return(status);

}   /* End of DateListFromFreq() */



/****************************************************************************/
/*                                                                          */
/*   FUNCTION    DrDatesInSchedule                                          */
/*                                                                          */
/*   Returns SUCCESS if all dates passed in are contained in the schedule   */
/*   specified by a start date, an end date, a frequency and a stub conv.   */
/*                                                                          */
/*   DatesIn must be an arrays of ordered, strictly increasing dates.       */
/*                                                                          */
 int DrDatesInSchedule
              (int        NbDatesIn,     /* (I) Number of dates to check   */
               long      *DatesIn,       /* (I) User dates                 */
               long       StartDate,     /* (I) Start of date list         */
               long       EndDate,       /* (I) End of date list           */
               char       Freq,          /* (I) Frequency of list          */
               char       StubConv)      /* (I)                            */
{


    int status = FAILURE;

    int        NbDatesLocal;
    long      *DateListLocal = NULL;    /* For local use only */
    
    int        i, j;
    


    /* Pre-checks */
    if (NbDatesIn == 0)
        return(SUCCESS); /* Nothing to check */

    if (DatesIn[0] < StartDate)
    {
        DR_Error("DrDatesInSchedule: first supplied date precedes start "
                    "of date list!");
        return(FAILURE);
    }
    if (DatesIn[NbDatesIn-1] > EndDate)
    {
        DR_Error("DrDatesInSchedule: last supplied date is beyond end "
                    "of date list!");
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
        DR_Error("DrDatesInSchedule: could not match all dates!");
        goto FREE_MEM_AND_RETURN;
    }


    status = SUCCESS;

    FREE_MEM_AND_RETURN:
    
    Free_DR_Array (DateListLocal, LONG, 0, NbDatesLocal-1);

    return(status);

}   /* End of DrDatesInSchedule() */



/****************************************************************************/
/*                                                                          */
/*   FUNCTION    DrSameDateSchedules                                        */
/*                                                                          */
/*   Returns SUCCESS if dates passed in are the same as dates generated by  */ 
/*   specifing a start date, an end date, a frequency, a stub conv,         */
/*   arrears reset and a value date, which is used to cut the list before it*/
/*   DatesIn must be an arrays of ordered, strictly increasing dates.       */
/*                                                                          */
 int DrSameDateSchedules
              (int       *NbDatesIn,     /* (I/0) Number of dates to check */
               long      *DatesIn,       /* (I) User dates                 */
               long       StartDate,     /* (I) Start of date list         */
               long       EndDate,       /* (I) End of date list           */
               long       ValueDate,     /* (I) Value date                 */
               char       Freq,          /* (I) Frequency of list          */
               char       StubConv,      /* (I)                            */
               char       Arrears)       /* (I) 'Y' if reset-in-arrears    */
{


    int status = FAILURE;

    int        NbDatesLocal;
    int        NbResetsLocal;
    long      *DateListLocal = NULL;
    long      *ResetListLocal = NULL;
    int        i;
    char    ErrorMsg[MAXBUFF];
    


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
            sprintf (ErrorMsg,
                    "DrSameDateSchedules: Fixing date %ld is missing!",
                     ResetListLocal[i]);
            DR_Error (ErrorMsg);
            goto FREE_MEM_AND_RETURN;
        }
        i++;
    }

    *NbDatesIn = i;
    
    status = SUCCESS;

    FREE_MEM_AND_RETURN:
    
    Free_DR_Array (DateListLocal, LONG, 0, NbDatesLocal-1);

    return(status);

}   /* End of DrSameDateSchedules() */



/****************************************************************************/
/*                                                                          */
/*   FUNCTION    DrSameDateSets                                             */
/*                                                                          */
/*   Returns SUCCESS if dates passed in are the same as dates generated by  */ 
/*   cuting a givne date list up to a value date.                           */
/*   DatesIn must be an arrays of ordered, strictly increasing dates.       */
/*                                                                          */
int DrSameDateSets
              (int       *NbDates1,      /* (I/0) Number of dates to check */
               long      *Dates1,        /* (I) User dates                 */
               int        NbDates2,      /* (I) Number of dates to cmp     */
               long      *Dates2,        /* (I) Cpm dates                  */
               long       ValueDate)     /* (I) Value date                 */
{
    int     i;
    int     status = FAILURE;
    char    ErrorMsg[MAXBUFF];
    


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
                sprintf (ErrorMsg,
                        "DrSameDateSets: Fixing date %ld is missing!",
                         Dates2[i]);
                DR_Error (ErrorMsg);
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
/*                                                                          */
/*   FUNCTION    DrSameDateSetsFlows                                        */
/*                                                                          */
/*   Returns SUCCESS if dates passed in are the same as dates generated by  */ 
/*   cuting a givne date list up to a value date.                           */
/*   DatesIn must be an arrays of ordered, strictly increasing dates.       */
/*                                                                          */
int DrSameDateSetsFlows
              (int       *NbDates1,      /* (I/0) Number of dates to check */
               long      *Dates1,        /* (I) User dates                 */
               int        NbDates2,      /* (I) Number of dates to cmp     */
               long      *Dates2,        /* (I) Cpm dates                  */
               long       ValueDate)     /* (I) Value date                 */
{
    int     i;
    int     status = FAILURE;
    char    ErrorMsg[MAXBUFF];
    


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
                sprintf (ErrorMsg,
                        "DrSameDateSets: Fixing date %ld is missing!",
                         Dates2[i]);
                DR_Error (ErrorMsg);
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
/*                                                                          */
/*   FUNCTION    DrDatesIn2FreqSchedule                                     */
/*                                                                          */
/*   Compares the passed-in date list to a date list generated in the       */
/*   current (with respect to the value date) pmt period according to some  */
/*   obs frequency. The current pmt period is found from the start date,    */
/*   end date, pmt frequency, stub convention and a value date              */
/*                                                                          */
int DrDatesIn2FreqSchedule
              (int        NbDatesIn,     /* (I) Number of dates to check   */
               long      *DatesIn,       /* (I) User dates                 */
               long       StartDate,     /* (I) Start of date list         */
               long       EndDate,       /* (I) End of date list           */
               long       ValueDate,     /* (I) Value date                 */
               char       PmtFreq,       /* (I) Frequency of pmts          */
               char       StubConv,      /* (I)                            */
               char       ObsFreq,       /* (I) Frequency of obs           */
               int       *NbDatesOut,    /* (O) Number of found dates      */
               long      *DatesOut)      /* (O) Found dates                */

{


    int        status = FAILURE;

    int        NbPmtDates = 0; /* always takes values >= 2 but initialized */
    int        NbObsDates = 0; /* here to avoid compiler warning           */
    int        startPmtIdx;
    int        endPmtIdx;
    int        startObsIdx;
    int        endObsIdx;
    int        i, j;

    long       *PmtDateList = NULL;
    long       *ObsDateList = NULL;
    long       startPmtDate;

    char       ErrorMsg[MAXBUFF];
    

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
            sprintf (ErrorMsg,
                     "DrDatesIn2FreqSchedule: Past fixing date %ld is missing!",
                     ObsDateList[i]);
            DR_Error (ErrorMsg);
            goto FREE_MEM_AND_RETURN;
        }
    }
    *NbDatesOut = endObsIdx - startObsIdx;


    status = SUCCESS;

    FREE_MEM_AND_RETURN:
    
    Free_DR_Array (ObsDateList, LONG, 0, NbObsDates-1);
    Free_DR_Array (PmtDateList, LONG, 0, NbPmtDates-1);

    return(status);

}   /* End of DrDatesIn2FreqSchedule */



/****************************************************************************/
/*                                                                          */
/*   FUNCTION    DrDatesInSet                                               */
/*                                                                          */
/*   Returns SUCCESS if a first set of dates is contained in the second set */
/*   of dates.                                                              */
/*                                                                          */
/*   DatesIn must be an arrays of ordered, strictly increasing dates.       */
/*                                                                          */
int DrDatesInSet
              (int        NbDates1,      /* (I) Number of dates to check   */
               long      *Dates1,        /* (I) Dates to check             */
               int        NbDates2,      /* (I) Number of org dates        */ 
               long      *Dates2)        /* (I) Org dates                  */
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
        /* 
         * TMX3: turn this error message off; as the same error message
         * is implemented by all the callers; unlike FIX3, failure of this
         * function is allowed. Note: this funtion should actually return
         * TRUE or FALSE ! 
         */
        /* DR_Error ("DrDatesInSet: could not match all dates!"); */
        goto RETURN;
    }


    status = SUCCESS;

    RETURN:

    return(status);

}   /* End of DrDatesInSet() */



/****************************************************************************/
/*                                                                          */
/*   FUNCTION   DrFreeEventList                                             */
/*                                                                          */
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
/*
 *      Allocates memory and add a date to the end of a list of date
 *      Returns SUCCESS or FAILURE
 */

int     AddDateToList
            (int      *NbDates,
             long    **DateList,
             long      NewDate)
{
    int     status = FAILURE;   /* Error status = FAILURE initially */
    long   *DateListLocal = NULL;
    int     Nb;  /* temp var for NbDates */

    if ((NbDates == NULL) || (DateList == NULL)) goto RETURN;

    /* Re-allocate memory */

    Nb = *NbDates + 1;
    if (Nb <= 0) goto RETURN;

    DateListLocal = (long *)DR_REALLOC(*DateList,Nb*sizeof(long));
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
/*
 *      Given a DateList and its size
 *      sort its contents in 1st: date ascending order, 
 *                           2nd: value ascending order. (if not NULL)
 *      Returns SUCCESS or FAILURE
 */

int     SortDateList(int    NbDates, 
                     long  *DateList,
                     long  *SuppValue)
{
    long   tmpDate, tmpValue;  /* temp storage */
    int    i, j;
    long  *Val = NULL;

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
    long *DateList,      /**< (I)List to be merged with MergedList            */
    int  *MergedListSize,/**< (I/O) MergedList size after duplicates removal  */
    long **MergedList)   /**< (I/O) Absorbing list                            */
{
    long    *TempList = NULL;
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

    TempList = (long*) DR_REALLOC (*MergedList,TempListSize * sizeof(long));

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
        TempList = (long*) DR_REALLOC(*MergedList,TempListSize * sizeof(long));
    
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


/****************************************************************************/
/*                                                                          */
/*   FUNCTION    hasStub                                                    */
/*                                                                          */
/*   Boolean function which returns 0 if the StartDate and EndDate define   */ 
/*   an integer number (>=0) of periods (i.e. no stub), 1 if there exists   */
/*   a stub.                                                                */
/*   Assumes StartDate > EndDate and Freq = (A,S,Q,M,W)                     */
/*                                                                          */
 int hasStub (long       StartDate,     /* (I) Start date                   */
              long       EndDate,       /* (I) End date list                */
              char       Freq)          /* (I) Frequency                    */
{

    char    IntvalPerType = 'A'; /* Period type corresponding to frequency  */
    long    IntvalNbPers  = 0L;  /* Nb of prds of type above corr to freq   */
    long    AuxDate;             /* Auxilary date                           */
    int     hasStub;             /* Boolean flag for return value           */

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
    }

    AuxDate = StartDate;
    while (AuxDate < EndDate)
    {
        if (IntvalPerType == 'M')
        {
            AuxDate = Nxtmth(AuxDate,
                             IntvalNbPers,
                             1L);
        }
        else
        {
            AuxDate = Nxtday(StartDate,
                             IntvalNbPers);
        }

    }

    hasStub = (AuxDate == EndDate) ? 0 : 1;

    return(hasStub);

}   /* End of hasStub */


/*****  NearestDate  ***************************************************/
/*
 *      Given a datelist and a targetDate,
 *      returns the offset of the element nearest to the targetdate
 */
int     NearestDate(int         NbDates,
                    long       *DL,
                    long        targetDate)
{
    int    offset = -999;
    int    i;
    double lower, upper;

    if ((NbDates <=0) || (DL == NULL)) goto RETURN;

    /* find the largest i s.t. DL[i] <= targetDate < DL[i+1]    */
    i = GetDLOffset(NbDates,
                    DL,
                    targetDate,
                    CbkLOWER);

    if (i == NbDates - 1)
    {
        offset = i;
    }
    else if (i >= 0)
    {
        lower = Daysact(DL[i],targetDate) / 365.;
        upper = Daysact(targetDate,DL[i+1]) / 365.;

        if (lower <= upper)
        {
            offset = i;
        }
        else
        {
            offset = i + 1;
        }
    }

RETURN:

    return (offset);

} /* NearestDate */


