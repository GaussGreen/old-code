/************************************************************************
 * Module:	DRL
 * Submodule:	TIME
 * File:	
 * Function:	DDate and Time Routines
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <stdarg.h>
#include <limits.h>

#ifdef DRL_CLIB
# include "dateconv.h"
#endif

#include "drlstr.h"		/* DrlStrRemoveChars() */
#include "drltime.h"

#define MIN_YEAR	1602
#define	MAX_YEAR	2070
#define DRL_YEAR_IS_LEAP(year) \
		(((year) % 4 == 0 && (year) % 100 != 0) || ((year) % 400 == 0))

static	char	monthName[14][5] = {"",
			"Jan", "Feb", "Mar", "Apr", "May", "Jun",
		 	"Jul", "Aug", "Sep", "Oct", "Nov", "Dec", 0};
static	int	daysInMonth[]  = { 31, 28, 31, 30, 31, 30,
			   31, 31, 30, 31, 30, 31};
static	int	daysInMonthLeap[] = { 31, 29, 31, 30, 31, 30,
			   31, 31, 30, 31, 30, 31};
static	char	shortDayNames[12][8]
		= {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"};

static  double 	C1_12 = 0.0833333333333 ;
static  double 	C1_30 = 0.0333333333333 ;

/* If double digit year less than
 * this, assume year 
 * belongs to 21 century.
 * as per Y2K Guidelines
 */
#define CROSSOVER_YEAR 60

static	void Dsplit(long date_i, long *mm_o, long *dd_o, long *yy_o);
static	int Isleap(long year);
static	int Isimm(long    date);
static	long   ThirdWed(long mth, long year);
static	int Dateok(long date);
static	long Y2toy4(long year_i);
static	long Datepack(long mm_i, long dd_i, long yy_i);
static	long Y4toy2(long year_i);
static	void Y2date_str(long date, char *string);
static	long eval_date(char *datest);
static	long eval_date2(char *datest);
static	long Dayofwk(long date);
static	long Days360(long date1_i, long date2_i);
static	long Months360(long date1_i, long date2_i);
static	long Daysact(long date1_i, long date2_i);
static	long Nxtday(long date, long days);
static	long Nxtmth(long date, long mths, long eom);
static	long Nxtimm(long  date, long nbPeriods);
static	long Nxtwkday(long date, long advance);

static	int strToMonth(char *cp, int *monthN);






/*f-------------------------------------------------------------
 * Date routines : check date valid.
 *                                                          
 * <br><br>
 * Returns SUCCESS (0) if "aDate" is a valid date between
 * 1900 and 2070, and FAILURE (1) otherwise.
 */

int
DrlDDateCheckValid(DDate aDate)
{
static	char	routine[] = "DrlDDateCheckValid";
	int		status = FAILURE;
	int		yy, mm, dd;
#ifdef  DRL_CLIB
	TMonthDayYear	mdy;

	if (GtoDateToMDY(aDate, &mdy) != SUCCESS) goto done;
	dd = (int) mdy.day;
	mm = (int) mdy.month;
	yy = (int) mdy.year;
#else
	IF_FAILED_DONE( DrlDDateSplit(aDate, &mm, &dd, &yy));
#endif

	if((yy < MIN_YEAR) || (yy > MAX_YEAR)) goto done;
	if((mm < 1) || (mm > 12)) goto done;
	if ((yy % 4 == 0 && (yy % 100 != 0)) || (yy % 400 == 0)) {
		/* year is leap */
		if((dd < 0) || (dd > daysInMonthLeap[mm-1])) goto done;
	} else {
		if((dd < 0) || (dd > daysInMonth[mm-1])) goto done;
	}

	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: invalid date %ld (m=%d,d=%d,y=%d)\n",
		routine, aDate, mm, dd, yy);
	}
	return(status);
}

/*f-------------------------------------------------------------
 * Date routines : split date.
 *                                                          
 * <br><br>
 * Breaks "aDate" in months "mm", days "dd", and years "yy".
 */

int
DrlDDateSplit(DDate aDate, int *mm, int *dd, int *yy)
{
#ifdef	DRL_CLIB
	int	errCode =0;
	TMonthDayYear	mdy;

	errCode = GtoDateToMDY(aDate, &mdy);
	*dd = (int) mdy.day;
	*mm = (int) mdy.month;
	*yy = (int) mdy.year;
	return(errCode);
#else
	long	mm_l, dd_l, yy_l;
	Dsplit(aDate, &mm_l, &dd_l, &yy_l);
	*mm = (int) mm_l;
	*dd = (int) dd_l;
	*yy = (int) yy_l;
	return(SUCCESS);
#endif
}


/*f-------------------------------------------------------------
 * Date routines : make date.
 *                                                          
 * <br><br>
 * Creates and returns a date corresponding to
 * month "m", day "d" and year "y". 
 */

int
DrlDDateMake(int m, int d, int y, DDate *date)
{
#ifdef	DRL_CLIB
	TMonthDayYear	mdy;

	mdy.month = (long) m;
	mdy.day = (long) d;
	mdy.year = (long) y;

	return GtoMDYToDate(&mdy, date);
#else
	*date = Datepack(m, d, y);
	return(SUCCESS);
#endif

}

/*f-------------------------------------------------------------
 * Date routines : scan date from string.
 *                                                          
 * <br><br>
 * Scan the string <i> string</i> for a date. Recognized formats
 * are <i> mm/dd/yy</i>, <i> mm/dd/yyyy</i>, or <i> dd-MMM-yyyy</i>
 * or <i> yyyy.mm.dd</i>.
 * Puts the result in <i> aDate</i>. Returns 0 if success.
 */

int
DrlDDateScan(char *string, DDate *aDate)
{
static	char	routine[] = "DrlDDateScan";
	int	status = FAILURE;

	char	dateCopy[32];
	int	mm, dd, yy;
	char	*cp;
	int	isMonth;


	/* special case of yyyy.mm.dd not handled below */
	if (strchr(string, '.') != NULL) {
		if (sscanf(string, "%d.%d.%d", &yy, &mm, &dd) != 3) {
			DrlErrMsg("%s: can't read date `%s'.\n", routine);
			return(FAILURE);
		}
		sprintf(dateCopy, "%04d%02d%02d", yy, mm, dd);
	} else {
		strncpy(dateCopy, string, sizeof(dateCopy)-1);
	}

	dateCopy[sizeof(dateCopy)-1] = '\0';
	for (cp = dateCopy; *cp != '\0'; cp++) {
		if (isspace(*cp)) *cp = '\0';
	}

#ifdef	DRL_CLIB
	/* Use Gto routine */
	return GtoStringToDate(dateCopy, aDate);
#else

        /* Support the international YYYYMMDD format.  This is needed
           for holiday files, regression test input and output files,
           etc. */
        if( strpbrk(dateCopy, "-/. ") == NULL )   /* Check for delimiter */
        {
            char smallBuf[5];
    
            if( strlen(dateCopy) != 8 )
                goto done;
            cp = dateCopy;
            
            /* Create the year. */
            smallBuf[0] = *cp++;
            smallBuf[1] = *cp++;
            smallBuf[2] = *cp++;
            smallBuf[3] = *cp++;
            smallBuf[4] = '\0';
            yy = atol(smallBuf);
            if( yy == 0 )
                goto done;
            /* No adjustments are needed for the year. */
    
            /* Create the month. */
            smallBuf[0] = *cp++;
            smallBuf[1] = *cp++;
            smallBuf[2] = '\0';
            mm = atol(smallBuf);
            if( mm == 0 )
                goto done;
    
            /* Create the day. */
            smallBuf[0] = *cp++;
            smallBuf[1] = *cp++;
            /* smallBuf[2] = '\0'; */
            dd = atol(smallBuf);
            if( dd == 0 )
                goto done;
        }
        else                                   /* Parse delimited format */
        {
            /* This block parses the following formats:
                  MM-DD-YYYY, Mth-DD-YYYY, DD-Mth-YYYY   */ 
            /* Get month (or possibly day)
             */
            isMonth = FALSE;
            cp = strtok(dateCopy, "-/. ");
            if (cp == NULL)
                goto done;
            if (sscanf(cp,"%ld",&mm) == 0)  /* If not integer try string */
            {
                if (strToMonth(cp,&mm) == FAILURE)
                {
                    goto done;
                }
                isMonth = TRUE;    /* Make sure input doesn't contain 2 months */
            }
    
            /* Get day (or possibly month)
             */
            cp = strtok(NULL, "-/. ");
            if (cp == NULL)
                goto done;
            if (sscanf(cp,"%ld",&dd) == 0)   /* If not integer try string */
            {
                dd = mm;              /* Month was really day */
                if (isMonth == TRUE || strToMonth(cp,&mm) == FAILURE)
                {
                    goto done;
                }
            }
            /* Get year and check
             */
            cp = strtok(NULL, "-/. ");
            if (cp == NULL)
                goto done;
            yy = atol(cp);
            if (yy > 99 && yy < 1601)
            {
                DrlErrMsg("%s: Year %ld out of range.\n", routine, yy);
    	    goto done;
            }       
    
            if (yy < CROSSOVER_YEAR)  /* Assume its past year 2000 */
                yy += 2000;           /* To get years since 1900 */
            else if (yy < 1061)
                yy += 1900;           /* Else assume we got only 2 digits */     
        }

	IF_FAILED_DONE( DrlDDateMake(mm, dd, yy, aDate));

	status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed reading `%s'.\n", routine, string);
	}
        return(status);    
#endif
}




/*f-------------------------------------------------------------
 * Date routines : print date to string.
 *                                                          
 * <br><br>
 * Prints the date "aDate" in the buffer "string".
 * If "string" is NULL, returns a pointer to a static char string.
 */

char*
DrlDDatePrint(char *string, DDate aDate)
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][64] ;
static	int	tmpIdx=0;
	char	*s ;
	int	mm, dd, yy;

	if (string == NULL) {
		s = tmp[tmpIdx];
		tmpIdx++;
		if (tmpIdx > MAX_IDX-1) tmpIdx=0;
	} else {
		s = string;
	}

	if (DrlDDateSplit(aDate, &mm, &dd, &yy) != SUCCESS)
	{
		sprintf(s, "(bad date)");
	} else {
		sprintf(s, "%02d-%s-%04d",
			dd,
			monthName[(int) mm],
			yy) ;
	}

	return(s)  ;
#undef	MAX_IDX
}

/*f-------------------------------------------------------------
 * Date routines : scan date in YYYYMDD format.
 *                                                          
 * <br><br>
 * Scan the string "string" for a date. Recognized format
 * is YYYYMMDD. Puts the
 * result in "aDate". Returns 0 if success.
 */

int
DrlDDateScanYMD(char *string, DDate *aDate)
{
	int	yyyy, mm, dd ;
	long	ymd;

	if (sscanf(string, "%ld", &ymd) != 1) return(-2);

	yyyy = (int) floor(ymd*1e-4);
	mm = (int) floor((ymd-10000*yyyy)*1e-2);
	dd = ymd - 10000*yyyy - 100*mm;

	if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
	if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

	DrlDDateMake(mm, dd, yyyy, aDate) ;
	return(DrlDDateCheckValid(*aDate)) ;
}

/*f-------------------------------------------------------------
 * Date routines : print date in YYYYMMDD format.
 *                                                          
 * <br><br>
 * Prints the date "aDate" in the buffer "string" in the format YYYMMDD.
 * If "string" is NULL, returns a pointer to a static char string.
 */

char*
DrlDDatePrintYMD(char *string, DDate aDate)
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][64] ;
static	int	tmpIdx=0;
	char	*s ;
	int	mm, dd, yy;

	if (string == NULL) {
		s = tmp[tmpIdx];
		tmpIdx++;
		if (tmpIdx > MAX_IDX-1) tmpIdx=0;
	} else {
		s = string;
	}

	if (DrlDDateSplit(aDate, &mm, &dd, &yy) != SUCCESS)
	{
		sprintf(s, "(bad date)");
	} else {
		sprintf(s, "%04d%02d%02d", yy, mm, dd);
	}

	return(s)  ;
#undef	MAX_IDX
}


/*f-------------------------------------------------------------
 * Date routines : return current system date.
 *                                                          
 * <br><br>
 * Returns the current system date or <i> -1L</i> of failed.
 */

int
DrlDDateToday(DDate *today)
{
static	char	routine[] = "DrlDDateToday";
	int	status = FAILURE;
	time_t		tp;
	struct tm	*hc;

	if (time(&tp) == -1) {
		DrlErrMsg("time: %s\n", routine);
		goto done;
	}
	if ((hc = localtime(&tp)) == NULL) {
		DrlErrMsg("localtime: %s\n", routine);
		goto done;
	}

	IF_FAILED_DONE( DrlDDateMake(
		hc->tm_mon+1,
		hc->tm_mday,
		hc->tm_year+1900,
		today));

	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



/*f-------------------------------------------------------------
 * Date routines : convert date to calendar time.
 *                                                          
 * <br><br>
 * Converts the date "date" to the time (double) format YYYY.XXXX
 * using the "daycount" convention.
 */


int
DrlDDateToTime(DDate date, DDayCount daycount, double *yearfrac)
{
static	char	routine[] = "DrlDDate2Time";
	int	status = FAILURE;

	int	m, d, y;
	double	x;

	IF_FAILED_DONE( DrlDDateSplit(date, &m, &d, &y));

	m -=  1 ;
	x = (double) y;
	x += (double) m * C1_12 ;

	switch (daycount) {
	case DRL_ACT_365F:
	case DRL_ACT_360:
		x += C1_12 * (d - 1) / daysInMonth[m] ;
		break ;
	case DRL_B30_360:
		d = (d >= 30 ? 30 : d) ;
		x += C1_12 * (d - 1) * C1_30 ;
		break ;
	default:
		DrlErrMsg("%s: bad day count.\n", routine);
		goto done;
	}
	*yearfrac = x;

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}

/*f-------------------------------------------------------------
 * Date routines : convert calendar time to date.
 *                                                          
 * <br><br>
 * Computes and returns tha date corresponding to a 
 * time "tt" in format YYYY.XXXX using "daycount" convention.
 */

int
DrlTimeToDDate(double tt, DDayCount daycount, DDate *date)
{
static	char	routine[] = "DrlTime2DDate";
	int	status = FAILURE;

	int	y1, m1, d1 ;

	tt += 1e-4;
	y1 = (int) floor(tt) ;
	m1 = (int) floor(12. * (tt - (double) y1)) ;

	switch (daycount) {
	case DRL_ACT_365F:
	case DRL_ACT_360:
		d1 = (int) floor(daysInMonth[m1] * 12. *
			(tt - (double) y1 - C1_12 * (double) m1)) + 1 ;
		break ;
	case DRL_B30_360:
		d1 = (int) (360. * (tt - (double) y1
			- C1_12 * (double) m1)) + 1 ;
		break ;
	default:
		DrlErrMsg("%s: bad day count.\n", routine);
		return(FAILURE);
	}

	return(DrlDDateMake(m1+1, d1, y1, date)) ;

}


/*f-------------------------------------------------------------
 * Date routines : convert to Julian day.
 *                                                          
 * <br><br>
 * Returns julian day number corresponding to "aDate".
 */

int
DrlDDateToJulday(DDate aDate, long *julday)
{
static	char	routine[] = "DrlDDate2Julday";
	int	status = FAILURE;

	int	ja, jy, jm, mm, id, iyyy ;
	long	jul ;
	long	IGREG  = (15+31L*(10+12L*1582)) ;

	IF_FAILED_DONE( DrlDDateSplit(aDate, &mm, &id, &iyyy));

	if (iyyy == 0) {
		DrlErrMsg("Dat2Julday: there is no year zero.\n");
		return LONG_MAX;
	}
	if (iyyy < 0) ++iyyy;
	if (mm > 2) {
		jy = iyyy;
		jm = mm+1;
	} else {
		jy = iyyy-1;
		jm = mm+13;
	}
	jul = (long) (floor(365.25*jy)+floor(30.6001*jm)+id+1720995);
	if (id+31L*(mm+12L*iyyy) >= IGREG) {
		ja = (int)(0.01*jy) ;
		jul += 2-ja+(int) (0.25*ja);
	}

	*julday = jul;

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}


/*f-------------------------------------------------------------
 * Date routines : convert from Julian day.
 *                                                          
 * <br><br>
 * Returns date corresponding to julian day number "julian".
 */

int
DrlJuldayToDDate(long julday, DDate *date)
{
static	char	routine[] = "DrlJulday2DDate";
	int	status = FAILURE;

	long int	ja, jalpha, jb, jc, jd, je ;
	int		mm, id, iyyy ;
const	long		IGREG  = 2299161 ;

	if (julday >= IGREG) {
		jalpha = (long int) (((float) (julday-1867216)-0.25)/36524.25);
		ja = julday+1+jalpha - (long)(0.25*jalpha);
	} else {
		ja = julday;
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

	if (DrlDDateMake(mm, id, iyyy, date) != SUCCESS) {
		PROGRAM_BUG();
	}

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}



/*f-------------------------------------------------------------
 * Date routines : convert from ALIB date to DR date format
 *                                                          
 * <br><br>
 * Converts a date in DDate type  to a DR date format
 * (encoded as YYYYMMDD in a long). Return -1 if failed.
 */

int
DrlAlibDateToDrDate(long alibDate, DDate *drDate)
{
static	char routine[] = "DrlAlibDate2DDate";
	int	status = FAILURE;

#ifdef	DRL_CLIB
	*drDate = alibDate;
#else
	long	julday;

	/*
	 * Date 1/1/2001 corresponds to ALIB date 146097
	 * and Julian day 2451911.
	 * To obtain Julian date from ALIB date therefore
	 * add the offset 2305814 (=2451911-156097)
	 */
	julday = alibDate + 2305814L;
	IF_FAILED_DONE( DrlJuldayToDDate(julday, drDate));
#endif

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}

/*f-------------------------------------------------------------
 * Date routines : convert from DR date to ALIB date format
 *                                                          
 * <br><br>
 * Converts a date in DR type  to a ALIB date format.
 * Return -1 if failed.
 */

int
DrlDrDateToAlibDate(DDate drDate, long *alibDate)
{
static	char routine[] = "DrlDrDateToAlibDate";
	int	status = FAILURE;

#ifdef	DRL_CLIB
	*alibDate = drDate;
#else
	long	julday;

	/*
	 * Date 1/1/2001 corresponds to ALIB date 146097
	 * and Julian day 2451911.
	 * To obtain ALib date from Julian day therefore
	 * Sunbtract the offset 2305814 (=2451911-156097)
	 */
	IF_FAILED_DONE( DrlDDateToJulday(drDate, &julday));
	*alibDate = julday - 2305814L;
#endif

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}       


/*f-------------------------------------------------------------
 * Date routines : convenience to convert from ALIB date array
 * to DR date format
 *                                                          
 * <br><br>
 * Converts a date in DDate type  to a DR date format
 * (encoded as YYYYMMDD in a long). Return -1 if failed.
 */

int
DrlAlibDateToDrDateArray(
	int numDates,		/* (I) number of dates in array */
	long *alibDates,	/* (I) array of alib dates */
	DDate **drDates)	/* (O) allocated array */
{
static	char routine[] = "DrlAlibDate2DDateArray";
	int	status = FAILURE;
	int	idx;

	*drDates = NULL;
	ASSERT_OR_DONE(numDates >= 0);
	if (numDates == 0) {
		return (SUCCESS);
	}

	*drDates = NEW_ARRAY(DDate, numDates);
	if (*drDates == NULL) goto done;
	for (idx=0; idx<numDates; idx++) {
		IF_FAILED_DONE( DrlAlibDateToDrDate(
			alibDates[idx], &(*drDates)[idx]));
	}

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
	    FREE (*drDates);
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}



/*f-------------------------------------------------------------
 * Date routines : return day of week.
 *                                                          
 * <br><br>
 * Returns the day of the week corresponding to "theDate"
 * (0 is Sunday, 1 is Monday, \dots).
 */

int
DrlWeekDay(DDate date, int *dayno)
{
static	char	routine[] = "DrlWeekDay";
	int	status = FAILURE;

	/*  Calculates day of week (0-6) of YYYYMMDD formatted date  */
static	int	noleap[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
static	int	leap[13] =  {0,31,29,31,30,31,30,31,31,30,31,30,31};

        int	mm, dd, yy, *daysin;

	IF_FAILED_DONE( DrlDDateSplit(date, &mm, &dd, &yy));

        daysin = (DRL_YEAR_IS_LEAP(yy) ? leap : noleap);
        while (mm > 1)
                dd += daysin[--mm];
        if (yy > 0) {
                --yy;
                dd += yy; 
                dd +=  yy/4 - yy/100 + yy/400; /* adjust for leap years */
        }
        *dayno = (dd % 7); /* 0 to 6*/

        /* made it through OK */
        status = SUCCESS;
done:
        if (status != SUCCESS) {
            DrlErrMsg("%s: failed.\n", routine);
        }
        return(status);
}

/*f-------------------------------------------------------------
 * Date routines : print day of week.
 *                                                          
 * <br><br>
 * Prints the day of the week corresponding to number "dayno"
 * (0 is Sunday, 1 is Monday, \dots)
 * in string "string".  Returns "string".
 */

char*
DrlWeekDayPrint(char *string, int dayno)
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][32] ;
static	int	tmpIdx=0;
	char	*s ;

	if (string == NULL) {
		s = tmp[tmpIdx];
		tmpIdx++;
		if (tmpIdx > MAX_IDX-1) tmpIdx=0;
	} else {
		s = string;
	}

	if (dayno < 7) {
		sprintf(s, "%s", shortDayNames[dayno]);
	} else {
		sprintf(s, "Err");
	}

	return(s)  ;
#undef	MAX_IDX
}


/*--------------------------------------------------------------
 * Date routines : print day count convention.
 *                                                          
 * <br><br>
 * Prints the day count convention "daycount" in the buffer "string".
 * If "string" is NULL, returns a pointer to a static char string.
 */

char*
DrlDDayCountPrint(char *string, DDayCount daycount)
{
static	char	*s, buf[16];

	s = (string == NULL ? buf : string);

	switch (daycount) {
		case DRL_ACT_ACT :
			strcpy(s, "ACT/ACT");
			break ;
		case DRL_ACT_365F:
			strcpy(s, "ACT/365F");
			break ;
		case DRL_ACT_360 :
			strcpy(s, "ACT/360");
			break ;
		case DRL_B30_360  :
			strcpy(s, "30/360");
			break ;
		case DRL_B30E_360  :
			strcpy(s, "30E/360");
			break ;
		default:
			sprintf(s, "<UNKNOWN-%ldL>", (long) daycount);
			break ;
	}

	return(s) ;
}



/*--------------------------------------------------------------
 * Date routines : scan day count convention.
 *                                                          
 * <br><br>
 * Scan the string "string" for a day count convention.
 * Recognized formats are such as <i> Act/360, Act/Act</i> etc\dots
 * Puts the result in "daycount". Returns 0 if success.
 */

int
DrlDDayCountScan(char *string, DDayCount *daycount)
{
static	char	routine[] = "DrlDDayCountScan";
	char	*p, buf[256];


	/* advance spaces, make copy  and convert to upper */
	strcpy(buf, string);
	for (p=buf; *p != '\0'; p++) *p = toupper(*p);
	for (p=buf; isspace(*p) && (*p != '\0'); p++);

#undef	COMPARE
#define	COMPARE(cs, ct)	(strncmp((cs), (ct), strlen(ct)) == 0)

	if (COMPARE(p, "ACT/ACT"))
	{
		*daycount = DRL_ACT_ACT;
		return(SUCCESS);
	}
	else if (COMPARE(p, "ACT/365F"))
	{
		*daycount = DRL_ACT_365F;
		return(SUCCESS);
	}
	else if (COMPARE(p, "ACT/360"))
	{
		*daycount = DRL_ACT_360;
		return(SUCCESS);
	}
	else if (COMPARE(p, "30/360"))
	{
		*daycount = DRL_B30_360;
		return(SUCCESS);
	}
	else if (COMPARE(p, "30E/360"))
	{
		*daycount = DRL_B30E_360;
		return(SUCCESS);
	}
	else
	{
		DrlErrMsg("%s: can't read dcc `%s'.\n",
			routine, string);
		return(FAILURE);
	}
#undef	COMPARE
}


/*--------------------------------------------------------------
 * Date routines : get day count convention from num/den.
 *                                                          
 * <br><br>
 * Determines the day count convention "daycount"
 * given a numerator "num" such as 'B' (Bonds) or 'A' (Actual)
 * and a denominator "den" (360, 365, 366).
 * Returns 0 is successful.
 */

int
DrlDDayCountFromNumDen(DDayCount *daycount, char num, int den)
{
static	char	routine[] = "DrlDDayCountFromDenNum";

	switch (num) {
	case 'B':
	case 'b':
	    switch (den) {
	    case 360:
		*daycount = DRL_B30_360;
		return(SUCCESS);
	    default:
		DrlErrMsg("%s: bad denominator %d.\n", routine, den);
		return(FAILURE);
	    }

	case 'A':
	case 'a':
	    switch (den) {
	    case 360:
		*daycount = DRL_ACT_360;
		return(SUCCESS);
	    case 365:
		*daycount = DRL_ACT_365F;
		return(SUCCESS);
	    case 366:
		*daycount = DRL_ACT_ACT;
		return(SUCCESS);
	    default:
		DrlErrMsg("%s: bad denominator %d.\n", routine, den);
		return(FAILURE);
	    }
	default:
	    break;
	}

	DrlErrMsg("%s: bad numerator `%c'.\n", routine, num);
	return(FAILURE);
}


/*--------------------------------------------------------------
 * Date routines : get num/den from day count convention.
 *                                                          
 * <br><br>
 * Given the day count convention "daycount",
 * returns the corresponding numerator "num" such as 'B' (Bonds)
 * or 'A' (Actual) * and denominator "den" (360, 365, 366).
 */

int
DrlDDayCountToNumDen(DDayCount daycount, char *num, int *den)
{
static	char	routine[] = "DrlDDayCountToNumDen";
	switch (daycount) {
	case DRL_ACT_ACT :
		*num = 'A'; *den = 366;
		break ;
	case DRL_ACT_365F :
		*num = 'A'; *den = 365;
		break ;
	case DRL_ACT_360 :
		*num = 'A'; *den = 360;
		break ;
	case DRL_B30_360  :
		*num = 'B'; *den = 360;
		break ;
#if !defined(DRL_CLIB)
	case DRL_B30E_360  :
		*num = 'C'; *den = 360;
		break ;
#endif
	default:
		DrlErrMsg("%s: bad day count conv.\n", routine);
		return(FAILURE);
	}
	return(SUCCESS);
}


#if !defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)

/*f-------------------------------------------------------------
 * Date routines : add number of days to a date.
 *                                                          
 * <br><br>
 * Adds an interval to a date.
 */

int
DrlDDateAddDays(
	DDate startDate,		/* (I) start date */
	int numDays,			/* (I) number of days */
	DDate *endDate)			/* (O) final date */
{
	*endDate = Nxtday(startDate, (long) numDays);
	return(SUCCESS);
}


/*f-------------------------------------------------------------
 * Date routines : add interval to a date.
 *                                                          
 * <br><br>
 * Adds an interval to a date.
 */

int
DrlDDateFwdAny(
	DDate startDate,		/* (I) start date */
	DInterval *interval,	/* (I) interval */
	DDate *endDate)			/* (O) final date */
{
static	char	routine[]="DrlDDateFwdAny";
	int	status = FAILURE;

	int	numPeriods;
	char	upperPrdTyp = (char) toupper(interval->prd_typ);

	int	mm, dd, yy;
	int	endOfMonth;

    
	switch ((char) toupper(interval->prd_typ))
	{
	case 'M':
	case 'A':
	case 'Y':
	case 'S':
	case 'Q':
		if (upperPrdTyp == 'M')  
			numPeriods = interval->prd;
		else if (upperPrdTyp == 'A' || interval->prd_typ == 'Y')
			numPeriods = interval->prd * 12;
		else if (upperPrdTyp == 'S')
			numPeriods = interval->prd * 6;
		else  
			numPeriods = interval->prd * 3;


		*endDate  = Nxtmth(startDate, (long) numPeriods, FALSE);

		break;

        
	case 'D':                         /* DAYly increments */
		numPeriods = interval->prd;
		*endDate = Nxtday(startDate, (long) numPeriods);

		break;

	case 'W':
		numPeriods = interval->prd * 7;
		*endDate = Nxtday(startDate, (long) numPeriods);

		break;

	case 'E':                         /* End of month forced */
		numPeriods = interval->prd;
		*endDate  = Nxtmth(startDate, (long) numPeriods, TRUE);

		break ;

	case 'I':                         /*  IMM dates (quarterly) */
		numPeriods = interval->prd;
		*endDate =Nxtimm(startDate, (long) numPeriods);

		break;
        
	default:
		DrlErrMsg("%s:  Period type %c unknown.\n",
			routine, interval->prd_typ);
		goto done;
	}


	status = SUCCESS;
done:
	if (status != SUCCESS) 
		DrlErrMsg("%s: failed.\n", routine);

	return(status);
}


/*f-------------------------------------------------------------
 * Date routines : computes and returns actual number of days bwteen two dates
 *                                                          
 * <br><br>
 */

int
DrlDDateDaysDiff(
	DDate startDate,	/* (I) start date */
	DDate endDate)		/* (I) end date   */ 
{
	return (int) Daysact(startDate, endDate);
}


/*f-------------------------------------------------------------
 * Date routines : compute day count fraction given a method
 *                                                          
 * <br><br>
 */

int
DrlDayCountFract(
	DDate date1,		/* (I) start date */
	DDate date2,		/* (I) end date   */ 
	DDayCount method,	/* (I) method */
	double *fract)		/* (O) day count fraction */
{
static	char	routine[] = "DrlDayCountFract";
	int	status = FAILURE;       /* Until successful */


	switch (method) {
	case DRL_ACT_365F:
		*fract = (double) DrlDDateDaysDiff(date1, date2) / 365e0;
		break;

	case DRL_ACT_360:
		*fract = (double) DrlDDateDaysDiff(date1, date2) / 360e0;
		break;
	
	case DRL_B30_360:
		*fract = (double) Days360(date1, date2) / 360e0;
		break;

	default:
		DrlErrMsg("%s: bad dcc %c.\n", routine, method);
		goto done;
	}

	status = SUCCESS;
done:
	if (status != SUCCESS) 
		DrlErrMsg("%s: failed.\n", routine);

	return(status);
}

#endif	/* !defined(DRL_CLIB) && !defined(DRL_TDATE_MDY) */



/*==============================================================*/
/*								*/
/*	SRM3 ROUTINES						*/
/*								*/
/*==============================================================*/

#define HOLIDAY_FILE "HOLIDAY"
#define OKAY 0
#define ERR -1


static	long noleap[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
static	long leap[13] =  {0,31,29,31,30,31,30,31,31,30,31,30,31};
static	long cumdays[12] = {0,31,59,90,120,151,181,212,243,273,304,334};


static	void Dsplit(long date_i, long *mm_o, long *dd_o, long *yy_o)
{
    /*  Returns month, day, year from an integer in                */
    /*  the format YYY(Y)MMDD as in 840701 or 20011115 or 1120104  */

    *yy_o = date_i/10000;
    *mm_o = (date_i - *yy_o*10000)/100;
    *dd_o = date_i - *yy_o*10000 - *mm_o*100;
    return;
}

static	int Isleap(long year)
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


static	int Isimm(long    date)
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


static	long   ThirdWed(long mth, long year)
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


static	int Dateok(long date)
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


static	long Y2toy4(long year_i)
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


static	long Datepack(long mm_i, long dd_i, long yy_i)
{
    /*  Packs mm_i,dd_i,yy_i into YYYYMMDD format  */
    /*  calls Y2toy4 on yy_i before packing        */

    long y4,packed_o;

    y4 = Y2toy4(yy_i);
    packed_o = y4*10000 + mm_i*100 + dd_i;
    return (packed_o);
}


static	long Y4toy2(long year_i)
{
    /*  Converts a 4 digit year to a 2 digit year;  */
    /*  years after 99 are 00, 01, etc.             */

    long y2_o;

    y2_o = year_i%100;
    return (y2_o);
}


static	void Y2date_str(long date, char *string)
{
    /*  Convert YYYYMMDD to MM/DD/YY  */

    long mm,dd,yy;

    Dsplit(date,&mm,&dd,&yy);
    yy = Y4toy2(yy);
    sprintf(string,"%02ld/%02ld/%02ld",mm,dd,yy);

    return;
}


static	long eval_date(char *datest)
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


static	long eval_date2(char *datest)
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
        month = 0;	                        /* Error */
                
    date = year * 10000 + month * 100 + day;

    return(date);
}


static	long Dayofwk(long date)
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
static	long Days360(long date1_i, long date2_i)
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


static	long Months360(long date1_i, long date2_i)
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


static	long Daysact(long date1_i, long date2_i)
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


static	long Nxtday(long date, long days)
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


static	long Nxtmth(long date, long mths, long eom)
{
    /*  Returns a date in YYYYMMDD format based on moving        */
    /*  forward or backwards (mths) from (date).                 */
    /*  if eom > 0 the returned does not spill over to the next  */
    /*  month. That is moving 1m from the 31rst of march gives   */
    /*  the 30th of april. If eom=0 it does spill over: 1m from  */
        /*  the 31rst of march is the 1rst of may.	         */
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




static	long  Nxtimm(long  date, long nbPeriods)
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



static	long Nxtwkday(long date, long advance)
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





static	int
strToMonth(char *cp, int *monthN)
{
        cp[0] = tolower(cp[0]);
        cp[1] = tolower(cp[1]);
        cp[2] = tolower(cp[2]);
        if (strcmp(cp,"jan") == 0)
        {
            *monthN = 1;
        }
        else if (strcmp(cp,"feb") == 0)
        {
            *monthN = 2;
        }
        else if (strcmp(cp,"mar") == 0)
        {
            *monthN = 3;
        }
        else if (strcmp(cp,"apr") == 0)
        {
            *monthN = 4;
        }
        else if (strcmp(cp,"may") == 0)
        {
            *monthN = 5;
        }
        else if (strcmp(cp,"jun") == 0)
        {
            *monthN = 6;
        }
        else if (strcmp(cp,"jul") == 0)
        {
            *monthN = 7;
        }
        else if (strcmp(cp,"aug") == 0)
        {
            *monthN = 8;
        }
        else if (strcmp(cp,"sep") == 0)
        {
            *monthN = 9;
        }
        else if (strcmp(cp,"oct") == 0)
        {
            *monthN = 10;
        }
        else if (strcmp(cp,"nov") == 0)
        {
            *monthN = 11;
        }
        else if (strcmp(cp,"dec") == 0)
        {
            *monthN = 12;
        }
        else
        {
            return (FAILURE);
        }
        return (SUCCESS);
}


