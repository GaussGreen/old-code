/************************************************************************
 * Module:	DRL
 * Submodule:	TIME
 * File:	
 * Function:	TDate and Time Routines
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#ifdef CLIB
# include "dateconv.h"
#endif

#include "drlstr.h"

#include "drlstr.h"		/* DrlStrRemoveChars() */

#include "drltime.h"

#define MIN_YEAR	1602
#define	MAX_YEAR	2070

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

#ifdef	_SKIP
static	int gDateDaysInMonth[2][13] =
{
	{31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
	{31, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

/* Cumulative to Dec 31 INclusive */
static	int gDateDaysToYearEnd[2][13] =	
{	
	{31, 365, 334, 306, 275, 245, 214, 184, 153, 122, 92, 61, 31},
	{31, 366, 335, 306, 275, 245, 214, 184, 153, 122, 92, 61, 31}
};

/* Cumulative from Dec 31 EXclusive */
static	int gDateDaysFromYearStart[2][13] =
{ 
	{334, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
	{334, 0, 31, 60, 91, 121, 152, 182, 213, 245, 274, 305, 335}
};

static	char	gDateLongDayNames[12][8]
		= {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday"
			, "Friday", "Saturday", "Sunday"};
static	char	gDateShortmonthNames[12][13]
		= {"Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"
			, "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
static	char	gDateLongmonthNames[12][13]
		= {"December", "January", "February", "March", "April"
			     , "May", "June", "July", "August"
			     , "September", "October", "November", "December"};
#endif


/*f-------------------------------------------------------------
 * Date routines : check date valid.
 *                                                          
 * <br><br>
 * Returns SUCCESS (0) if "aDate" is a valid date between
 * 1900 and 2070, and FAILURE (1) otherwise.
 */

DLL_EXPORT(int)
DrlTDateCheckValid(TDate aDate)
{
static	char	routine[] = "DrlTDateCheckValid";
	int		status = FAILURE;
	int		yy, mm, dd;
#ifdef  CLIB
	TMonthDayYear	mdy;

	if (GtoDateToMDY(aDate, &mdy) != SUCCESS) goto done;
	dd = (int) mdy.day;
	mm = (int) mdy.month;
	yy = (int) mdy.year;
#else
	dd = (int) DrlDAY(aDate);
	mm = (int) DrlMONTH(aDate);
	yy = (int) DrlYEAR(aDate);
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
	    GtoErrMsg("%s: invalid date %ld (m=%d,d=%d,y=%d)\n",
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

DLL_EXPORT(int)
DrlTDateSplit(TDate aDate, int *mm, int *dd, int *yy)
{
#ifdef	CLIB
	int	errCode =0;
	TMonthDayYear	mdy;

	errCode = GtoDateToMDY(aDate, &mdy);
	*dd = (int) mdy.day;
	*mm = (int) mdy.month;
	*yy = (int) mdy.year;
	return(errCode);
#else
	*dd = (int) DrlDAY(aDate);
	*mm = (int) DrlMONTH(aDate);
	*yy = (int) DrlYEAR(aDate);
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

DLL_EXPORT(TDate)
DrlTDateMake(int m, int d, int y)
{
#ifdef	CLIB
	TMonthDayYear	mdy;
	TDate		theDate;

	mdy.month = (long) m;
	mdy.day = (long) d;
	mdy.year = (long) y;

	if (GtoMDYToDate(&mdy, &theDate) != 0) {
		return DRL_DATE_NULL;
	} else {
		return theDate;
	}
#else
	return DRL_DATE_NULL;
#endif

}

/*f-------------------------------------------------------------
 * Date routines : get day from date.
 *                                                          
 * <br><br>
 * Returns the day corresponding to date "dt".
 */

DLL_EXPORT(int)
DrlDAY(TDate dt)
{
	TMonthDayYear	mdy;

	GtoDateToMDY(dt, &mdy);
	return (int) mdy.day;
}

/*f-------------------------------------------------------------
 * Date routines : get month from date.
 *                                                          
 * <br><br>
 * Returns the month corresponding to date "dt".
 */

DLL_EXPORT(int)
DrlMONTH(TDate dt)
{
	TMonthDayYear	mdy;

	GtoDateToMDY(dt, &mdy);
	return (int) mdy.month;
}



/*f-------------------------------------------------------------
 * Date routines : get year from date.
 *                                                          
 * <br><br>
 * Returns the year corresponding to date "dt".
 */

DLL_EXPORT(int)
DrlYEAR(TDate dt)
{
	TMonthDayYear	mdy;

	GtoDateToMDY(dt, &mdy);
	return (int) mdy.year;
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

DLL_EXPORT(int)
DrlTDateScan(char *string, TDate *aDate)
{
static	char	routine[] = "DrlTDateScan";
	char	dateCopy[32];
	int	yyyy, mm, dd;


	/* special case of yyyy.mm.dd not handled by GTO */
	if (strchr(string, '.') != NULL) {
		if (sscanf(string, "%d.%d.%d", &yyyy, &mm, &dd) != 3) {
			GtoErrMsg("%s: can't read date `%s'.\n", 
				  routine, string);
			return(FAILURE);
		}
		sprintf(dateCopy, "%04d%02d%02d", yyyy, mm, dd);
	} else {
		strncpy(dateCopy, string, sizeof(dateCopy)-1);
	}

	/* Use Gto routine */
	return GtoStringToDate(dateCopy, aDate);


#ifdef	_SKIP	/* Do not use */
static	char	routine[] = "DrlTDateScan";

	int	yyyy, mm, dd ;
	char	s0[255], *s1, *s2 ;
static	char	toks[] = "-/ ,";

#define	FMT_EXTENDED
#ifdef	FMT_EXTENDED

	if ((s1 = DrlStrToken(string, toks, s0, &s2)) == NULL) return(__LINE__);

	if (sscanf(s1, "%d", &mm) == 1) {
		/*
		 * scan for `dd-Mon-yyyy' or `mm-dd-yyyy'
		 */
		if ((s1 = DrlStrToken(NULL, toks, s0, &s2)) == NULL)
			return(__LINE__);
		if (sscanf(s1, "%d", &dd) == 1) {
			/*
			 * format `mm-dd-yyyy'
			 */
			if ((s1 = DrlStrToken(NULL, toks, s0, &s2)) == NULL)
				return(__LINE__);
			if (sscanf(s1, "%d", &yyyy) != 1)
				return(__LINE__);
		} else {
			/*
			 * format `dd-Mon-yyyy'
			 */
			dd = mm ;
			mm = 0;
			while ((++mm <= 12)
				&& (strncmp(s1, monthName[mm], 3) != 0));
			if (mm >= 13) return(__LINE__);

			if ((s1 = DrlStrToken(NULL, toks, s0, &s2)) == NULL)
				return(__LINE__);
			if (sscanf(s1, "%d", &yyyy) != 1)
				return(__LINE__);
		}
	} else {
		/*
		 * format `Mon dd, yyy'
		 */
		mm = 0;
		while ((++mm <= 12) && (strncmp(s1, monthName[mm], 3) != 0));
		if (mm >= 13) return(__LINE__);

		if ((s1 = DrlStrToken(NULL, toks, s0, &s2)) == NULL)
			return(__LINE__);
		if (sscanf(s1, "%d", &dd) != 1)
			return(__LINE__);
		if ((s1 = DrlStrToken(NULL, toks, s0, &s2)) == NULL)
			return(__LINE__);
		if (sscanf(s1, "%d", &yyyy) != 1)
			return(__LINE__);
	}

#else	/* FMT_EXTENDED */

	if (sscanf(string, "%d/%d/%d", &mm, &dd, &yyyy) == 3) {
		if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
		if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

		*aDate = DrlTDateMake(mm, dd, yyyy) ;
		return(DrlTDateCheckValid(*aDate)) ;
	}
	if (sscanf(string, "%d-%s-%d", &dd, s0, &yyyy) == 3) {
		if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
 		if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

		mm = 0;
		while ((++mm <= 12) && (strncmp(s0, monthName[mm], 3) != 3));
 		*aDate = DrlTDateMake(mm, dd, yyyy) ;
		return(DrlTDateCheckValid(*aDate)) ;
	}
#endif	/* FMT_EXTENDED */


	if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
	if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

 	*aDate = DrlTDateMake(mm, dd, yyyy) ;
	return(DrlTDateCheckValid(*aDate)) ;
#endif
}


/*f-------------------------------------------------------------
 * Date routines : print date to string.
 *                                                          
 * <br><br>
 * Prints the date "aDate" in the buffer "string".
 * If "string" is NULL, returns a pointer to a static char string.
 */

DLL_EXPORT(char*)
DrlTDatePrint(char *string, TDate aDate)
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][64] ;
static	int	tmpIdx=0;
	char	*s ;

	if (string == NULL) {
		s = tmp[tmpIdx];
		tmpIdx++;
		if (tmpIdx > MAX_IDX-1) tmpIdx=0;
	} else {
		s = string;
	}

	/*sprintf(s, "%02d/%02d/%04d",
		DrlMONTH(aDate), DrlDAY(aDate), DrlYEAR(aDate)) ;*/
	sprintf(s, "%02d-%s-%04d",
		DrlDAY(aDate),
		monthName[(int) DrlMONTH(aDate)],
		DrlYEAR(aDate)) ;

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

DLL_EXPORT(int)
DrlTDateScanYMD(char *string, TDate *aDate)
{
	int	yyyy, mm, dd ;
	long	ymd;

	if (sscanf(string, "%ld", &ymd) != 1) return(-2);

	yyyy = (int) floor(ymd*1e-4);
	mm = (int) floor((ymd-10000*yyyy)*1e-2);
	dd = ymd - 10000*yyyy - 100*mm;

	if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
	if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

	*aDate = DrlTDateMake(mm, dd, yyyy) ;
	return(DrlTDateCheckValid(*aDate)) ;
}

/*f-------------------------------------------------------------
 * Date routines : print date in YYYYMMDD format.
 *                                                          
 * <br><br>
 * Prints the date "aDate" in the buffer "string" in the format YYYMMDD.
 * If "string" is NULL, returns a pointer to a static char string.
 */

DLL_EXPORT(char*)
DrlTDatePrintYMD(char *string, TDate aDate)
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][64] ;
static	int	tmpIdx=0;
	char	*s ;

	if (string == NULL) {
		s = tmp[tmpIdx];
		tmpIdx++;
		if (tmpIdx > MAX_IDX-1) tmpIdx=0;
	} else {
		s = string;
	}

	sprintf(s, "%04d%02d%02d", DrlYEAR(aDate), DrlMONTH(aDate), DrlDAY(aDate));

	return(s)  ;
#undef	MAX_IDX
}


/*f-------------------------------------------------------------
 * Date routines : return current system date.
 *                                                          
 * <br><br>
 * Returns the current system date or <i> -1L</i> of failed.
 */

DLL_EXPORT(TDate)
DrlTDateToday(void)
{
static	char	routine[] = "DrlTDateToday";
	int	status = FAILURE;
	time_t		tp;
	struct tm	*hc;
	TDate	today = -1L;

	if (time(&tp) == -1) {
		GtoErrMsg("time: %s\n", routine);
		goto done;
	}
	if ((hc = localtime(&tp)) == NULL) {
		GtoErrMsg("localtime: %s\n", routine);
		goto done;
	}

	today = DrlTDateMake(hc->tm_mon+1, hc->tm_mday,
		hc->tm_year+1900);

	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
		return(-1L);
	}
	return(today);
}



/*f-------------------------------------------------------------
 * Date routines : convert date to calendar time.
 *                                                          
 * <br><br>
 * Converts the date "date" to the time (double) format YYYY.XXXX
 * using the "daycount" convention.
 */


DLL_EXPORT(double)
DrlTDate2Time(TDate date, TDayCount daycount)
{
static	double	x ;
static	int	m, d ;
	m = DrlMONTH(date) - 1 ;
	x = (double) DrlYEAR(date) ;

	x += (double) m * C1_12 ;

	switch (daycount) {
	case GTO_ACT_365:
	case GTO_ACT_365F:
	case GTO_ACT_360:
		x += C1_12 * (DrlDAY(date) - 1) / daysInMonth[m] ;
		break ;
	case GTO_B30_360:
		d = ((d = DrlDAY(date)) >= 30 ? 30 : d) ;
		x += C1_12 * (d - 1) * C1_30 ;
		break ;
	default:
		GtoErrMsg("Date2Time: not implemented\n");
		x = DBL_MAX;
	}


	return(x) ;
}

/*f-------------------------------------------------------------
 * Date routines : convert calendar time to date.
 *                                                          
 * <br><br>
 * Computes and returns tha date corresponding to a 
 * time "tt" in format YYYY.XXXX using "daycount" convention.
 */

DLL_EXPORT(TDate)
DrlTime2TDate(double tt, TDayCount daycount)
{
	int	y1, m1, d1 ;

	tt += 1e-4;
	y1 = (int) floor(tt) ;
	m1 = (int) floor(12. * (tt - (double) y1)) ;

	switch (daycount) {
	case GTO_ACT_365:
	case GTO_ACT_365F:
	case GTO_ACT_360:
		d1 = (int) floor(daysInMonth[m1] * 12. *
			(tt - (double) y1 - C1_12 * (double) m1)) + 1 ;
		break ;
	case GTO_B30_360:
		d1 = (int) (360. * (tt - (double) y1
			- C1_12 * (double) m1)) + 1 ;
		break ;
	default:
		GtoErrMsg("Time2Date: not implemented\n");
		return(DRL_DATE_NULL);
	}

	return(DrlTDateMake(m1+1, d1, y1)) ;

}


/*f-------------------------------------------------------------
 * Date routines : convert to Lotus date.
 *                                                          
 * <br><br>
 * Converts a date "aDate" to Lotus 123 integer format.
 */

DLL_EXPORT(long)
DrlTDate2Lotus(TDate aDate)
{
	long	n123, n123ref;

	/* Lotus 123 date 367 corresponds to 1/1/1901 */
	n123ref = DrlTDateMake(1,1,1901);
	n123    = aDate;
	return((long) (n123 - n123ref + 367L));
}


/*f-------------------------------------------------------------
 * Date routines : convert from Lotus date.
 *                                                          
 * <br><br>
 * Converts a Lotus 123 date in integer format to a 
 * <i> TDate</i> format.
 */

DLL_EXPORT(TDate)
DrlLotus2TDate(long n123)
{
	long	n123ref;
	TDate	aDate;

	/* Lotus 123 date 367 corresponds to 1/1/1901 */
	n123ref = DrlTDateMake(1,1,1901);

	if (n123 < 367L) {
		GtoErrMsg("DrlLotus2TDate: invalid Lotus date %ld\n", n123);
		return(DRL_DATE_NULL);
	}
	aDate = (TDate) n123 + n123ref - 367L;
	return(aDate);
}



/*f-------------------------------------------------------------
 * Date routines : convert to Julian day.
 *                                                          
 * <br><br>
 * Returns julian day number corresponding to "aDate".
 */

DLL_EXPORT(long)
DrlTDate2Julday(TDate aDate)
{
	int	ja, jy, jm, mm, id, iyyy ;
	long	jul ;
	long	IGREG  = (15+31L*(10+12L*1582)) ;

	id = DrlDAY(aDate) ;
	mm = DrlMONTH(aDate) ;
	iyyy = DrlYEAR(aDate) ;

	if (iyyy == 0) {
		GtoErrMsg("Dat2Julday: there is no year zero.\n");
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
	return(jul) ;
}


/*f-------------------------------------------------------------
 * Date routines : convert from Julian day.
 *                                                          
 * <br><br>
 * Returns date corresponding to julian day number "julian".
 */

DLL_EXPORT(TDate)
DrlJulday2TDate(long julday)
{
	long int	ja, jalpha, jb, jc, jd, je ;
	int		mm, id, iyyy ;
const	long		IGREG  = 2299161 ;
	TDate		aDate ;

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

	aDate = DrlTDateMake(mm, id, iyyy);
	return(aDate) ;
}



/*f-------------------------------------------------------------
 * Date routines : convert from YYYYMMDD long format.
 *                                                          
 * <br><br>
 * Converts a DR date type (encoded as YYYYMMDD in a long)
 * to a TDate. Returns -1 if failed.
 */

DLL_EXPORT(TDate)
DrlDrDate2TDate(long drDate)
{
static	char routine[] = "DrlDrDate2TDate";
	TDate		tDate;
	TMonthDayYear mdy;

	mdy.year  =  drDate/10000; /* Truncation forced */
	mdy.month = (drDate - 10000*mdy.year)/100;
	mdy.day   = (drDate - 10000*mdy.year - 100*mdy.month);

	if (GtoMDYToDate(&mdy, &tDate) != SUCCESS) {
		GtoErrMsg("%s: failed.\n",routine);
		return(-1L);
	}
	return(tDate);
}       


/*f-------------------------------------------------------------
 * Date routines : convert to YYYYMMDD long format.
 *                                                          
 * <br><br>
 * Converts a date in TDate type  to a DR date format
 * (encoded as YYYYMMDD in a long). Return -1 if failed.
 */

DLL_EXPORT(long)
DrlTDate2DrDate(TDate tDate)
{
static	char routine[] = "DrlTDate2DrDate";
	long		drDate;
	TMonthDayYear	mdy;

	if (GtoDateToMDY(tDate, &mdy) != SUCCESS) {
		GtoErrMsg("%s: failed.\n",routine);
		return(-1L);
	}

	drDate  = (long) mdy.year * 10000L;
	drDate += (long) mdy.month * 100L;
	drDate += (long) mdy.day;

	return(drDate);
}       


/*f-------------------------------------------------------------
 * Date routines : return day of week.
 *                                                          
 * <br><br>
 * Returns the day of the week corresponding to "theDate"
 * (0 is Sunday, 1 is Monday, \dots).
 */

DLL_EXPORT(int)
DrlWeekDay(TDate date)
{
#define	_USE_LIONNELS_ROUTINE_
#ifdef	_USE_LIONNELS_ROUTINE_

	/*  Calculates day of week (0-6) of YYYYMMDD formatted date  */
static	int	noleap[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
static	int	leap[13] =  {0,31,29,31,30,31,30,31,31,30,31,30,31};

        int	mm, dd, yy, *daysin;
	dd = DrlDAY(date);
	mm = DrlMONTH(date);
	yy = DrlYEAR(date);

        daysin = (DRL_YEAR_IS_LEAP(yy) ? leap : noleap);
        while (mm > 1)
                dd += daysin[--mm];
        if (yy > 0) {
                --yy;
                dd += yy; 
                dd +=  yy/4 - yy/100 + yy/400; /* adjust for leap years */
        }
        return (dd % 7); /* 0 to 6*/

#else
	long	jul ;
	jul = DrlTDate2Julday(date);
	return (int) ((jul+1) % 7);
#endif
}

/*f-------------------------------------------------------------
 * Date routines : print day of week.
 *                                                          
 * <br><br>
 * Prints the day of the week corresponding to number "wDay"
 * (0 is Sunday, 1 is Monday, \dots)
 * in string "string".  Returns "string".
 */

DLL_EXPORT(char*)
DrlWeekDayPrint(char *string, int wDay)
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

	if (wDay < 7) {
		sprintf(s, "%s", shortDayNames[wDay]);
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

DLL_EXPORT(char*)
DrlTDayCountPrint(char *string, TDayCount daycount)
{
static	char	*s, buf[16];

	s = (string == NULL ? buf : string);

	switch (daycount) {
		case GTO_ACT_365 :
			strcpy(s, "ACT/ACT");
			break ;
		case GTO_ACT_365F:
			strcpy(s, "ACT/365F");
			break ;
		case GTO_ACT_360 :
			strcpy(s, "ACT/360");
			break ;
		case GTO_B30_360  :
			strcpy(s, "30/360");
			break ;
		case GTO_B30E_360  :
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

DLL_EXPORT(int)
DrlTDayCountScan(char *string, TDayCount *daycount)
{
	char	*p, buf[256];


	/* advance spaces, make copy  and convert to upper */
	strcpy(buf, string);
	for (p=buf; *p != '\0'; p++) *p = toupper(*p);
	for (p=buf; isspace(*p) && (*p != '\0'); p++);

#undef	COMPARE
#define	COMPARE(cs, ct)	(strncmp((cs), (ct), strlen(ct)) == 0)

	if (COMPARE(p, "ACT/ACT"))
	{
		*daycount = GTO_ACT_365;
		return(SUCCESS);
	}
	else if (COMPARE(p, "ACT/365F"))
	{
		*daycount = GTO_ACT_365F;
		return(SUCCESS);
	}
	else if (COMPARE(p, "ACT/360"))
	{
		*daycount = GTO_ACT_360;
		return(SUCCESS);
	}
	else if (COMPARE(p, "30/360"))
	{
		*daycount = GTO_B30_360;
		return(SUCCESS);
	}
	else if (COMPARE(p, "30E/360"))
	{
		*daycount = GTO_B30E_360;
		return(SUCCESS);
	}
	else
	{
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

DLL_EXPORT(int)
DrlTDayCountFromNumDen(TDayCount *daycount, char num, int den)
{
static	char	routine[] = "DrlTDayCountFromDenNum";

	switch (num) {
	case 'B':
	case 'b':
	    switch (den) {
	    case 360:
		*daycount = GTO_B30_360;
		return(SUCCESS);
	    default:
		GtoErrMsg("%s: bad denominator %d.\n", routine, den);
		return(FAILURE);
	    }

	case 'A':
	case 'a':
	    switch (den) {
	    case 360:
		*daycount = GTO_ACT_360;
		return(SUCCESS);
	    case 365:
		*daycount = GTO_ACT_365F;
		return(SUCCESS);
	    case 366:
		*daycount = GTO_ACT_365;
		return(SUCCESS);
	    default:
		GtoErrMsg("%s: bad denominator %d.\n", routine, den);
		return(FAILURE);
	    }
	default:
	    break;
	}

	GtoErrMsg("%s: bad numerator `%c'.\n", routine, num);
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

DLL_EXPORT(int)
DrlTDayCountToNumDen(TDayCount daycount, char *num, int *den)
{
static	char	routine[] = "DrlTDayCountToNumDen";
	switch (daycount) {
	case GTO_ACT_365 :
		*num = 'A'; *den = 366;
		break ;
	case GTO_ACT_365F :
		*num = 'A'; *den = 365;
		break ;
	case GTO_ACT_360 :
		*num = 'A'; *den = 360;
		break ;
	case GTO_B30_360  :
		*num = 'B'; *den = 360;
		break ;
	case GTO_B30E_360  :
		*num = 'C'; *den = 360;
		break ;
	default:
		GtoErrMsg("%s: bad day count conv.\n", routine);
		return(FAILURE);
	}
	return(SUCCESS);
}



