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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include "date_sup.h"		/* C analytics */

#include "drltime.h"		/* prototype consistency */


static	TDate	drlAddTDateInterval(TDate aDate, TDateInterval tlaps);
static	double	drlTDateIntervalDiff(TDate aDate, TDateInterval tlapse);


/*f-------------------------------------------------------------
 * Date routines : create date interval.
 *                                                          
 * <br><br>
 * Creates a number "nDays" of days, a number "nMonths" of months
 * and a number "nYears" of years in a date interval.
 * Returns the date interval created.
 */

DLL_EXPORT(TDateInterval)
DrlTDateIntervalFromDMY(int nDays, int nMonths, int nYears)
{
	TDateInterval	tlaps;

	if (nYears != 0) {
		if (nMonths != 0) {
			GtoMakeDateInterval(12*nYears+nMonths, 'M', &tlaps);
			return (tlaps);
		} else {
			GtoMakeDateInterval(nYears, 'A', &tlaps);
			return (tlaps);
		}
	}
	if (nMonths != 0) {
		GtoMakeDateInterval(nMonths, 'M', &tlaps);
		return (tlaps);
	}

	GtoMakeDateInterval(nDays, 'D', &tlaps);
	return (tlaps);

}

#ifdef	_SKIP
/*f-------------------------------------------------------------
 * Converts a (possibly fractional) number of years "tt"
 * in a date interval by rounding it to the closest number of months.
 * For example, 0.25 returns 3M, 1.4 is 18M.
 */

DLL_EXPORT(TDateInterval)
DrlTDateIntervalFromDouble(double tt)
{
	TDateInterval	tlaps;
	int	ny, nm, nd;
#define	HALFDAY	(1.3698e-3)
	/*
	 *
	 */
	ny = (int) floor(tt + HALFDAY);
	nm = (int) floor(12e0 * (tt - (double) ny + HALFDAY));
	nd = (int) (365e0 * (tt - (double) ny - (double) nm / 12e0
		+ HALFDAY));

	tlaps = DrlTDateIntervalFromDMY(nd, nm, ny);
	return(tlaps);
#undef	HALFDAY
}
#endif

DLL_EXPORT(TDateInterval)
DrlTDateIntervalNil(void)
{
	TDateInterval	tlaps;
	tlaps.prd_typ = '\0';
	return(tlaps);
}


/*--------------------------------------------------------------
 * Allocates an array of date intervals of size "size" and returns
 * a pointer to it. Returns NULL if failure.
 */

DLL_EXPORT(TDateInterval*)
DrlTDateIntervalArrayAlloc(int size)
{
	return (TDateInterval*) MALLOC(size*sizeof(TDateInterval));
}

/*--------------------------------------------------------------
 * Frees an array previousely allocated by <i> DrlTDateIntervalArrayAlloc</i>.
 */

DLL_EXPORT(int)
DrlTDateIntervalArrayFree(TDateInterval* ptr)
{
	FREE((void*)ptr);
	return(0);
}


/*f-------------------------------------------------------------
 * Date routines : compare date intervals.
 *                                                          
 * <br><br>
 * Comparison test on date intervals. A reference date needs to
 * be provided for ``date-dependent'' intervals (e.g. IMM).
 */

DLL_EXPORT(int)
DrlTDateIntervalIsEqual(TDate refDate, TDateInterval t1, TDateInterval t2)
{
	return (fabs(drlAddTDateInterval(refDate, t1) - drlAddTDateInterval(refDate, t2)) <= 1e-6);
}

/*f-------------------------------------------------------------
 * Date routines : compare date intervals.
 *                                                          
 * <br><br>
 * Comparison test on date intervals. A reference date needs to
 * be provided for ``date-dependent'' intervals (e.g. IMM).
 */

DLL_EXPORT(int)
DrlTDateIntervalIsLowerOrEqual(TDate refDate, TDateInterval t1, TDateInterval t2)
{
	return (drlAddTDateInterval(refDate, t1) <= drlAddTDateInterval(refDate, t2));
}

/*f-------------------------------------------------------------
 * Date routines : compare date intervals.
 *                                                          
 * <br><br>
 * Comparison test on date intervals. A reference date needs to
 * be provided for ``date-dependent'' intervals (e.g. IMM).
 */

DLL_EXPORT(int)
DrlTDateIntervalIsGreaterOrEqual(TDate refDate, TDateInterval t1, TDateInterval t2)
{
	return (drlAddTDateInterval(refDate, t1) >= drlAddTDateInterval(refDate, t2));
}

/*f-------------------------------------------------------------
 * Date routines : compare date intervals.
 *                                                          
 * <br><br>
 * Comparison test on date intervals. A reference date needs to
 * be provided for ``date-dependent'' intervals (e.g. IMM).
 */

DLL_EXPORT(int)
DrlTDateIntervalIsLower(TDate refDate, TDateInterval t1, TDateInterval t2)
{
	return (drlAddTDateInterval(refDate, t1) < drlAddTDateInterval(refDate, t2));
}

/*f-------------------------------------------------------------
 * Date routines : compare date intervals.
 *                                                          
 * <br><br>
 * Comparison test on date intervals. A reference date needs to
 * be provided for ``date-dependent'' intervals (e.g. IMM).
 */

DLL_EXPORT(int)
DrlTDateIntervalIsGreater(TDate refDate, TDateInterval t1, TDateInterval t2)
{
	return (drlAddTDateInterval(refDate, t1) > drlAddTDateInterval(refDate, t2));
}


/*--------------------------------------------------------------
 * Computes and returns the date obtained by adding
 * the time interval "tlaps" to the date "aDate".
 */

static	TDate
drlAddTDateInterval(TDate aDate, TDateInterval tlaps)
{
	TDate	toTDate;
	GtoDtFwdAny(aDate, &tlaps, &toTDate);
	return toTDate;
}

/*--------------------------------------------------------------
 * Returns the time length in years (Act/365) between
 * a date "aDate" and the date obtained by adding
 * the date interval "tlapse" to it.
 */

static	double
drlTDateIntervalDiff(TDate aDate, TDateInterval tlapse)
{
	TDate	toTDate;
	double	yearFrac;

	switch (toupper(tlapse.prd_typ)) {
       		case 'A':
			return (double) tlapse.prd;
       		case 'S':
			return (double) tlapse.prd / 2.0;
       		case 'Q':
			return (double) tlapse.prd / 4.0;
       		case 'M':
			return (double) tlapse.prd / 12.0;
		default:
			GtoDtFwdAny(aDate, &tlapse, &toTDate);
			GtoDayCountFraction(aDate, toTDate,
				GTO_ACT_365F, &yearFrac);
			return yearFrac;
		}
}


/*f-------------------------------------------------------------
 * Date routines : add years to date interval.
 *                                                          
 * <br><br>
 * Adds a (possibly fractional) number of years "tt"
 * to a date by rounding it to the closest number of days.
 * 
 */

DLL_EXPORT(int)
DrlTDateAddYears(TDate startDate, double tt, TDate *endDate)
{
#define	HALFDAY	(1.3698e-3)
static	TDateInterval	interval;
	int		ny, nm, nd;
	/*
	 *
	 */
	ny = (int) floor(tt + HALFDAY);
	nm = (int) floor(12e0 * (tt - (double) ny + HALFDAY));
	nd = (int) (365e0 * (tt - (double) ny - (double) nm / 12e0
		+ HALFDAY));

	*endDate = startDate;

	if (ny != 0) {
		GtoMakeDateInterval(ny, 'A', &interval);
		GtoDtFwdAny(*endDate, &interval, endDate);
	}

	if (nm != 0) {
		GtoMakeDateInterval(nm, 'M', &interval);
		GtoDtFwdAny(*endDate, &interval, endDate);
	}

	if (nd != 0) {
		endDate += (long) nd;
		/*GtoMakeDateInterval(ny, 'A', &interval);
		GtoDtFwdAny(*endDate, &interval, endDate);*/
	}

	return(SUCCESS);
#undef	HALFDAY
}

#ifdef	_SKIP
/*--------------------------------------------------------------
 * Adds a (possibly fractional) number of years "tt"
 * to a date interval by rounding it to the closest number of days.
 */

DLL_EXPORT(TDate)
DrlTDateAddTime(TDate startDate, double tt)
{
	TDateInterval	interval;
	TDate	endDate = startDate;
	int	ny, nm, nd;
#define	HALFDAY	(1.3698e-3)
	/*
	 *
	 */
	ny = (int) floor(tt + HALFDAY);
	nm = (int) floor(12e0 * (tt - (double) ny + HALFDAY));
	nd = (int) (365e0 * (tt - (double) ny - (double) nm / 12e0
		+ HALFDAY));

	if (ny != 0) {
		GtoMakeDateInterval(ny, 'A', &interval);
		GtoDtFwdAny(endDate, &interval, &endDate);
	}

	if (nm != 0) {
		GtoMakeDateInterval(nm, 'M', &interval);
		GtoDtFwdAny(endDate, &interval, &endDate);
	}

	if (nd != 0) {
		endDate += (long) nd;
		/*GtoMakeDateInterval(ny, 'A', &interval);
		GtoDtFwdAny(endDate, &interval, &endDate);*/
	}

	return(endDate);
#undef	HALFDAY
}
#endif


#ifdef	_SKIP

/*--------------------------------------------------------------
 * Date routines : convert date interval to year fraction.
 *                                                          
 * <br><br>
 * Converts a date interval "tlapse" into a number of years.
 * This function works for ``date-independent'' date intervals
 * (such as days, months, years, etc\dots). It will not
 * work for periods such as IMM, end of months,\dots
 */

DLL_EXPORT(double)
DrlTDateIntervalToDouble(TDateInterval tlapse)
{
	TDate	fromTDate, toTDate;
	double	yearFrac;

	switch (toupper(tlapse.prd_typ)) {
       		case 'A':
			return (double) tlapse.prd;
       		case 'S':
			return (double) tlapse.prd / 2.0;
       		case 'Q':
			return (double) tlapse.prd / 4.0;
       		case 'M':
			return (double) tlapse.prd / 12.0;
		default:
			fromTDate = DrlTDateMake(1,1,1995);
			GtoDtFwdAny(fromTDate, &tlapse, &toTDate);
			GtoDayCountFraction(fromTDate, toTDate,
				GTO_B30_360, &yearFrac);
			return yearFrac;
		}
}
#endif

/*f-------------------------------------------------------------
 * Date routines : advance date interval by years.
 *                                                          
 * <br><br>
 * Converts a date interval <i> interval</i> into a number of years
 * using <i> baseDate</i> as reference date.
 */

DLL_EXPORT(int)
DrlTDateAdvanceYears(
	TDate baseDate,		/* (I) reference date */
	double years,		/* (I) year fraction */
	TDayCount dayCount,	/* (I) day count for year fraction */
	TDate *newDate)		/* (O) */
{
static	char	routine[] = "DrlTDateAdvanceYears";
	TDateInterval	interval;
	int		months, days;

#define	HALF_DAY	0.001370

	switch (dayCount) {
	case GTO_ACT_365F:
		*newDate = baseDate + (long) floor(years*365e0 + 0.5e0);
		return (SUCCESS);
	case GTO_B30_360:

		/* get number of months */
		months = (int) floor(years * 12e0 + HALF_DAY);
		GtoMakeDateInterval(months, 'M', &interval);
		if (GtoDtFwdAny(baseDate, &interval, newDate) ISNT SUCCESS)
			goto done;

		/* advance days */
		days = (int)floor((years - (double)months/12e0)
			* 365e0 + 0.5e0);
		GtoMakeDateInterval(days, 'D', &interval);
		if (GtoDtFwdAny(*newDate, &interval, newDate) ISNT SUCCESS)
			goto done;
		return(SUCCESS);
	default:
		GtoErrMsg("%s: day count not available.\n", routine);
		return(FAILURE);
	}
#undef	HALF_DAY
done:
	*newDate = -1L;
	return(FAILURE);
}


/*f-------------------------------------------------------------
 * Date routines : convert years to date interval.
 *                                                          
 * <br><br>
 * Converts a number of years (expressed as ACT/365) to date interval
 * using <i> baseDate</i> as reference date.
 * This function also works for ``date-dependent'' intervals
 * (such as IMM dates, etc\dots).
 */

DLL_EXPORT(int)
DrlTDateIntervalFromYears(
	TDate baseDate,		/* (I) reference date  */
	double years,		/* (I) year fraction */
	TDayCount dayCount,	/* (I) day count for year fraction */
	TDateInterval *interval)/* (O) */
{
	TDate		expDate, newDate;
	int		numInt;


	if (DrlTDateAdvanceYears(baseDate, years, dayCount, &expDate)
		!= SUCCESS)
			goto error;

	/* check for A */
	numInt = (int) floor(years + 0.5e0);
	GtoMakeDateInterval(numInt, 'A', interval);
	GtoDtFwdAny(baseDate, interval, &newDate);
	if (expDate == newDate) return(SUCCESS);

	/* check for M */
	numInt = (int) floor(years*12e0 + 0.5e0);
	GtoMakeDateInterval(numInt, 'M', interval);
	GtoDtFwdAny(baseDate, interval, &newDate);
	if (expDate == newDate) return(SUCCESS);

	/* check for I */
	numInt = (int) ceil(years*4e0+ 0.5e0);
	GtoMakeDateInterval(numInt, 'I', interval);
	GtoDtFwdAny(baseDate, interval, &newDate);
	if (expDate == newDate) return(SUCCESS);

	/* default is D */
	numInt = (int) floor(years*365e0 + 0.5e0);
	GtoMakeDateInterval(numInt, 'D', interval);
	GtoDtFwdAny(baseDate, interval, &newDate);
	if (expDate == newDate) return(SUCCESS);


	/* Error value */
error:
	GtoErrMsg("DrlTDateIntervalFromYears: failed.\n");
	GtoMakeDateInterval(-100, 'A', interval);
	return(FAILURE);
}


/*f-------------------------------------------------------------
 * Date routines : convert date interval to years.
 *                                                          
 * <br><br>
 * Converts a date interval to a number of years (expressed as ACT/365).
 * This function works for ``date-dependent'' intervals
 * (such as IMM dates, etc\dots).
 */

DLL_EXPORT(int)
DrlTDateIntervalToYears(
	TDate baseDate,		/* (I) reference date  */
	TDateInterval interval,	/* (I) input interval */
	TDayCount dayCount,	/* (I) day count for year fraction */
	double *years)		/* (O) */
{
	TDate		expDate;
	if (GtoDtFwdAny(baseDate, &interval, &expDate) != SUCCESS)
		goto error;
	if (GtoDayCountFraction(baseDate, expDate, dayCount, years) != SUCCESS)
		goto error;
	return(SUCCESS);
error:
	GtoErrMsg("DrlTDateIntervalToYears: failed.\n");
	return(FAILURE);
}



/*f-------------------------------------------------------------
 * Date routines : scan date interval in char string.
 *                                                          
 * <br><br>
 * Scans for a date interval in the string "s". If scan is
 * successful, returns 0 and puts the result in "tlapse".
 * If failure, returns a non-zero value.
 * Recognized formats are of the form ``MXNY..'' where
 * M, N are number of periods (non-negative) and X, Y
 * a character defining a period (D,W,M,Q,etc\dots).
 */

DLL_EXPORT(int)
DrlTDateIntervalScan(char *s, TDateInterval *tlaps)
{
static	char	routine[]="DrlTDateIntervalScan";
	char		*p, *q,  buf[256], buf2[256];
	int		nPer;
	TDateInterval	intval;

	/* Special case */
	if ((!strcmp(s, "O/N")) ||
	    (!strcmp(s, "ON"))  ||
	    (!strcmp(s, "SN"))  ||
	    (!strcmp(s, "S/N"))) {
		if (GtoMakeDateInterval(1, 'D', tlaps) != SUCCESS)
			goto done;
		return(SUCCESS);
	}

	/* Scan for IMMn intervals */
	if ((p = strstr(s, "IMM")) != NULL) {
		p +=3;
		if (sscanf(p, "%d", &nPer) != 1)
			goto done;
		if (GtoMakeDateInterval(nPer, 'I', tlaps) != SUCCESS)
			goto done;
		return(SUCCESS);
	}

	/* replace '3y' by '3A' for GtoString... */
	strcpy(buf, s);
	for (p = buf; *p != '\0'; p++) {
		if (toupper(*p) == 'Y') *p = 'A';
	}

#define	_OLDVERSION
#ifdef	_OLDVERSION
	return (GtoStringToDateInterval(buf, routine, tlaps));
#endif

	if (GtoMakeDateInterval(0, 'M', tlaps) != SUCCESS)
		goto done;

	p = buf;
	while (*p != '\0') {
		q = buf2;
		while (isdigit(*p)) *q++ = *p++;
		if (*p == '\0') goto done;
		while (isalpha(*p)) *q++ = *p++;
		*q = '\0';

		if (GtoStringToDateInterval(buf2, routine, &intval) != SUCCESS)
			goto done;

		if (GtoIntervalAddApprox(tlaps, &intval, tlaps) != SUCCESS)
			goto done;
	}
	return(SUCCESS);


	return(SUCCESS);
done:
	GtoErrMsg("DrlTDateIntervalScan: can't scan `%s'.\n", s);
	return(FAILURE);
}


/*f-------------------------------------------------------------
 * Date routines : print date interval to char string.
 *                                                          
 * <br><br>
 * Prints the date interval "tlapse" in the buffer "s".
 * If "buf" is NULL, returns a pointer to a static char string.
 */

DLL_EXPORT(char*)
DrlTDateIntervalPrint(char *string, TDateInterval tlapse)
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

	strcpy(s, GtoFormatDateInterval(&tlapse));
	return(s) ;
#undef	MAX_IDX
}


/*f-------------------------------------------------------------
 * Date routines : convert frequency to date interval.
 *                                                          
 * <br><br>
 * Converts a frequency "freq" (1,2,4,12) to a date
 * interval (1A,1S,1Q,1M).
 */

DLL_EXPORT(int)
DrlFreq2TDateInterval(int freq, TDateInterval *interval)
{
static	char	routine[] = "DrlFreq2TDateInterval";
	switch (freq) {
	case 1:
		*interval = DrlTDateIntervalFromDMY(0, 0, 1);
		return(SUCCESS);
	case 2:
		*interval = DrlTDateIntervalFromDMY(0, 6, 0);
		return(SUCCESS);
	case 4:
		*interval = DrlTDateIntervalFromDMY(0, 3, 0);
		return(SUCCESS);
	case 12:
		*interval = DrlTDateIntervalFromDMY(0, 1, 0);
		return(SUCCESS);
	default:
		GtoErrMsg("%s: bad freq %d\n", routine, freq);
		return(FAILURE);
	}
}

/*f-------------------------------------------------------------
 * Date routines : convert date interval to frequency.
 *                                                          
 * <br><br>
 * Converts a date interval "laps" (1A,1S,1Q,1M)
 * to frequency (1,2,4,12).
 */

DLL_EXPORT(int)
DrlTDateInterval2Freq(TDateInterval laps, int *freq)
{
static	char	routine[] = "DrlTDateInterval2Freq";
	TDate	refDate = DrlTDateMake(1,1,1990);

	*freq = (int) floor(1/drlTDateIntervalDiff(refDate, laps) + 0.5);
	switch (*freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		return(SUCCESS);
	default:
		GtoErrMsg("%s: bad freq (%s->%d)\n",
			routine, DrlTDateIntervalPrint(NULL, laps), freq);
		return(FAILURE);
	}
}

/*f-------------------------------------------------------------
 * Date routines : check valid frequency.
 *                                                          
 * <br><br>
 * Checks if frequency  "freq" is valid
 * Valid frequencies are 1, 2, 4 or 12.
 * Returns 0 (SUCCESS) if frequency is valid.
 */

DLL_EXPORT(int)
DrlCheckFreqValid(int freq)
{
	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
	    return 0;
	default:
	    GtoErrMsg("CheckFreqValid: invalid frequency %d\n",
		freq);
	    return 1;
	}
}


/*f---------------------------------------------------------------------
 * Date routines : advance date by date interval.
 *                                                          
 * <br><br>
 * Finds the first date after <i> refDate</i> plus <i> numNotifDays</i>
 * that is obtained by adding an integer number of <i> interval</i>
 * to <i> oldDate</i>. Puts the result in <i> newDate</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlTDateAdvanceToRefDate(
	TDate oldDate,		/* (I) original date */
	int numNotifDays,	/* (I) # of notification days */
	TDateInterval interval,	/* (I) offset interval */
	TDate refDate,		/* (I) referance date */
	TDate *newDate)		/* (O) new date */
{
	refDate += (long) numNotifDays;
	*newDate = oldDate;
	while (*newDate < refDate) {
	    if (GtoDtFwdAny(*newDate, &interval, newDate) != SUCCESS)
		return(FAILURE);
	}
	return(SUCCESS);
}



