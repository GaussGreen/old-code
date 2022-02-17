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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#ifdef	DRL_CLIB
#include "date_sup.h"		/* C analytics */
#endif

#include "drltime.h"		/* prototype consistency */




#if !defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)
/*f-------------------------------------------------------------
 * Date routines : create date interval.
 *                                                          
 * <br><br>
 * Creates a number "nDays" of days, a number "nMonths" of months
 * and a number "nYears" of years in a date interval.
 * Returns the date interval created.
 */


int
DrlDIntervalSet(
	int numPeriods,                    /* (I) Number of periods */
	char periodType,                   /* (I) Period type */
	DInterval *interval)           /* (O) Value read from file */
{
    static char routine[]="DrlDIntervalSet";    

    interval->flag = 0;
    
    switch(toupper(periodType))
    {
        case 'A':                       /* Year */
        case 'Y':                       /* Year */
            interval->prd_typ = 'M';
            interval->prd = numPeriods * 12;
            break;

        case 'S':                       /* Semi annual */
            interval->prd_typ = 'M';
            interval->prd = numPeriods * 6;
            break;
            
        case 'Q':                       /* Quarter */
            interval->prd_typ = 'M';
            interval->prd = numPeriods * 3;
            break;

        case 'W':                       /* Week */          
            interval->prd_typ = 'D';
            interval->prd = numPeriods * 7;
            break;

        case 'D':                       /* Day */
        case 'M':                       /* Normal Month */
        case 'F':                       /* Flexible End of month */ 
        case 'E':                       /* End of month unconditional */
        case 'G':                       /* 29th of month */
        case 'H':                       /* 30th of month */      
        case 'I':                       /* IMM date */
        case 'J':                       /* monthly IMM date */
        case 'K':                       /* Aussie quarterly IMM date */
            interval->prd_typ = (char) toupper(periodType);         
            interval->prd = numPeriods;
            break;
            
        default:
            if (isprint((int)periodType))    /* If periodType is printable */
            {
                DrlErrMsg("%s: Interval type %c not valid.\n",
                          routine, periodType);
            }
            else
            {
                DrlErrMsg("%s: Interval type (unprintable) not valid.\n",
                          routine);
            }
            return(FAILURE);
    }
    return(SUCCESS);
}

#endif /*!defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)*/

DInterval
DrlDIntervalNil(void)
{
	DInterval	tlaps;
	tlaps.prd_typ = '\0';
	return(tlaps);
}


/*--------------------------------------------------------------
 * Allocates an array of date intervals of size "size" and returns
 * a pointer to it. Returns NULL if failure.
 */

DInterval*
DrlDIntervalArrayAlloc(int size)
{
	return (DInterval*) MALLOC(size*sizeof(DInterval));
}

/*--------------------------------------------------------------
 * Frees an array previousely allocated by <i> DrlDIntervalArrayAlloc</i>.
 */

int
DrlDIntervalArrayFree(DInterval* ptr)
{
	FREE((void*)ptr);
	return(0);
}


#if !defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)
/*f-------------------------------------------------------------
 * Date routines : advance date interval by years.
 *                                                          
 * <br><br>
 * Converts a date interval <i> interval</i> into a number of years
 * using <i> baseDate</i> as reference date.
 */

int
DrlDDateAdvanceYears(
	DDate baseDate,		/* (I) reference date */
	double years,		/* (I) year fraction */
	DDate *newDate)		/* (O) */
{
static	char	routine[] = "DrlDDateAdvanceYears";
	DInterval	interval;
	int		months, days;

#define	HALF_DAY	0.001370

	/* get number of months */
	months = (int) floor(years * 12e0 + HALF_DAY);
	IF_FAILED_DONE( DrlDIntervalSet(months, 'M', &interval));
	IF_FAILED_DONE( DrlDDateFwdAny(baseDate, &interval, newDate));

	/* advance days */
	days = (int)floor((years - (double)months/12e0)
		* 365e0 + 0.5e0);
	IF_FAILED_DONE( DrlDIntervalSet(days, 'D', &interval));
	IF_FAILED_DONE( DrlDDateFwdAny(*newDate, &interval, newDate));
	return(SUCCESS);
#undef	HALF_DAY
done:
	DrlErrMsg("%s: failed.\n", routine);
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

int
DrlDIntervalFromYears(
	DDate baseDate,		/* (I) reference date  */
	double years,		/* (I) year fraction */
	DInterval *interval)/* (O) */
{
	DDate		expDate, newDate;
	int		numInt;


	if (DrlDDateAdvanceYears(baseDate, years, &expDate) != SUCCESS)
			goto done;

	/* check for A */
	numInt = (int) floor(years + 0.5e0);
	IF_FAILED_DONE( DrlDIntervalSet(numInt, 'A', interval));
	IF_FAILED_DONE( DrlDDateFwdAny(baseDate, interval, &newDate));
#if !defined(DRL_TDATE_DMY)
	if (expDate == newDate) return(SUCCESS);
#else
	NOT_IMPLEMENTED;
#endif

	/* check for M */
	numInt = (int) floor(years*12e0 + 0.5e0);
	IF_FAILED_DONE( DrlDIntervalSet(numInt, 'M', interval));
	IF_FAILED_DONE( DrlDDateFwdAny(baseDate, interval, &newDate));
#if !defined(DRL_TDATE_DMY)
	if (expDate == newDate) return(SUCCESS);
#else
	NOT_IMPLEMENTED;
#endif

	/* check for I */
	numInt = (int) ceil(years*4e0+ 0.5e0);
	IF_FAILED_DONE( DrlDIntervalSet(numInt, 'I', interval));
	IF_FAILED_DONE( DrlDDateFwdAny(baseDate, interval, &newDate));
#if !defined(DRL_TDATE_DMY)
	if (expDate == newDate) return(SUCCESS);
#else
	NOT_IMPLEMENTED;
#endif

	/* default is D */
	numInt = (int) floor(years*365e0 + 0.5e0);
	IF_FAILED_DONE( DrlDIntervalSet(numInt, 'D', interval));
	IF_FAILED_DONE( DrlDDateFwdAny(baseDate, interval, &newDate));
#if !defined(DRL_TDATE_DMY)
	if (expDate == newDate) return(SUCCESS);
#else
	NOT_IMPLEMENTED;
#endif


	/* Error value */
done:
	DrlErrMsg("DrlDIntervalFromYears: failed.\n");
	return(FAILURE);
}
#endif /*!defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)*/


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

int
DrlDIntervalScan(char *s, DInterval *tlaps)
{
static	char	routine[]="DrlDIntervalScan";
	char		*p, *q,  buf[256], buf2[256];

        char		*inp=s;                    /* Pointer to input */
#define MAX_DIGITS 32
        char		numberBuff[MAX_DIGITS];
        char		*nump = numberBuff;            /* Pointer to number */
        char		periodType;
        int		numPeriods;

	/* Special case */
	if ((!strcmp(s, "O/N")) ||
	    (!strcmp(s, "ON"))  ||
	    (!strcmp(s, "SN"))  ||
	    (!strcmp(s, "S/N"))) {
		if (DrlDIntervalSet(1, 'D', tlaps) != SUCCESS)
			goto done;
		return(SUCCESS);
	}

	/* Scan for IMMn intervals */
	if ((p = strstr(s, "IMM")) != NULL) {
		p +=3;
		if (sscanf(p, "%d", &numPeriods) != 1)
			goto done;
		if (DrlDIntervalSet(numPeriods, 'I', tlaps) != SUCCESS)
			goto done;
		return(SUCCESS);
	}

	/* replace '3y' by '3A' for GtoString... */
	strcpy(buf, s);
	for (p = buf; *p != '\0'; p++) {
		if (toupper(*p) == 'Y') *p = 'A';
	}

#ifdef	DRL_CLIB
	return (GtoStringToDateInterval(buf, routine, tlaps));
#else

    
        /* Copy sign,if any, to numberBuff */
        if ((*inp == '-') || (*inp == '+'))
        {
            *nump++ = *inp++;
        }
    
        /* Copy digits, if any, to numberBuff */
        while (isdigit((int)*inp))   
            *nump++ = *inp++;
        *nump = '\0';                       /* Null terminate */
    
        if (inp != s)                 /* Found some digits */
            numPeriods = atoi(numberBuff);
        else                                /* Found on digits */
        {
            numPeriods = 1;
        }
    
        periodType = (char)toupper(*inp);   /* To upper case */
        
        IF_FAILED_DONE( DrlDIntervalSet(numPeriods, periodType, tlaps));

#endif


	return(SUCCESS);
done:
	DrlErrMsg("DrlDIntervalScan: can't scan `%s'.\n", s);
	return(FAILURE);
#undef	MAX_DIGITS
}


/*f-------------------------------------------------------------
 * Date routines : print date interval to char string.
 *                                                          
 * <br><br>
 * Prints the date interval "tlapse" in the buffer "s".
 * If "buf" is NULL, returns a pointer to a static char string.
 */

char*
DrlDIntervalPrint(char *string, DInterval *interval)
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][64] ;
static	int	tmpIdx=0;
	char	*s ;
        char periodType;
        int numPeriods;

	if (string == NULL) {
		s = tmp[tmpIdx];
		tmpIdx++;
		if (tmpIdx > MAX_IDX-1) tmpIdx=0;
	} else {
		s = string;
	}
#ifdef	DRL_CLIB
	strcpy(s, GtoFormatDateInterval(interval));
#else
        switch(interval->prd_typ)
        {
            case 'M':                       /* Months */
                if (interval->prd % 12 == 0)
                {
                    periodType = 'A';
                    numPeriods = interval->prd/12;
                }
                else if (interval->prd % 6 == 0)
                {
                    periodType = 'S';
                    numPeriods = interval->prd/6;
                }
                else if (interval->prd % 3 == 0)
                {
                    periodType = 'Q';
                    numPeriods = interval->prd/3;
                }
                else
                {
                    periodType = 'M';
                    numPeriods = interval->prd;
                }
                break;
                
            case 'D':                       /* Days */
                if (interval->prd % 7 == 0)
                {
                    periodType = 'W';
                    numPeriods = interval->prd/7;
                }
                else
                {
                    periodType = 'D';
                    numPeriods = interval->prd;
                }
                
            default:
                periodType = interval->prd_typ;
                numPeriods = interval->prd;
        }
    
	sprintf(s, "%d%c", numPeriods, periodType);
	/*sprintf(s, "%d%c", tlaps->prd, tlaps->prd_typ);*/
#endif
	return(s) ;
#undef	MAX_IDX
}


#if !defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)
/*f-------------------------------------------------------------
 * Date routines : convert date interval to years.
 *                                                          
 * <br><br>
 * Converts a date interval to years (remark: IMM periods
 * are treated as quarters, etc.).
 */

int
DrlDIntervalToYears(
	DInterval *interval,	/* (I) */
	double *years)			/* (O) # times per year */
{
static	char	routine[]="DrlDIntervalToYears";

    switch(toupper(interval->prd_typ))
    {
        case 'A':
        case 'Y':
            *years = (double) interval->prd;
            break;
        case 'S':
            *years = (double) interval->prd / 2e0;
            break;
        case 'Q':
        case 'I':
        case 'K':
            *years = (double) interval->prd / 4e0;
            break;
        case 'M':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'J':
            *years = (double) interval->prd / 12e0;
            break;
        case 'W':
            *years = (double) interval->prd * 7e0/365e0;
            break;
        case 'D':
            *years = (double) interval->prd / 365e0;
            break;
        default:
            DrlErrMsg("%s: unknown interval type %c.\n",
                      routine, interval->prd_typ);
            return(FAILURE);
    }
    return(SUCCESS);
}


/*f-------------------------------------------------------------
 * Date routines : convert year o date interval.
 *                                                          
 * <br><br>
 * Converts a date interval to years.
 */

int
DrlYearsToDInterval(
	double years,			/* (I) # times per year */
	DInterval *interval)	/* (O) */
{
	double	monthsDouble;
	int	months, days;
const	double	HALF_DAY = .001370e0;	/* In terms of years */

	/* Get right number of months */
	monthsDouble = years * 12e0;

	/* Check if this can be expressed in terms of months.  */
	if (fabs(monthsDouble - floor(monthsDouble)) < HALF_DAY)
	{
		months = (int)floor(monthsDouble + HALF_DAY);
		IF_FAILED_DONE( DrlDIntervalSet(months, 'M', interval));
	}
	else                                /* Instead use days */
	{
		days = (int)floor(years * 365e0 + .5);
		IF_FAILED_DONE( DrlDIntervalSet(days, 'D', interval));
	}
	return(SUCCESS);
done:
	return(FAILURE);
}






/*f-------------------------------------------------------------
 * Date routines : convert frequency to date interval.
 *                                                          
 * <br><br>
 * Converts a frequency "freq" (1,2,4,12) to a date
 * interval (1A,1S,1Q,1M).
 */

int
DrlFreqToDInterval(int freq, DInterval *interval)
{
static	char	routine[] = "DrlFreqToDInterval";
	switch (freq) {
	case 1:
		return DrlDIntervalSet(1, 'A', interval);
	case 2:
		return DrlDIntervalSet(1, 'S', interval);
	case 4:
		return DrlDIntervalSet(1, 'Q', interval);
	case 12:
		return DrlDIntervalSet(1, 'M', interval);
	default:
		DrlErrMsg("%s: bad freq %d.\n", routine, freq);
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

int
DrlDIntervalToFreq(DInterval *interval, int *freq)
{
static	char	routine[] = "DrlDIntervalToFreq";
	double	years;
  
	if (DrlDIntervalToYears(interval, &years) != SUCCESS)
	{
		DrlErrMsg("%s: failed.\n", routine);
		return(FAILURE);
	}

	if (IS_ALMOST_ZERO(years))
	{
		DrlErrMsg("%s: interval is zero.\n", routine);
		return(FAILURE);
	}

	*freq = (int) floor(1e0 / years);
	switch (*freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		return(SUCCESS);
	default:
		DrlErrMsg("%s: bad freq (%s->%d)\n",
			routine, DrlDIntervalPrint(NULL, interval), freq);
		return(FAILURE);
	}
}

#endif /*!defined(DRL_CLIB) && !defined(DRL_TDATE_MDY)*/

/*f-------------------------------------------------------------
 * Date routines : check valid frequency.
 *                                                          
 * <br><br>
 * Checks if frequency  "freq" is valid
 * Valid frequencies are 1, 2, 4 or 12.
 * Returns 0 (SUCCESS) if frequency is valid.
 */

int
DrlCheckFreqValid(int freq)
{
	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
	    return 0;
	default:
	    DrlErrMsg("CheckFreqValid: invalid frequency %d\n",
		freq);
	    return 1;
	}
}


