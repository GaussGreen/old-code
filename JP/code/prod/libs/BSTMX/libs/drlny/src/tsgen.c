/****************************************************************
 * Module:	DRL
 * Submodule:	TS
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <float.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

# include "date_sup.h"
# include "yearfrac.h"
# include "convert.h"
# include "check.h"
# include "zr2coup.h"		/* Analytics C Library */
# include "zr2simp.h"		/* Analytics C Library */
# include "duration.h"
# include "convex.h"
# include "swapadj.h"

#include "drlio.h"
#include "drlvtype.h"	/* LIL interface */
#include "drltime.h"	/* DrlTDateScanYMD */
#include "drlstr.h"	/* for sscanf */
#include "drlproc.h"	/* DrlStrerror */

#include "drlts.h"


/*f--------------------------------------------------------------
 * Reads a TCurve in file "fnam". The result is
 * placed in the structure "that".
 * The format of the file is
 * \begin{verbatim}
 * <value date>
 * <number of points>
 * <compounding>
 * <day count>
 * <date1>  <rate1>
 * ...
 * <dateN>  <rateN>
 * \end{verbatim}
 * An example of valid format is
 * \begin{verbatim}
 * 09/03/1994
 * 71
 * 1.000000
 * Act/365
 *         09/04/1994      4.667885
 *         10/03/1994      4.659252
 *         11/03/1994      4.847664
 *         12/03/1994      4.969433
 *         01/03/1995      5.038225
 *         02/03/1995      5.156677
 *         ...[list of 71 points]
 *         09/03/2026      8.263582
 * \end{verbatim}
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveFileRead(TCurve **that, char *fnam, long fmt)
{
	FILE		*fp;
	int		errCode ;
static	char		routine[] = "DrlTCurveFileRead";


	if ((fp = fopen(fnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open file `%s' (%s)\n",
			routine, fnam, DrlStrerror());
		errCode = -1;
		return errCode;
	} else {
		errCode = DrlTCurveFpRead(that, fp, fmt);
		if (errCode != 0)
			GtoErrMsg("%s: can't read file `%s' (code %d)\n",
				routine, fnam, errCode);
		fclose(fp);
		return(errCode);
	}
}

/*f--------------------------------------------------------------
 * Writes a TCurve "that" in file "fnam".
 * The format of the file is the same as for <i> DrlTCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveFileWrite(TCurve *that, char *fnam, long fmt)
{
	FILE		*fp;

	if ((fp = fopen(fnam, "w")) == NULL) {
		GtoErrMsg("DrlTCurveFileWrite: can't open file `%s' (%s)\n",
			fnam, DrlStrerror());
		return(-1);
	} else {
		DrlTCurveFpWrite(that, fp, fmt);
		fclose(fp);
		return(0);
	}
}



/*f--------------------------------------------------------------
 * Reads a TCurve "that" from a file pointer "fp".
 * The format of the file is the same as for <i> DrlTCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveFpRead(TCurve **that, FILE *fp, long fmt)
{
static	char		routine[] = "DrlTCurveFpRead";
	char		buf[255], buf2[255];
	int		i, line = 0, c,
			errCode = FAILURE;
	TDate		baseDate;
	int		nPts, nDays;
	double		basis;
	TDayCount	dayCountConv;

	*that = NULL;


	/*
	 * Identify format
	 */
	if ((c = fgetc(fp)) == EOF) goto done;
	ungetc(c, fp);
	if (c != '#') goto cdfmt;

	/*
	 * (1) London format:
	 * # Start date
	 * 19970411
	 * # Money Market basis (360 or 365)
	 * 360
	 * # Annual or semi-annual curve ("A" or "S")
	 * S
	 * # Year basis for benchmark swaps ("ACT", "365" or "360")
	 * ACT
	 * # No of entries
	 * 84
	 * #zero maturity yyyymmdd rates (ACT/365F annual)
	 * 19970412 5.734349908886
	 * 19970511 5.913551526763
	 * etc
	 */
	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL)
		goto done;
	if (DrlTDateScanYMD(buf, &baseDate) != 0) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (sscanf(buf, "%d", &nPts) != 1) 
		goto done;

	basis = 1L;
	dayCountConv = GTO_ACT_365F;

	/* create new TCurve */
	*that = GtoNewTCurve(baseDate, nPts, basis, dayCountConv);
	if (*that == NULL) goto done;

	if (nPts > 0) {
	    for (i=0; i<=nPts-1; i++) {
		if ((DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) ||
		    (sscanf(buf, "%s %lf", buf2, &DrlTCurveRate(*that,i)) != 2) ||
		    (DrlTDateScanYMD(buf2, &DrlTCurveDate(*that, i)) != 0)) {
			goto done;
		}
		DrlTCurveRate(*that, i) *= 1e-2;
	    }
	}

	/* made it through */
	errCode = SUCCESS;
	goto done;


	/*
	 * (2) CD format
	 */
cdfmt:
	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL)
		goto done;
	if (GtoStringToDate(buf, &baseDate) != 0) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (sscanf(buf, "%d", &nPts) != 1) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (sscanf(buf, "%lf", &basis) != 1) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (GtoStringToDayCountConv(buf, &dayCountConv) != SUCCESS) 
		goto done;


	/* create new TCurve */
	*that = GtoNewTCurve(baseDate, nPts, basis, dayCountConv);
	if (*that == NULL) goto done;

	if (nPts > 0) {
	    for (i=0; i<=nPts-1; i++) {
		line++;

		/* scan '#days offset   value' */
		if (DrlFScanInt(fp, &nDays) != SUCCESS)
			goto done;
		if (nDays < 0) {
			GtoErrMsg("%s: item %d has negative number of "
				"days from value date.\n", 
				routine, i+1);
			goto done;
		}
		DrlTCurveDate(*that, i) = baseDate + (TDate)nDays;


		if (DrlFScanDouble(fp, &DrlTCurveRate(*that,i)) != SUCCESS)
			goto done;

		switch (fmt) {
		case DRL_TCURVE_FMT_PERCENT:
			DrlTCurveRate(*that, i) *= 1e-2;
			break;
		case DRL_TCURVE_FMT_STD:
		case DRL_TCURVE_FMT_CUR:
			break;
		default:
			GtoErrMsg("%s: invalid format.\n", routine);
			goto done;
		}
	    }
	}

	/* made it through */
	errCode = SUCCESS;
done:
	if (errCode != SUCCESS) {
	    if (*that != NULL) GtoFreeTCurve(*that);
	    *that = NULL;
	    GtoErrMsg("%s: failed (line %d).\n", routine, line);
	}
	return(errCode);
}


/*f--------------------------------------------------------------
 * Writes a TCurve "that" to a file pointer "fp".
 * The format of the file is the same as for <i> DrlTCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveFpWrite(TCurve *that, FILE *fp, long fmt)
{
static	char	routine[] = "DrlTCurveFpWrite";
	int	i, cnt;


	DrlFPrintf(fp, "%s\n",  GtoFormatDate(DrlTCurveBaseDate(that)));
	DrlFPrintf(fp, "%d\n",  DrlTCurveNumItems(that));
	DrlFPrintf(fp, "%lf\n", DrlTCurveBasis(that));
	DrlFPrintf(fp, "%s\n",  GtoFormatDayCountConv(
				DrlTCurveDayCountConv(that)));

	switch (fmt) {
	case DRL_TCURVE_FMT_PERCENT:
		for (i=0; i<=DrlTCurveNumItems(that)-1; i++) {
		    DrlFPrintf(fp, "\t%4d\t%8.6f    # %10s\n",
			(int) (DrlTCurveDate(that, i) - that->fBaseDate),
			DrlTCurveRate(that, i) * 1e2,
			DrlTDatePrint(NULL, DrlTCurveDate(that, i)));
		}
	break;
	case DRL_TCURVE_FMT_STD:
		for (i=0; i<=DrlTCurveNumItems(that)-1; i++) {
		    DrlFPrintf(fp, "\t%4d\t%8.6f\n",
			(int) (DrlTCurveDate(that, i) - that->fBaseDate),
			DrlTCurveRate(that, i));
		}
	break;
	case DRL_TCURVE_FMT_CUR:
		for (i=0; i<=DrlTCurveNumItems(that)-1; i++) {
		    DrlFPrintf(fp, "\t%4d\t%15s\n",
			(int) (DrlTCurveDate(that, i) - that->fBaseDate),
			DrlCurPrint(NULL, DrlTCurveRate(that, i), 2));
		}
	break;

	case DRL_TCURVE_FMT_REGTEST:
		DrlFPrintf(fp, "ZC_BASE_DATE: %s ",
			GtoFormatDate(DrlTCurveBaseDate(that)));
		DrlFPrintf(fp, "\t%%ARRAY_END%% %%ARG_END%%\n");

		DrlFPrintf(fp, "ZC_DATES:\n");
		cnt = 0;
		for (i=0; i<=DrlTCurveNumItems(that)-1; i++) {
		    DrlFPrintf(fp, "\t%s",
			GtoFormatDate(DrlTCurveDate(that, i)));
		    if (++cnt == 5) { DrlFPrintf(fp, "\n"); cnt=0;}
		}
		DrlFPrintf(fp, "\t%%ARRAY_END%% %%ARG_END%%\n");

		DrlFPrintf(fp, "ZC_RATES:\n");
		cnt = 0;
		for (i=0; i<=DrlTCurveNumItems(that)-1; i++) {
		    DrlFPrintf(fp, "\t%8.6f",
			DrlTCurveRate(that, i));
		    if (++cnt == 5) { DrlFPrintf(fp, "\n"); cnt=0;}
		}
		DrlFPrintf(fp, "\t%%ARRAY_END%% %%ARG_END%%\n");

	break;
	case 9999:
		for (i=0; i<=DrlTCurveNumItems(that)-1; i++) {
		    DrlFPrintf(fp, "\t%10s\t%15.8f\n",
			DrlTDatePrint(NULL, DrlTCurveDate(that, i)),
			DrlTCurveRate(that, i));
		}
	break;

	default:
		GtoErrMsg("%s: invalid format.\n", routine);
		return(FAILURE);
	}

	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * Reads a TCurve "that" from a wrapper interface.
 * On entry, "refDate" is the array of length 1 containg the
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates (annual compounding, Act/365).
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveWrapRead(
	TCurve	**that,		/* (O) new curve */
	TDateL	*refDate,	/* (I) */
	TDateL	*zDate,		/* (I) */
	FloatL	*zRate)		/* (I) */
{
static	char	routine[] = "DrlTCurveWrapRead";
	int	status = FAILURE;
	int	i, n1, n2,
		nZero;


	if (((int) refDate[0] != 1) ||
	    ((n1 = (int) zDate[0]) < 1) ||
	    ((n2 = (int) zRate[0]) < 1)) {
		GtoErrMsg("%s: bad array length\n", routine);
		goto done;
	}

	/* get number of nonzero dates */
	nZero = -1;
	for (i=0; (i<=MIN(n1,n2)-1) && (zDate[i+1] > 1); i++)
		nZero = i+1;
	if (nZero < 1) {
	    GtoErrMsg("%s: not enough valid dates "
	    "(range sizes %d and %d)\n",
			routine, n1, n2);
	     goto done;
	}

	*that = GtoMakeTCurve(
			refDate[1],
			zDate+1,
			zRate+1,
			nZero,
			(double) 1L,
			GTO_ACT_365F);
	if (*that == NULL) goto done;

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	    *that = NULL;
	}
	return(status);

}


/*f--------------------------------------------------------------
 * Writes a TCurve "that" to a wrapper interface.
 * On exit, "refDate" is the array of length 1 containg thhe
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveWrapWrite(
	TCurve	*that,		/* (I) */
	TDateL	*refDate,	/* (O) (can be NULL) */
	TDateL	*zDate,		/* (O) (can be NULL) */
	FloatL	*zRate)		/* (O) (can be NULL) */
{
static	char	routine[] = "DrlTCurveWrapWrite";
	int	status = FAILURE;
	int	i, nZero;
	/* */
	nZero = that->fNumItems;

	if (refDate != NULL) {
		WRAP_CHECK_VECTOR(refDate);
		refDate[1] = that->fBaseDate;
	}

	if (zDate != NULL) {
		WRAP_CHECK_VECTOR(zDate);
		ASSERT_OR_DONE(ARGSIZE(zDate) >= nZero);
		for (i=0; i<=nZero-1; i++) {
			zDate[i+1] = that->fArray[i].fDate;
		}
	}

	if (zRate != NULL) {
		WRAP_CHECK_VECTOR(zRate);
		ASSERT_OR_DONE(ARGSIZE(zRate) >= nZero);
		for (i=0; i<=nZero-1; i++) {
			zRate[i+1] = that->fArray[i].fRate;
		}
	}

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}

#ifdef	_SKIP	/*$$$ now in drl vtype */
/*f--------------------------------------------------------------
 * Reads a TCurve "that" from a wrapper interface.
 * On entry, "refDate" is the array of length 1 containg the
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates (annual compounding, Act/365).
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveWrapRead(
	TCurve	**that,		/* (O) new curve */
	TDateL	*refDate,	/* (I) */
	TDateL	*zDate,		/* (I) */
	FloatL	*zRate)		/* (I) */
{
static	char	routine[] = "DrlTCurveWrapRead";
	int	status = FAILURE;
	int	i,
		n1,
		n2,
		nZero;


	if (((int) refDate[0] != 1) ||
	    ((n1 = (int) zDate[0]) < 1) ||
	    ((n2 = (int) zRate[0]) < 1)) {
		GtoErrMsg("%s: bad array length\n", routine);
		goto done;
	}

	/* get number of nonzero dates */
	nZero = -1;
	for (i=0; (i<=MIN(n1,n2)-1) && (zDate[i+1] > 1); i++)
		nZero = i+1;
	if (nZero < 1) {
		GtoErrMsg("%s: not enough valid dates (range sizes %d and %d)\n",
			routine, n1, n2);
		goto done;
	}

	*that = GtoMakeTCurve(
			refDate[1],
			zDate+1,
			zRate+1,
			nZero,
			(double) 1L,
			GTO_ACT_365F);
	if (*that == NULL) goto done;

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	    *that = NULL;
	}
	return(status);

}


/*f--------------------------------------------------------------
 * Writes a TCurve "that" to a wrapper interface.
 * On exit, "refDate" is the array of length 1 containg thhe
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveWrapWrite(
	TCurve	*that,		/* (I) */
	TDateL	*refDate,	/* (O) */
	TDateL	*zDate,		/* (O) */
	FloatL	*zRate)		/* (O) */
{
static	char	routine[] = "DrlTCurveWrapWrite";
	int	n1 = -1,
		n2 = -1,
		i,
		nZero,
		errCode=0;
#undef	CHECK
#define	CHECK(str)	if (errCode != 0) {GtoErrMsg(\
			"%s: failed %s (%d)\n", routine, (str), errCode);\
			 goto done;}


	/*
	 *
	 */
	nZero = DrlTCurveNumItems(that);

	errCode = ((refDate != NULL) && ((int) refDate[0] != 1) ? 1 : 0);
	CHECK("bad refDate argument");

	errCode = ((zDate != NULL) && ((n1 = (int) zDate[0]) < 1) ? 1 : 0);
	CHECK("bad zDate argument");

	errCode = ((zRate != NULL) && ((n2 = (int) zRate[0]) < 1) ? 1 : 0);
	CHECK("bad zRate argument");


	errCode = ((zDate != NULL) && (n1 < nZero) ? 1 : 0);
	CHECK("zDate array not long enough");

	errCode = ((zRate != NULL) && (n2 < nZero) ? 1 : 0);
	CHECK("zRate array not long enough");



	if (refDate != NULL) {
		refDate[1] = DrlTCurveBaseDate(that);
	}

	if (zDate != NULL) {
		for (i=0; i<=nZero-1; i++) {
			zDate[i+1] = DrlTCurveDate(that, i);
		}
	}

	if (zRate != NULL) {
		for (i=0; i<=nZero-1; i++) {
			zRate[i+1] = DrlTCurveRate(that, i);
		}
	}

done:
	return (errCode);
}

#endif	/*_SKIP*/	/*$$$ now in drl vtype */



/*f--------------------------------------------------------------
 * Forwards a TCurve "that".
 */

DLL_EXPORT(int)
DrlTCurveForwardDateTCurve(
	TCurve *that,		/* (I) zero curve */
	TDate resetDate,	/* (I) reset date */
	TCurve **outCurve)	/* (O) output zero curve */
{
static	char	routine[] = "DrlTCurveForwardTCurve";
	int	status = FAILURE;
	long	diffDays;
	int	idx;
	double	z0, z1;
	TCurve	*newCurve = NULL;
	


	if ((diffDays = resetDate - that->fBaseDate) < 0) {
	    GtoErrMsg("%s: reset date (%s) < base Date (%s).\n",
		GtoFormatDate(resetDate),
		GtoFormatDate(that->fBaseDate));
	    goto done;
	}


	if ((newCurve = GtoCopyCurve(that)) == NULL)
	    goto done;


	if (GtoDiscountDate(
		resetDate,
		that,
		GTO_LINEAR_INTERP,
		&z0) != SUCCESS)
			goto done;

	newCurve->fBaseDate = resetDate;

	for (idx=0; idx<=that->fNumItems-1; idx++) {
	    newCurve->fArray[idx].fDate += diffDays;

	    if (GtoDiscountDate(
		newCurve->fArray[idx].fDate,
		that,
		GTO_LINEAR_INTERP,
		&z1) != SUCCESS)
			goto done;

	    if (GtoDiscountToRate(
		z1/z0,
		newCurve->fBaseDate,
		newCurve->fArray[idx].fDate,
		newCurve->fDayCountConv,
		(long)newCurve->fBasis,

		&newCurve->fArray[idx].fRate) != SUCCESS)
			goto done;

	}

	*outCurve = newCurve;

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoFreeTCurve(newCurve);
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*f--------------------------------------------------------------
 * Forwards a TCurve "that".
 */

DLL_EXPORT(int)
DrlTCurveShiftRates(
	TCurve *that,		/* (B) zero curve */
	double amount)		/* (I) amount */
{
	int	idx;

	for (idx=0; idx<=that->fNumItems-1; idx++) {
	    that->fArray[idx].fRate += amount;
	}
	return(SUCCESS);
}


/*f--------------------------------------------------------------
 * Forwards a TCurve "that".
 */

DLL_EXPORT(int)
DrlTCurveShiftDates(
	TCurve *that,		/* (B) zero curve */
	TDate newBaseDate)	/* (I) new base date */
{
	int	idx;
	long	nDays = (newBaseDate - that->fBaseDate);

	that->fBaseDate = newBaseDate;

	for (idx=0; idx<=that->fNumItems-1; idx++) {
	    that->fArray[idx].fDate += nDays;
	}
	return(SUCCESS);
}


/*f--------------------------------------------------------------
 * Copies the curve <i> templateCurve</i> but replaces
 * all values by interpolating the curve <i> valueCurve</i>.\\
 * <b> Remark:</b> if in the template the first date is equal
 * to the base date, the first date is removed.
 */

DLL_EXPORT(TCurve*)
DrlTCurveInterpTCurve(
	TCurve *templateCurve,		/* (I) template */
	TCurve *valueCurve)		/* (I) to be interpolated */
{
static	char	routine[] = "DrlTCurveInterpTCurve";
	int	status = FAILURE;

	TCurve	*newCurve = NULL;
	int	idx;

	if (templateCurve->fBaseDate == templateCurve->fArray[0].fDate) {
	    newCurve = GtoNewTCurve(
		templateCurve->fBaseDate,
		templateCurve->fNumItems-1,
		templateCurve->fBasis,
		templateCurve->fDayCountConv);
	    if (newCurve == NULL) goto done;

	    for (idx=0; idx<=newCurve->fNumItems-1; idx++) {

		    newCurve->fArray[idx].fDate = 
		    	templateCurve->fArray[idx+1].fDate;

	    	    if (GtoInterpRate(
			newCurve->fArray[idx].fDate,
			valueCurve,
			GTO_LINEAR_INTERP,
			&newCurve->fArray[idx].fRate) != SUCCESS)
				goto done;
	    }



	} else {
		if ((newCurve = GtoCopyCurve(templateCurve)) == NULL)
			goto done;

		for (idx=0; idx<=templateCurve->fNumItems-1; idx++) {
	    	    if (GtoInterpRate(
			newCurve->fArray[idx].fDate,
			valueCurve,
			GTO_LINEAR_INTERP,
			&newCurve->fArray[idx].fRate) != SUCCESS)
				goto done;
		}
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoFreeTCurve(newCurve);
		newCurve = NULL;
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(newCurve);
}




/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

/*f--------------------------------------------------------------
 * Reads a TCurve "that" from a file pointer "fp"
 * as a bae volatility curve with London format.
 * The basis is set to 4 (quarterly).
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveLondonBaseVolFpRead(TCurve **that, FILE *fp)
{
static	char		routine[] = "DrlTCurveLondonBaseVolFpRead";
	char		buf[255], buf2[255];
	int		i, line = 0,
			errCode = FAILURE;
	TDate		baseDate;
	int		nPts;
	double		basis;
	TDayCount	dayCountConv;

	*that = NULL;

	/*
	 * (1) London format:
	 * # Base volatility frequency ("S" or "Q")
	 * Q
	 * # Nb of base volatilities
	 * 20
	 * # Base volatility dates and volatilities in %
	 * 19970411         7.300
	 */
	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	switch (toupper(buf[0])) {
	case 'A':
		basis = 1L; break;
	case 'S':
		basis = 2L; break;
	case 'Q':
		basis = 4L; break;
	default:
	    GtoErrMsg("%s: bad fequency `%c'\n", routine, buf[0]);
	}


	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (sscanf(buf, "%d", &nPts) != 1) 
		goto done;


	baseDate = 0L;
	dayCountConv = GTO_ACT_365F;

	/* create new TCurve */
	*that = GtoNewTCurve(baseDate, nPts, basis, dayCountConv);
	if (*that == NULL) goto done;

	if (nPts > 0) {
	    for (i=0; i<=nPts-1; i++) {
		if ((DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) ||
		    (sscanf(buf, "%s %lf", buf2, &DrlTCurveRate(*that,i)) != 2) ||
		    (DrlTDateScanYMD(buf2, &DrlTCurveDate(*that, i)) != 0)) {
			goto done;
		}
		DrlTCurveRate(*that, i) *= 1e-2;
	    }

	    /*for (i=0; i<=nPts-1; i++) 
			printf("[%3d/%3d]  %s   %lf\n", i, nPts,
				TDatePrint(NULL, dates[i]), rates[i]);*/
	}
	(*that)->fBaseDate = (*that)->fArray[0].fDate;

	/* made it through */
	errCode = SUCCESS;
done:
	if (errCode != SUCCESS) {
	    if (*that != NULL) GtoFreeTCurve(*that);
	    GtoErrMsg("%s: failed (line %d)\n", routine, line);
	}
	return(errCode);
}

/*f--------------------------------------------------------------
 * Reads a TCurve "that" from a file pointer "fp"
 * as a base smile curves with London format.
 * The basis is set to 4 (quarterly).
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveLondonBaseVolFpRead_Mod(TCurve ***that, FILE *fp, int nbCurves)
{
static  char    routine[] = "DrlTCurveLondonBaseVolFpRead_Mod";
    char        buf[255], buf2[255];
    int         i, j, line = 0,
                errCode = FAILURE;
    TDate       baseDate;
    int         nPts;
    double      basis;
    TDayCount   dayCountConv;
    TCurve      *CurrCurve = NULL;
    
    *that = NULL;

    /*
     * (1) London format:
     * # Base volatility frequency ("S" or "Q")
     * Q
     * # Nb of base volatilities
     * 20
     * # skew
     * 19970411         1.700
     */
    if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
        goto done;
    switch (toupper(buf[0])) {
    case 'A':
        basis = 1L; break;
    case 'S':
        basis = 2L; break;
    case 'Q':
        basis = 4L; break;
    default:
        GtoErrMsg("%s: bad fequency `%c'\n", routine, buf[0]);
    }


    if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
        goto done;
    if (sscanf(buf, "%d", &nPts) != 1) 
        goto done;


    baseDate = 0L;
    dayCountConv = GTO_ACT_365F;

    if( (*that = NEW_ARRAY(
                 TCurve*,
                 nbCurves)) == NULL) goto done;

    for (j = 0; j < nbCurves; j++)
    {
        /* create new TCurve */
        CurrCurve = GtoNewTCurve(baseDate, nPts, basis, dayCountConv);
        if (CurrCurve == NULL) goto done;
        *((*that) + j) = CurrCurve;

        if (nPts > 0) 
        {
            for (i=0; i<=nPts-1; i++) 
            {
                if ((DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) ||
                (sscanf(buf, "%s %lf", buf2, &DrlTCurveRate(CurrCurve,i)) != 2) ||
                (DrlTDateScanYMD(buf2, &DrlTCurveDate(CurrCurve, i)) != 0)) 
                {
                    goto done;
                }
            }
        }
        CurrCurve->fBaseDate = CurrCurve->fArray[0].fDate;
        CurrCurve = NULL;
    }
    /* made it through */
    errCode = SUCCESS;
done:
    if (errCode != SUCCESS) 
    {
        if (*that != NULL) 
        {
            for (i = 0; i < nbCurves; i++)
            {
                GtoFreeTCurve(*((*that)+i));
            }
            FREE((void*) *that);
        }
        GtoErrMsg("%s: failed (line %d)\n", routine, line);
    }
    return(errCode);
}


/*---------------------------------------------------------------
 *
 */

DLL_EXPORT(int)
DrlTCurveLondonBaseVolFileRead(TCurve **that, char *fnam)
{
	FILE		*fp;
	int		errCode ;
static	char		routine[] = "DrlTCurveLondonBaseVolFileRead";


	if ((fp = fopen(fnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open file `%s' (%s)\n",
			routine, fnam, DrlStrerror());
		errCode = -1;
		return errCode;
	} else {
		errCode = DrlTCurveLondonBaseVolFpRead(that, fp);
		if (errCode != 0)
			GtoErrMsg("%s: can't read file `%s' (code %d)\n",
				routine, fnam, errCode);
		fclose(fp);
		return(errCode);
	}
}

/*---------------------------------------------------------------
 *
 */

DLL_EXPORT(int)
DrlTCurveLondonBaseVolFileRead_Mod(TCurve ***that, char *fnam, int nbCurves)
{
	FILE		*fp;
	int		errCode ;
static	char		routine[] = "DrlTCurveLondonBaseVolFileRead_Mod";


	if ((fp = fopen(fnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open file `%s' (%s)\n",
			routine, fnam, DrlStrerror());
		errCode = -1;
		return errCode;
	} else {
		errCode = DrlTCurveLondonBaseVolFpRead_Mod(that, fp, nbCurves);
		if (errCode != 0)
			GtoErrMsg("%s: can't read file `%s' (code %d)\n",
				routine, fnam, errCode);
		fclose(fp);
		return(errCode);
	}
}

/****************************************************************/
/*---------------------------------------------------------------
 * Prints a 1-line report of the par yields of the TCurve "that"
 * in the string "s". Returns "s".
 * If "s" is NULL, prints in a static string.
 */

DLL_EXPORT(char*)
DrlTCurvePrintYields(TCurve *that, char *s)
{
static	char	buf[256];
	char	buf2[32];
	double	fwdRate;

	s = (s == NULL ? buf : s);


#undef	PRINTRATE
#define	PRINTRATE(sMat, tMat, freq, dcb) {sprintf(buf2, " %s=%7.4f%%", sMat,\
	(1e2*fwdRate, DrlTCurveFwdRate(that, 0e0, tMat, freq, dcb, &fwdRate))); \
	strcat(s, buf2);}


	sprintf(s, "[%10s, %s%2d] ",
		GtoFormatDate(that->fBaseDate),
		GtoFormatDayCountConv(that->fDayCountConv),
		(int) that->fBasis);

	PRINTRATE("ON",  2e0/365e0, 0, GTO_ACT_360)
	PRINTRATE("1M",  1e0/12e0,  0, GTO_ACT_360)
	PRINTRATE("3M",  0.25,      0, GTO_ACT_360)
	PRINTRATE("6M",  0.50,      0, GTO_ACT_360)
	PRINTRATE("12M", 1.00,      0, GTO_ACT_360)

	PRINTRATE("18M",  1.50, 2, GTO_B30_360)
	PRINTRATE("2Y",   2.00, 2, GTO_B30_360)
	PRINTRATE("3Y",   3.00, 2, GTO_B30_360)
	PRINTRATE("4Y",   4.00, 2, GTO_B30_360)
	PRINTRATE("5Y",   5.00, 2, GTO_B30_360)
	PRINTRATE("6Y",   6.00, 2, GTO_B30_360)
	PRINTRATE("7Y",   7.00, 2, GTO_B30_360)
	PRINTRATE("8Y",   8.00, 2, GTO_B30_360)
	PRINTRATE("9Y",   9.00, 2, GTO_B30_360)
	PRINTRATE("10Y", 10.00, 2, GTO_B30_360)
	PRINTRATE("12Y", 12.00, 2, GTO_B30_360)
	PRINTRATE("15Y", 15.00, 2, GTO_B30_360)
	PRINTRATE("20Y", 20.00, 2, GTO_B30_360)
	PRINTRATE("25Y", 25.00, 2, GTO_B30_360)
	PRINTRATE("30Y", 30.00, 2, GTO_B30_360)


#undef	PRINTRATE

	return(s);
}


/*f--------------------------------------------------------------
 * Convenience routine to compute a forward rate
 * on a given zero curve.
 */

DLL_EXPORT(int)
DrlTCurveFwdRate(
	TCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	TDayCount dayCountConv,	/* (I) day count convention */
	double *yield)		/* (O) forward rate */
{
static	char	routine[] = "DrlTCurveFwdRate";
	int	status = FAILURE;

	TDateInterval	payInterval;
	TDate		startDate,
			maturityDate;


	if (GtoTDateAdvanceYears(that->fBaseDate, tExp, &startDate)
		!= SUCCESS) goto done;

	if (GtoTDateAdvanceYears(startDate, tMat, &maturityDate)
		!= SUCCESS) goto done;

	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (GtoFreq2TDateInterval((long) freq, &payInterval)
			!= SUCCESS) goto done;

		if (GtoZerosToCouponsPoint(that,
			GTO_LINEAR_INTERP,
			startDate,
			&payInterval,
			maturityDate,
			dayCountConv,
			GTO_STUB_BOND,
			FALSE,
			yield) != SUCCESS) goto done;
		break;

	case 0:
		if (GtoZerosToSimplePoint(that,
			GTO_LINEAR_INTERP,
			startDate,
			maturityDate,
			dayCountConv,
			yield) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad frequency.\n", routine);
		goto done;
	}

	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}

/*f--------------------------------------------------------------
 * Convenience routine to compute a forward sensitivities
 * on a given zero curve.
 */

DLL_EXPORT(int)
DrlTCurveFwdSwapsSens(
	TCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	TDayCount dayCountConv,	/* (I) day count convention */
	char *what,		/* (I) "Y"=yield, "D"=duration */
	double *retVal)		/* (O) output fwd sensitivity */
{
static	char	routine[] = "DrlTCurveFwdDuration";
	TDateInterval	payInterval;
	TDate		startDate,
			maturityDate;
	double		yield;
	int		status = FAILURE;


	if (GtoTDateAdvanceYears(that->fBaseDate, tExp, &startDate)
		!= SUCCESS) goto done;

	if (GtoTDateAdvanceYears(startDate, tMat, &maturityDate)
		!= SUCCESS) goto done;

	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (GtoFreq2TDateInterval((long) freq, &payInterval)
			!= SUCCESS) goto done;

		if (GtoZerosToCouponsPoint(that,
			GTO_LINEAR_INTERP,
			startDate,
			&payInterval,
			maturityDate,
			dayCountConv,
			GTO_STUB_BOND,
			FALSE,
			&yield) != SUCCESS) goto done;
		break;

	case 0:
		if (GtoZerosToSimplePoint(that,
			GTO_LINEAR_INTERP,
			startDate,
			maturityDate,
			dayCountConv,
			&yield) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad frequency (%d).\n", routine, freq);
		goto done;
	}


	switch (toupper(what[0])) {
	case 'Y':
		/* fwd yield */
		*retVal = yield;
		break;
	case 'D':
		/* fwd duration */
		switch (freq) {
		case 1:
		case 2:
		case 4:
		case 12:
			if (GtoBondModDuration(
				yield,
				yield,
				(long) freq,
				tMat,
				GTO_STUB_SIMPLE,
				retVal) != SUCCESS)
					goto done;
			break;
		default:
			GtoErrMsg("%s: bad frequency (%d).\n", routine, freq);
			goto done;
		}
		break;
	default:
		GtoErrMsg("%s: bad sensitivity (%c).\n", routine, what[0]);
		goto done;
	}


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Allocates and returns an array containing the dates 
 * of the curve <i> that</i>.
 */

DLL_EXPORT(int)
DrlTCurveDatesArray(
	TCurve *that,		/* (I) curve */
	int *nDates,		/* (O) # of dates */
	TDate **dates)		/* (O) array of dates */
{
static	char	routine[] = "DrlTCurveDatesArray";
	int	i;

	*dates = NULL;
	*nDates = 0;
	if (DrlTCurveNumItems(that) <= 0) {
		GtoErrMsg("%s: Empty curve.\n", routine);
		return(FAILURE);
	}

	if ((*dates = NEW_ARRAY(TDate, DrlTCurveNumItems(that))) == NULL) {
		GtoErrMsg("%s: Malloc failed.\n", routine);
		return(FAILURE);
	}
	*nDates = DrlTCurveNumItems(that);

	for (i=0; i<=*nDates-1; i++) {
		(*dates)[i] = DrlTCurveDate(that, i);
	}
	return(SUCCESS);
}




/*f--------------------------------------------------------------
 * Convenience function to compute a forward rate on a zero
 * coupon curve.
 * It calls the C Analytics Library routines
 * <i> GtoZerosToCouponsPoint</i> and
 * <i> GtoZerosToSimplePoint</i>.
 */

DLL_EXPORT(int)
DrlTCurveForwardRate(
	TCurve* that,		/* (I) zero curve */
	TDate resetDate,	/* (I) reset date */
	TDateInterval lMat,	/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	TDayCount dayCountConv,	/* (I) day count convention */
	double *yield)		/* (I) forward yield */
{
static	char	routine[] = "DrlTCurveForwardRate";
	int		status = FAILURE;
	TDateInterval	payInterval;
	TDate		maturityDate;


	if (GtoDtFwdAny(resetDate, &lMat, &maturityDate)
		!= SUCCESS) goto done;

	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (GtoFreq2TDateInterval((long) freq, &payInterval)
			!= SUCCESS) goto done;

		if (GtoZerosToCouponsPoint(that,
			GTO_LINEAR_INTERP,
			resetDate,
			&payInterval,
			maturityDate,
			dayCountConv,
			GTO_STUB_BOND,
			FALSE,
			yield) != SUCCESS) goto done;
		break;

	case 0:
		if (GtoZerosToSimplePoint(that,
			GTO_LINEAR_INTERP,
			resetDate,
			maturityDate,
			dayCountConv,
			yield) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad frequency %d \n", routine, freq);
		goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}





/*f-------------------------------------------------------------
 * Computes the convexity adjustement of a forward rate
 * on a zero curve.
 */


DLL_EXPORT(int)
DrlTCurveForwardRateCvxAdj(
	TCurve* that,		/* (I) zero curve (base date only used) */
	TDate resetDate,	/* (I) reset date */
	TDateInterval lMat,	/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	TDayCount dayCountConv,	/* (I) day count convention */
	double yield,		/* (I) forward yield (use DrlTCurveForwardRate) */
	double yldVol,		/* (I) yield volatility */
	int cvxAdjType,		/* (I) 0=std, 1=brute force */
	double *cvxAdj)		/* (I) cvx adjustment */
{
static	char	routine[] = "DrlTCurveForwardRateCvxAdj";
	int	status = FAILURE;
	TDate	maturityDate;
	long	numDays;
	long	den;
	double	tExp, tMat, dur, cvx;

	/*
	 *
	 */
#ifdef	_SKIP
	if(DrlTCurveForwardRate(
		that,
		resetDate,
		lMat,
		freq,
		dayCountConv,
		&yield) != SUCCESS) goto done;
#endif

	if (GtoDayCountFraction(that->fBaseDate, resetDate,
		GTO_ACT_365F, &tExp) != SUCCESS) goto done;
	if (GtoDtFwdAny(resetDate, &lMat, &maturityDate)
		!= SUCCESS) goto done;

	if (GtoDayCountFraction(resetDate, maturityDate,
		GTO_B30_360, &tMat) != SUCCESS) goto done;


	switch (freq) {
	case 0:
		/* Simple Rate */
		if (GtoDaysInYearFromDayCountConv(dayCountConv, &den)
			!= SUCCESS) goto done;

		if (GtoDaysDiff(resetDate, maturityDate, dayCountConv,
			&numDays) != SUCCESS) goto done;

		if (GtoMoneyMarketModDuration(
			yield,
			numDays,
			(long) den,
			&dur) != SUCCESS) goto done;

		if (GtoMoneyMarketConvexity(
			yield,
			numDays,
			(long) den,
			&cvx) != SUCCESS) goto done;
		break;
	case 1:
	case 2:
	case 4:
	case 12:
		/* Compounded Rate */
		if (GtoBondModDuration(
			yield,
			yield,
			(long) freq,
			tMat,
			GTO_STUB_BOND,
			&dur) != SUCCESS) goto done;

		if (GtoBondConvexity(
			yield,
			yield,
			(long) freq,
			tMat,
			GTO_STUB_BOND,
			&cvx) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad frequency %d\n", routine, freq);
		goto done;
	}


	switch (cvxAdjType) {
	case 0:
		if (GtoSwapConvexityAdj(
			yield,
			dur,
			cvx,
			yldVol,
			tExp,
			cvxAdj) != SUCCESS) goto done;
		break;
#ifdef	_SKIP
	case 1:
		if (DrlNumCvxAdj(
			tExp,
			tMat,
			yield,
			freq,
			yldVol,
			1e0,
			cvxAdj) != SUCCESS) goto done;
		break;
#endif
	default:
		GtoErrMsg("%s: bad adjustment type\n", routine);
		goto done;
	}




	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



/*f--------------------------------------------------------------
 * Computes the zero curve sensitivities (duration,
 * convexity, cubicity) of a bond.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlTCurveBondSens(
	TCurve *zcCurve,
	TDate expDate,		/* (I) time to reset */
	TDate matDate,		/* (I) fwd maturity */
	int freq,		/* (I) yield frequency (1,2,4,12) */
	double *dur,		/* (O) regular duration */
	double *zDur,		/* (O) zCurve duration */
	double *zCvx,		/* (O) zCurve convexity */
	double *zP)		/* (O) zCurve cubicity */
{
static	char	routine[] = "DrlTCurveBondSens";
	int		status = FAILURE;
	int		nc, i;
	TDate		refDate;
	TDateInterval	payInterval;
	double		Rb, f,
			DR, Br, Dr, Cr, Pr,
			Ze, Zm, Ym, tm, te;

	/*
	 *
	 */
	refDate = zcCurve->fBaseDate;
	if (checkFrequency(routine, freq) == FALSE) goto done;
	if (GtoFreq2TDateInterval((long) freq, &payInterval) != SUCCESS)
		goto done;

	if (GtoDayCountFraction(refDate, expDate, GTO_B30_360,
		&te) != SUCCESS) goto done;
	if (GtoDayCountFraction(expDate, matDate, GTO_B30_360,
		&tm) != SUCCESS) goto done;
	if (GtoDiscountDate(expDate, zcCurve, GTO_LINEAR_INTERP,
		&Ze) != SUCCESS) goto done;

	nc = freq * (int) tm;
	f = (double) freq;


	/* compute fwd yield */
	if (GtoZerosToCouponsPoint(zcCurve,
			GTO_LINEAR_INTERP,
			refDate,
			&payInterval,
			matDate,
			GTO_B30_360,
			GTO_STUB_BOND,
			FALSE,
			&Rb) != SUCCESS) goto done;

	/* duration of par swap */
	DR = (1e0 - pow(1e0 + Rb/f, (double) -nc)) / Rb;

	/*DrlFPrintf(stdout, "Rb=%lf DR=%lf nc=%d f=%lf te=%lf Ze=%lf\n",
		Rb, DR, nc, f, te, Ze);*/

	/* compute Z durations */
	Br = Dr = Cr = Pr = 0e0;
	matDate = expDate;
	for (i=1; i<=nc; i++)
	{
		if (GtoDtFwdAny(matDate, &payInterval, &matDate)
			!= SUCCESS) goto done;
		if (GtoDiscountDate(matDate, zcCurve, GTO_LINEAR_INTERP,
			&Zm) != SUCCESS) goto done;
		Zm /= Ze;
		if (GtoDayCountFraction(expDate, matDate, GTO_B30_360,
			&tm) != SUCCESS) goto done;

		Ym = f * (pow(Zm, -1e0/tm/f) - 1e0);

		/*printf("[%2d] tm=%lf Zm=%lf Ym=%lf\n", i, tm, Zm, Ym);*/

		Br += Zm;
		Dr += Zm * Ym * tm;
		Cr += Zm * Ym*Ym * tm*tm;
		Pr += Zm * Ym*Ym * tm*tm*tm;

		/* if maturity */
		if (i == nc)
		{
			Br = (Rb / f) * Br + Zm;
			Dr = Dr / f
			    + Zm * Ym * tm / Rb;
			Cr = Cr / (f*Rb)
			    + Zm * Ym*Ym * tm*tm / (Rb*Rb);
			Pr = Pr / (f*Rb)
			    + Zm * Ym*Ym * tm*tm*tm / (Rb*Rb);
			break;
		}
	}

	*dur  = DR;
	*zDur = Dr;
	*zCvx = Cr;
	*zP   = Pr;


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}




