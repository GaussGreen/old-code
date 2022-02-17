/****************************************************************
 * Module:	DRL
 * Submodule:	TS
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <float.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#ifdef	DRL_CLIB
# include "date_sup.h"
# include "yearfrac.h"
# include "convert.h"
# include "check.h"
# include "zr2coup.h"		/* Analytics C Library */
# include "zr2simp.h"		/* Analytics C Library */
# include "duration.h"
# include "convex.h"
# include "swapadj.h"
#endif

#include "drlio.h"
#include "drlvtype.h"	/* LIL interface */
#include "drltime.h"	/* DrlDDateScanYMD */
#include "drlstr.h"	/* for sscanf */
#include "drlproc.h"	/* DrlStrerror */
#include "drlinter.h"	/* linear interpolation */

#include "drlts.h"


/*f--------------------------------------------------------------
 * Creates and returns a new DCurve.
 */

DCurve*
DrlDCurveNew(
	DDate valueDate,		/* (I) value date */
	int numItems,			/* (I) number of zero dates */
	double basis,			/* (I) rate frequency (1,2,4, 12) */
	DDayCount dayCountConv)		/* (I) curve DCC */
{
static	char	routine[] = "DrlDCurveNew";
	int	status = FAILURE;
	DCurve	*that = NULL;

#ifdef	DRL_CLIB
	that = GtoNewTCurve(valueDate, numItems, basis, dayCountConv);
	if (that == NULL) goto done;
#else
	that = NEW(DCurve);
	if (that == NULL) goto done;
	DRLTCURVE_BASEDATE(that) = valueDate;
	DRLTCURVE_NUMITEMS(that) = numItems;
	DRLTCURVE_DCC(that) = dayCountConv;
	if (IS_ALMOST_ZERO(basis-1)) {
		that->CurveFreq = 'A';
	} else if (IS_ALMOST_ZERO(basis-2)) {
		that->CurveFreq = 'S';
	} else if (IS_ALMOST_ZERO(basis-4)) {
		that->CurveFreq = 'Q';
	} else if (IS_ALMOST_ZERO(basis-12)) {
		that->CurveFreq = 'M';
	} else {
		DrlErrMsg("%s: bad basis %lf\n", routine, (double) basis);
		goto done;
	}

#endif

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine);
	    return (NULL);
	}
	return(that);
}


/*f--------------------------------------------------------------
 * Copies a DCurve.
 */

DCurve*
DrlDCurveNewCopy(DCurve* that)
{
static	char	routine[] = "DrlDCurveNew";
	int	status = FAILURE;
	DCurve	*copy = NULL;
	int	idx;

#ifdef	DRL_CLIB
	copy = GtoCopyCurve(that);
	if (copy == NULL) goto done;
#else
	copy = NEW(DCurve);
	if (copy == NULL) goto done;

	copy->Today = that->Today;
	copy->SpotDays = that->SpotDays;
	copy->ValueDate = that->ValueDate;

	copy->CurveDCC = that->CurveDCC;
	copy->CurveFreq = that->CurveFreq;
	copy->NbZero = that->NbZero;

	for (idx=0; idx<copy->NbZero; idx++) {
		copy->ZeroDate[idx] = that->ZeroDate[idx];
		copy->Zero[idx] = that->Zero[idx];
	}

#endif
	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine);
	    return (NULL);
	}
	return(copy);
}

/*f--------------------------------------------------------------
 * Frees a DCurve allocated by DrlDCurveNew.
 */

int
DrlDCurveFree(DCurve* that)			/* (I) curve o free */
{
#ifdef	DRL_CLIB
	if (that != NULL) GtoFreeTCurve(that);
#else
	if (that != NULL) FREE(that);
#endif
	return (SUCCESS);
}


/*f--------------------------------------------------------------
 * Reads a DCurve in file "fnam". The result is
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

int
DrlDCurveFileRead(DCurve **that, char *fnam, long fmt)
{
	FILE		*fp;
	int		errCode ;
static	char		routine[] = "DrlDCurveFileRead";


	if ((fp = fopen(fnam, "r")) == NULL) {
		DrlErrMsg("%s: can't open file `%s' (%s)\n",
			routine, fnam, DrlStrerror());
		errCode = -1;
		return errCode;
	} else {
		errCode = DrlDCurveFpRead(that, fp, fmt);
		if (errCode != 0)
			DrlErrMsg("%s: can't read file `%s' (code %d)\n",
				routine, fnam, errCode);
		fclose(fp);
		return(errCode);
	}
}


/*f--------------------------------------------------------------
 * Writes a DCurve "that" in file "fnam".
 * The format of the file is the same as for <i> DrlDCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

int
DrlDCurveFileWrite(DCurve *that, char *fnam, long fmt)
{
	FILE		*fp;

	if ((fp = fopen(fnam, "w")) == NULL) {
		DrlErrMsg("DrlDCurveFileWrite: can't open file `%s' (%s)\n",
			fnam, DrlStrerror());
		return(-1);
	} else {
		DrlDCurveFpWrite(that, fp, fmt);
		fclose(fp);
		return(0);
	}
}



/*f--------------------------------------------------------------
 * Reads a DCurve "that" from a file pointer "fp".
 * The format of the file is the same as for <i> DrlDCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

int
DrlDCurveFpRead(DCurve **that, FILE *fp, long fmt)
{
static	char		routine[] = "DrlDCurveFpRead";
	char		buf[255], buf2[255];
	int		i, line = 0, c,
			errCode = FAILURE;
	DDate		baseDate;
	int		nPts, nDays;
	double		basis;
	DDayCount	dayCountConv;

	*that = NULL;


	/*
	 * Identify format
	 */
	if ((c = fgetc(fp)) == EOF) goto done;
	ungetc(c, fp);
	if ((fmt != DRL_TCURVE_FMT_WRAPZC) && (c != '#'))
		goto cdfmt;

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
	if (DrlDDateScan(buf, &baseDate) != 0) 
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
	dayCountConv = DRL_ACT_365F;

	/* create new DCurve */
	*that = DrlDCurveNew(baseDate, nPts, basis, dayCountConv);
	if (*that == NULL) goto done;

	if (nPts > 0) {
	    for (i=0; i<=nPts-1; i++) {
		if ((DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) ||
		    (sscanf(buf, "%s %lf", buf2,
					&DRLTCURVE_RATE(*that,i)) != 2) ||
		    (DrlDDateScan(buf2, &DRLTCURVE_DATE(*that, i)) != 0))
		{
			goto done;
		}
		DRLTCURVE_RATE(*that, i) *= 1e-2;
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
	if (DrlDDateScan(buf, &baseDate) != 0) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (sscanf(buf, "%d", &nPts) != 1) 
		goto done;

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
#ifdef	DRL_CLIB
	if (sscanf(buf, "%lf", &basis) != 1) 
		goto done;
#else
	if (sscanf(buf, "%s", buf2) != 1) 
		goto done;
	switch (toupper(buf2[0])) {
		case 'A': basis =  1e0; break;
		case 'S': basis =  2e0; break;
		case 'Q': basis =  4e0; break;
		case 'M': basis = 12e0; break;
		default:
			DrlErrMsg("%s: bad basis %c\n", routine, buf2[0]);
			goto done;
	}
#endif

	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (DrlDDayCountScan(buf, &dayCountConv) != SUCCESS) 
		goto done;


	/* create new DCurve */
	*that = DrlDCurveNew(baseDate, nPts, basis, dayCountConv);
	if (*that == NULL) goto done;

	if (nPts > 0) {
	    for (i=0; i<=nPts-1; i++) {
		line++;

		/* scan '#days offset   value' */
		if (DrlFScanInt(fp, &nDays) != SUCCESS)
			goto done;
		if (nDays < 0) {
			DrlErrMsg("%s: item %d has negative number of "
				"days from value date.\n", 
				routine, i+1);
			goto done;
		}

#ifdef	DRL_CLIB
		DRLTCURVE_DATE(*that, i) = baseDate + (DDate)nDays;
#else
		IF_FAILED_DONE( DrlDDateAddDays(
			baseDate,
			nDays,
			&DRLTCURVE_DATE(*that, i)));
#endif


		if (DrlFScanDouble(fp, &DRLTCURVE_RATE(*that,i)) != SUCCESS)
			goto done;

		switch (fmt) {
		case DRL_TCURVE_FMT_PERCENT:
			DRLTCURVE_RATE(*that, i) *= 1e-2;
			break;
		case DRL_TCURVE_FMT_STD:
		case DRL_TCURVE_FMT_CUR:
			break;
		default:
			DrlErrMsg("%s: invalid format.\n", routine);
			goto done;
		}
	    }
	}

	/* made it through */
	errCode = SUCCESS;
done:
	if (errCode != SUCCESS) {
	    if (*that != NULL) DrlDCurveFree(*that);
	    *that = NULL;
	    DrlErrMsg("%s: failed (line %d).\n", routine, line);
	}
	return(errCode);
}

/*f--------------------------------------------------------------
 * Writes a DCurve "that" to a file pointer "fp".
 * The format of the file is the same as for <i> DrlDCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

int
DrlDCurveFpWrite(DCurve *that, FILE *fp, long fmt)
{
static	char	routine[] = "DrlDCurveFpWrite";
	int	i, cnt, ndays;


	switch (fmt) {
	case DRL_TCURVE_FMT_PERCENT:
		DrlFPrintf(fp, "%s\n",  DrlDDatePrint(NULL,
				DRLTCURVE_BASEDATE(that)));
		DrlFPrintf(fp, "%d\n",  DRLTCURVE_NUMITEMS(that));
#ifdef	DRL_CLIB
		DrlFPrintf(fp, "%lf\n", DRLTCURVE_BASIS(that));
#else
		DrlFPrintf(fp, "%c\n", that->CurveFreq);
#endif
		DrlFPrintf(fp, "%s\n",
			DrlDDayCountPrint(NULL, DRLTCURVE_DCC(that)));

		for (i=0; i<=DRLTCURVE_NUMITEMS(that)-1; i++) {
		    ndays = DrlDDateDaysDiff(
			DRLTCURVE_BASEDATE(that),
			DRLTCURVE_DATE(that, i));

		    DrlFPrintf(fp, "\t%4d\t%8.6f    # %10s\n",
			ndays,
			DRLTCURVE_RATE(that, i) * 1e2,
			DrlDDatePrint(NULL, DRLTCURVE_DATE(that, i)));
		}
	break;
	case DRL_TCURVE_FMT_STD:
		DrlFPrintf(fp, "%s\n",  DrlDDatePrint(NULL,
				DRLTCURVE_BASEDATE(that)));
		DrlFPrintf(fp, "%d\n",  DRLTCURVE_NUMITEMS(that));
#ifdef	DRL_CLIB
		DrlFPrintf(fp, "%lf\n", DRLTCURVE_BASIS(that));
#else
		DrlFPrintf(fp, "%c\n", that->CurveFreq);
#endif
		DrlFPrintf(fp, "%s\n",
			DrlDDayCountPrint(NULL, DRLTCURVE_DCC(that)));

		for (i=0; i<=DRLTCURVE_NUMITEMS(that)-1; i++) {
		    ndays = DrlDDateDaysDiff(
			DRLTCURVE_BASEDATE(that),
			DRLTCURVE_DATE(that, i));

		    DrlFPrintf(fp, "\t%4d\t%8.6f    # %10s\n",
			ndays,
			DRLTCURVE_RATE(that, i),
			DrlDDatePrint(NULL, DRLTCURVE_DATE(that, i)));
		}
	break;
	case DRL_TCURVE_FMT_CUR:
		DrlFPrintf(fp, "%s\n",  DrlDDatePrint(NULL,
				DRLTCURVE_BASEDATE(that)));
		DrlFPrintf(fp, "%d\n",  DRLTCURVE_NUMITEMS(that));
#ifdef	DRL_CLIB
		DrlFPrintf(fp, "%lf\n", DRLTCURVE_BASIS(that));
#else
		DrlFPrintf(fp, "%c\n", that->CurveFreq);
#endif
		DrlFPrintf(fp, "%s\n",
			DrlDDayCountPrint(NULL, DRLTCURVE_DCC(that)));

		for (i=0; i<=DRLTCURVE_NUMITEMS(that)-1; i++) {
		    ndays = DrlDDateDaysDiff(
			DRLTCURVE_BASEDATE(that),
			DRLTCURVE_DATE(that, i));

		    DrlFPrintf(fp, "\t%4d\t%15s\n",
			ndays,
			DrlCurPrint(NULL, DRLTCURVE_RATE(that, i), 2));
		}
	break;

	case DRL_TCURVE_FMT_REGTEST:
		DrlFPrintf(fp, "ZC_BASE_DATE: %s ",
			DrlDDatePrint(NULL, DRLTCURVE_BASEDATE(that)));
		DrlFPrintf(fp, "\t%%ARRAY_END%% %%ARG_END%%\n");

		DrlFPrintf(fp, "ZC_DATES:\n");
		cnt = 0;
		for (i=0; i<=DRLTCURVE_NUMITEMS(that)-1; i++) {
		    DrlFPrintf(fp, "\t%s",
			DrlDDatePrint(NULL, DRLTCURVE_DATE(that, i)));
		    if (++cnt == 5) { DrlFPrintf(fp, "\n"); cnt=0;}
		}
		DrlFPrintf(fp, "\t%%ARRAY_END%% %%ARG_END%%\n");

		DrlFPrintf(fp, "ZC_RATES:\n");
		cnt = 0;
		for (i=0; i<=DRLTCURVE_NUMITEMS(that)-1; i++) {
		    DrlFPrintf(fp, "\t%8.6f",
			DRLTCURVE_RATE(that, i));
		    if (++cnt == 5) { DrlFPrintf(fp, "\n"); cnt=0;}
		}
		DrlFPrintf(fp, "\t%%ARRAY_END%% %%ARG_END%%\n");

	break;

	case DRL_TCURVE_FMT_WRAPZC:
		/*
		 * Wrapper format
		 * $$$ WARNING: assumes 360 MM and ACT for swaps
		 * since misses information.
		 */
		DrlFPrintf(fp, "# Start date\n%ld\n",
			that->ValueDate);

		DrlFPrintf(fp, "# Money Market basis (360 or 365)\n");
		DrlFPrintf(fp, "360\n");

		DrlFPrintf(fp, "# Annual or semi-annual curve"
			" (\"A\" or \"S\")\n");
		DrlFPrintf(fp, "S\n");

		DrlFPrintf(fp, "# Year basis for benchmark swaps "
			"(\"ACT\", \"365\" or \"360\")\n");
		DrlFPrintf(fp, "ACT\n");

		DrlFPrintf(fp, "# No of entries\n%d\n",
			that->NbZero);

		DrlFPrintf(fp, "#zero maturity yyyymmdd rates "
			"(ACT/365F annual)");

		for (i=0; i<that->NbZero; i++) {
			DrlFPrintf(fp, "%ld %12.8f\n",
    				that->ZeroDate[i],
    				that->Zero[i]*1e2);
		}
	break;

	case DRL_TCURVE_FMT_WRAPBV:
		/*
		 * Wrapper format for base vol matrix
		 */
		DrlFPrintf(fp, "# Base volatility frequency "
			"(\"M\" or \"S\" or \"Q\" or \"A\")\n");
		DrlFPrintf(fp, "%c\n",
    			that->CurveFreq);


		DrlFPrintf(fp, "# Nb of base volatilities\n%d\n",
			that->NbZero);

		DrlFPrintf(fp, "# Base volatility dates and volatilities"
			" in percentage\n");
		for (i=0; i<that->NbZero; i++) {
			DrlFPrintf(fp, "%ld %12.8f\n",
    				that->ZeroDate[i],
    				that->Zero[i]*1e2);
		}

	break;

	case 9999:
		DrlFPrintf(fp, "%s\n",  DrlDDatePrint(NULL,
				DRLTCURVE_BASEDATE(that)));
		DrlFPrintf(fp, "%d\n",  DRLTCURVE_NUMITEMS(that));
#ifdef	DRL_CLIB
		DrlFPrintf(fp, "%lf\n", DRLTCURVE_BASIS(that));
#else
		DrlFPrintf(fp, "%c\n", that->CurveFreq);
#endif
		DrlFPrintf(fp, "%s\n",
			DrlDDayCountPrint(NULL, DRLTCURVE_DCC(that)));

		for (i=0; i<=DRLTCURVE_NUMITEMS(that)-1; i++) {
		    DrlFPrintf(fp, "\t%10s\t%15.8f\n",
			DrlDDatePrint(NULL, DRLTCURVE_DATE(that, i)),
			DRLTCURVE_RATE(that, i));
		}
	break;

	default:
		DrlErrMsg("%s: invalid format.\n", routine);
		return(FAILURE);
	}

	return(SUCCESS);
done:
	DrlErrMsg("%s: failed.\n", routine);
	return(FAILURE);
}



/*f--------------------------------------------------------------
 * Reads a DCurve "that" from a wrapper interface.
 * On entry, "refDate" is the array of length 1 containg the
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates (annual compounding, Act/365).
 * Returns 0 iff OK.
 */

int
DrlDCurveWrapRead(
	DCurve	**that,		/* (O) new curve */
	long	*refDate,	/* (I) */
	long	*zDate,		/* (I) */
	double	*zRate)		/* (I) */
{
static	char	routine[] = "DrlDCurveWrapRead";
	int	status = FAILURE;
	int	i, n1, n2,
		nZero;
	DDate	baseDate;


	if (((int) refDate[0] != 1) ||
	    ((n1 = (int) zDate[0]) < 1) ||
	    ((n2 = (int) zRate[0]) < 1)) {
		DrlErrMsg("%s: bad array length\n", routine);
		goto done;
	}

	/* get number of nonzero dates */
	nZero = -1;
	for (i=0; (i<=MIN(n1,n2)-1) && (zDate[i+1] > 1); i++)
		nZero = i+1;
	if (nZero < 1) {
	    DrlErrMsg("%s: not enough valid dates "
	    "(range sizes %d and %d)\n",
			routine, n1, n2);
	     goto done;
	}


#ifdef	DRL_CLIB
	*that = GtoMakeTCurve(
			refDate[1],
			zDate+1,
			zRate+1,
			nZero,
			(double) 1L,
			DRL_ACT_365F);
	if (*that == NULL) goto done;
#else
	IF_FAILED_DONE( DrlAlibDateToDrDate(
		refDate[1],
		&baseDate));

	*that = DrlDCurveNew(
			baseDate,
			nZero,
			(double) 1L,
			DRL_ACT_365F);
	if (*that == NULL) goto done;
	for (i=0; i<nZero; i++) {
		IF_FAILED_DONE( DrlAlibDateToDrDate(
			zDate[i+1],
			&DRLTCURVE_DATE(*that, i)));
		DRLTCURVE_RATE(*that, i) = zRate[i+1];
	}
#endif


	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed\n", routine);
	    *that = NULL;
	}
	return(status);

}


/*f--------------------------------------------------------------
 * Writes a DCurve "that" to a wrapper interface.
 * On exit, "refDate" is the array of length 1 containg thhe
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates.
 * Returns 0 iff OK.
 */

int
DrlDCurveWrapWrite(
	DCurve	*that,		/* (I) */
	long	*refDate,	/* (O) (can be NULL) */
	long	*zDate,		/* (O) (can be NULL) */
	double	*zRate)		/* (O) (can be NULL) */
{
static	char	routine[] = "DrlDCurveWrapWrite";
	int	status = FAILURE;
	int	i, nZero;
	/* */
	nZero = DRLTCURVE_NUMITEMS(that);

	if (refDate != NULL) {
		DRL_WRAP_CHECK_VECTOR(refDate);
#ifdef	DRL_CLIB
		refDate[1] = DRLTCURVE_BASEDATE(that);
#else
		NOT_IMPLEMENTED;
#endif
	}

	if (zDate != NULL) {
		DRL_WRAP_CHECK_VECTOR(zDate);
		ASSERT_OR_DONE(DRL_ARGSIZE(zDate) >= nZero);
		for (i=0; i<=nZero-1; i++) {
#ifdef	DRL_CLIB
			zDate[i+1] = DRLTCURVE_DATE(that, i);
#else
			NOT_IMPLEMENTED;
#endif
		}
	}

	if (zRate != NULL) {
		DRL_WRAP_CHECK_VECTOR(zRate);
		ASSERT_OR_DONE(DRL_ARGSIZE(zRate) >= nZero);
		for (i=0; i<=nZero-1; i++) {
			zRate[i+1] = DRLTCURVE_RATE(that, i);
		}
	}

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}

#ifdef	_SKIP	/*$$$ now in drl vtype */
/*f--------------------------------------------------------------
 * Reads a DCurve "that" from a wrapper interface.
 * On entry, "refDate" is the array of length 1 containg the
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates (annual compounding, Act/365).
 * Returns 0 iff OK.
 */

int
DrlDCurveWrapRead(
	DCurve	**that,		/* (O) new curve */
	DDateL	*refDate,	/* (I) */
	DDateL	*zDate,		/* (I) */
	FloatL	*zRate)		/* (I) */
{
static	char	routine[] = "DrlDCurveWrapRead";
	int	status = FAILURE;
	int	i,
		n1,
		n2,
		nZero;


	if (((int) refDate[0] != 1) ||
	    ((n1 = (int) zDate[0]) < 1) ||
	    ((n2 = (int) zRate[0]) < 1)) {
		DrlErrMsg("%s: bad array length\n", routine);
		goto done;
	}

	/* get number of nonzero dates */
	nZero = -1;
	for (i=0; (i<=MIN(n1,n2)-1) && (zDate[i+1] > 1); i++)
		nZero = i+1;
	if (nZero < 1) {
		DrlErrMsg("%s: not enough valid dates "
			"(range sizes %d and %d)\n",
			routine, n1, n2);
		goto done;
	}

	*that = GtoMakeDCurve(
			refDate[1],
			zDate+1,
			zRate+1,
			nZero,
			(double) 1L,
			DRL_ACT_365F);
	if (*that == NULL) goto done;

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed\n", routine);
	    *that = NULL;
	}
	return(status);

}


/*f--------------------------------------------------------------
 * Writes a DCurve "that" to a wrapper interface.
 * On exit, "refDate" is the array of length 1 containg thhe
 * value date, and "zDate" and "zRate" are arrays of respectively
 * the dates and zero rates.
 * Returns 0 iff OK.
 */

int
DrlDCurveWrapWrite(
	DCurve	*that,		/* (I) */
	DDateL	*refDate,	/* (O) */
	DDateL	*zDate,		/* (O) */
	FloatL	*zRate)		/* (O) */
{
static	char	routine[] = "DrlDCurveWrapWrite";
	int	n1 = -1,
		n2 = -1,
		i,
		nZero,
		errCode=0;
#undef	CHECK
#define	CHECK(str)	if (errCode != 0) {DrlErrMsg(\
			"%s: failed %s (%d)\n", routine, (str), errCode);\
			 goto done;}


	/*
	 *
	 */
	nZero = DRLTCURVE_NUMITEMS(that);

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
		refDate[1] = DRLTCURVE_BASEDATE(that);
	}

	if (zDate != NULL) {
		for (i=0; i<=nZero-1; i++) {
			zDate[i+1] = DRLTCURVE_DATE(that, i);
		}
	}

	if (zRate != NULL) {
		for (i=0; i<=nZero-1; i++) {
			zRate[i+1] = DRLTCURVE_RATE(that, i);
		}
	}

done:
	return (errCode);
}

#endif	/*_SKIP*/	/*$$$ now in drl vtype */




/*f--------------------------------------------------------------
 * Reads a DCurve "that" from a file pointer "fp"
 * as a bae volatility curve with London format.
 * The basis is set to 4 (quarterly).
 * Returns 0 iff OK.
 */

int
DrlDCurveLondonBaseVolFpRead(DCurve **that, DDate baseDate, FILE *fp)
{
static	char		routine[] = "DrlDCurveLondonBaseVolFpRead";
	char		buf[255], buf2[255];
	int		i, line = 0,
			errCode = FAILURE;
	int		nPts;
	double		basis;
	DDayCount	dayCountConv;

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
	    DrlErrMsg("%s: bad fequency `%c'\n", routine, buf[0]);
	}


	if (DrlFGetLine(buf, sizeof(buf),fp, &line) == NULL) 
		goto done;
	if (sscanf(buf, "%d", &nPts) != 1) 
		goto done;


	dayCountConv = DRL_ACT_365F;

	/* create new DCurve */
	*that = DrlDCurveNew(baseDate, nPts, basis, dayCountConv);
	if (*that == NULL) goto done;

	if (nPts > 0) {
	    for (i=0; i<=nPts-1; i++) {
		if ((DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) ||
		    (sscanf(buf, "%s %lf", buf2, &DRLTCURVE_RATE(*that,i)) != 2) ||
		    (DrlDDateScanYMD(buf2, &DRLTCURVE_DATE(*that, i)) != 0)) {
			goto done;
		}
		DRLTCURVE_RATE(*that, i) *= 1e-2;
	    }

	}

	/* made it through */
	errCode = SUCCESS;
done:
	if (errCode != SUCCESS) {
	    if (*that != NULL) DrlDCurveFree(*that);
	    DrlErrMsg("%s: failed (line %d)\n", routine, line);
	}
	return(errCode);
}

/*---------------------------------------------------------------
 *
 */

int
DrlDCurveLondonBaseVolFileRead(DCurve **that, DDate baseDate, char *fnam)
{
	FILE		*fp;
	int		errCode ;
static	char		routine[] = "DrlDCurveLondonBaseVolFileRead";


	if ((fp = fopen(fnam, "r")) == NULL) {
		DrlErrMsg("%s: can't open file `%s' (%s)\n",
			routine, fnam, DrlStrerror());
		errCode = -1;
		return errCode;
	} else {
		errCode = DrlDCurveLondonBaseVolFpRead(that, baseDate, fp);
		if (errCode != 0)
			DrlErrMsg("%s: can't read file `%s' (code %d)\n",
				routine, fnam, errCode);
		fclose(fp);
		return(errCode);
	}
}


/****************************************************************/
/*---------------------------------------------------------------
 * Prints a 1-line report of the par yields of the DCurve "that"
 * in the string "s". Returns "s".
 * If "s" is NULL, prints in a static string.
 */

char*
DrlDCurvePrintYields(DCurve *that, char *s)
{
static	char	buf[256];
	char	buf2[32];
	double	fwdRate;

	s = (s == NULL ? buf : s);


#undef	PRINTRATE
#define	PRINTRATE(sMat, tMat, freq, dcb) { \
	sprintf(buf2, " %s=%7.4f%%", \
		sMat,\
		(DrlDCurveForwardRate3(that,0e0,tMat,freq,dcb, &fwdRate), \
		1e2*fwdRate));  \
		strcat(s, buf2);}


	sprintf(s, "[%10s, %s%] ",
		DrlDDatePrint(NULL, DRLTCURVE_BASEDATE(that)),
		DrlDDayCountPrint(NULL, DRLTCURVE_DCC(that)));

	PRINTRATE("ON",  2e0/365e0, 0, DRL_ACT_360)
	PRINTRATE("1M",  1e0/12e0,  0, DRL_ACT_360)
	PRINTRATE("3M",  0.25,      0, DRL_ACT_360)
	PRINTRATE("6M",  0.50,      0, DRL_ACT_360)
	PRINTRATE("12M", 1.00,      0, DRL_ACT_360)

	PRINTRATE("18M",  1.50, 2, DRL_B30_360)
	PRINTRATE("2Y",   2.00, 2, DRL_B30_360)
	PRINTRATE("3Y",   3.00, 2, DRL_B30_360)
	PRINTRATE("4Y",   4.00, 2, DRL_B30_360)
	PRINTRATE("5Y",   5.00, 2, DRL_B30_360)
	PRINTRATE("6Y",   6.00, 2, DRL_B30_360)
	PRINTRATE("7Y",   7.00, 2, DRL_B30_360)
	PRINTRATE("8Y",   8.00, 2, DRL_B30_360)
	PRINTRATE("9Y",   9.00, 2, DRL_B30_360)
	PRINTRATE("10Y", 10.00, 2, DRL_B30_360)
	PRINTRATE("12Y", 12.00, 2, DRL_B30_360)
	PRINTRATE("15Y", 15.00, 2, DRL_B30_360)
	PRINTRATE("20Y", 20.00, 2, DRL_B30_360)
	PRINTRATE("25Y", 25.00, 2, DRL_B30_360)
	PRINTRATE("30Y", 30.00, 2, DRL_B30_360)


#undef	PRINTRATE

	return(s);
}

#ifdef	_SKIP
/*---------------------------------------------------------------
 * Convenience routine to compute a forward sensitivities
 * on a given zero curve.
 */

int
DrlDCurveFwdSwapsSens(
	DCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	char *what,		/* (I) "Y"=yield, "D"=duration */
	double *retVal)		/* (O) output fwd sensitivity */
{
static	char	routine[] = "DrlDCurveFwdDuration";
	DInterval	payInterval;
	DDate		startDate,
			maturityDate;
	double		yield;
	int		status = FAILURE;


	if (DrlDDateAdvanceYears(DRLTCURVE_BASEDATE(that), tExp, &startDate)
		!= SUCCESS) goto done;

	if (DrlDDateAdvanceYears(startDate, tMat, &maturityDate)
		!= SUCCESS) goto done;

	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (DrlFreqToDInterval((long) freq, &payInterval)
			!= SUCCESS) goto done;

#ifdef	DRL_CLIB
		if (GtoZerosToCouponsPoint(that,
			GTO_LINEAR_INTERP,
			startDate,
			&payInterval,
			maturityDate,
			dayCountConv,
			GTO_STUB_BOND,
			FALSE,
			&yield) != SUCCESS) goto done;
#else
		PROGRAM_BUG();
#endif
		break;

	case 0:
#ifdef	DRL_CLIB
		if (GtoZerosToSimplePoint(that,
			GTO_LINEAR_INTERP,
			startDate,
			maturityDate,
			dayCountConv,
			&yield) != SUCCESS) goto done;
#else
		PROGRAM_BUG();
#endif
		break;
	default:
		DrlErrMsg("%s: bad frequency (%d).\n", routine, freq);
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
#ifdef	DRL_CLIB
			if (GtoBondModDuration(
				yield,
				yield,
				(long) freq,
				tMat,
				GTO_STUB_SIMPLE,
				retVal) != SUCCESS)
					goto done;
#else
		PROGRAM_BUG();
#endif
			break;
		default:
			DrlErrMsg("%s: bad frequency (%d).\n", routine, freq);
			goto done;
		}
		break;
	default:
		DrlErrMsg("%s: bad sensitivity (%c).\n", routine, what[0]);
		goto done;
	}


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}
#endif


/*f--------------------------------------------------------------
 * Allocates and returns an array containing the dates 
 * of the curve <i> that</i>.
 */

int
DrlDCurveDatesArray(
	DCurve *that,		/* (I) curve */
	int *nDates,		/* (O) # of dates */
	DDate **dates)		/* (O) array of dates */
{
static	char	routine[] = "DRLTCURVE_DATEsArray";
	int	i;

	*dates = NULL;
	*nDates = 0;
	if (DRLTCURVE_NUMITEMS(that) <= 0) {
		DrlErrMsg("%s: Empty curve.\n", routine);
		return(FAILURE);
	}

	if ((*dates = NEW_ARRAY(DDate, DRLTCURVE_NUMITEMS(that))) == NULL) {
		DrlErrMsg("%s: Malloc failed.\n", routine);
		return(FAILURE);
	}
	*nDates = DRLTCURVE_NUMITEMS(that);

	for (i=0; i<=*nDates-1; i++) {
		(*dates)[i] = DRLTCURVE_DATE(that, i);
	}
	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * Returns the curve frequency as an integer (0=cc, 1, 2, 4, 12).
 */

int DrlDCurveFreq(
	DCurve *that,		/* (I) zero curve */
	int *freq)		/* (O) basis (0, 1, 2, 4, 12) */
{
static	char	routine[] = "DrlDCurveFreq";

#if !defined(DRL_CLIB)
	switch (toupper(that->CurveFreq)) {
	case 'A': *freq =  1; break;
	case 'S': *freq =  2; break;
	case 'Q': *freq =  4; break;
	case 'M': *freq = 12; break;
	default:
		DrlErrMsg("%s: bad freq `%s'.\n", routine,
			that->CurveFreq);
		return (FAILURE);
	}
	return(SUCCESS);
#else
	switch ((int) (that->fBasis)) {
	case 0: *freq =  0; break;
	case 1: *freq =  1; break;
	case 2: *freq =  2; break;
	case 4: *freq =  4; break;
	case 12: *freq = 12; break;
	default:
		DrlErrMsg("%s: bad freq `%lf'.\n", routine,
			that->fBasis);
		return (FAILURE);
	}
	return(SUCCESS);
#endif
}

/*f--------------------------------------------------------------
 * Sets the curve frequency as an integer (0=cc, 1, 2, 4, 12).
 */

int DrlDCurveFreqSet(
	DCurve *that,		/* (O) zero curve */
	int freq)		/* (I) basis (0, 1, 2, 4, 12) */
{
static	char	routine[] = "DrlDCurveFreq";

#if !defined(DRL_CLIB)
	switch (freq) {
	case 1:
		that->CurveFreq = 'A'; break;
	case 2:
		that->CurveFreq = 'S'; break;
	case 4:
		that->CurveFreq = 'Q'; break;
	case 12:
		that->CurveFreq = 'M'; break;
	default:
		DrlErrMsg("%s: bad freq `%d'.\n", routine, freq);
		return (FAILURE);
	}
	return(SUCCESS);
#else
	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		that->fBasis = double(freq);
		break;
	default:
		DrlErrMsg("%s: bad freq `%d'.\n", routine, freq);
		return (FAILURE);
	}
	return(SUCCESS);
#endif
}

/*f--------------------------------------------------------------
 * Interpolates a rate off from a curve.
 */

int
DrlDCurveInterp(
	DCurve* that,		/* (I) zero curve */
	DDate matDate,		/* (I) reset date */
	double *discRate)	/* (O) forward yield */
{
static	char	routine[] = "DrlDCurveDiscFact";
	int	status = FAILURE;

#ifdef	DRL_CLIB
	IF_FAILED_DONE( GtoInterpRate(
		matDate,
		that,
		GTO_LINEAR_INTERP,
		discRate));
#else
	if (that->ValueDate > matDate) {
		DrlErrMsg("%s: value date (%s) > maturity date %s\n",
			routine,
			DrlDDatePrint(NULL, that->ValueDate),
			DrlDDatePrint(NULL, matDate));
		return (FAILURE);
	}
	IF_FAILED_DONE( DrlDDateLinearInterp1d(
    		that->ZeroDate,
    		that->Zero,
    		that->NbZero,
		matDate,
		discRate));
#endif

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Computes a discount factor off a zero curve
 */

int
DrlDCurveDiscFact(
	DCurve* that,		/* (I) zero curve */
	DDate matDate,		/* (I) reset date */
	double *discFact)	/* (O) forward yield */
{
static	char	routine[] = "DrlDCurveDiscFact";
	int	status = FAILURE;
	double	zeroYield,
		basis = 1e0,
		yearFract;

#ifdef	DRL_CLIB
	IF_FAILED_DONE( GtoDiscountDate(
		matDate,
		that,
		GTO_LINEAR_INTERP,
		discFact));
#else
	if (that->ValueDate > matDate) {
		DrlErrMsg("%s: value date (%s) > maturity date %s\n",
			routine,
			DrlDDatePrint(NULL, that->ValueDate),
			DrlDDatePrint(NULL, matDate));
		return (FAILURE);
	}
	if (that->ValueDate ==  matDate) {
		*discFact = 1e0;
		return (SUCCESS);
	}
	IF_FAILED_DONE( DrlDayCountFract(
		that->ValueDate,
		matDate,
		that->CurveDCC,
		&yearFract));

	IF_FAILED_DONE( DrlDDateLinearInterp1d(
    		that->ZeroDate,
    		that->Zero,
    		that->NbZero,
		matDate,
		&zeroYield));


	*discFact = pow(1e0+zeroYield / basis , -yearFract * basis);

	/*printf("Z: E=%s YRFRACT=%lf ZYLD=%lf ZPR=%lf\n",
		DrlDDatePrint(NULL, matDate),
		yearFract, zeroYield, *discFact);*/
#endif
	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * DCurve : compute a froward rate
 * <br><br>
 * Convenience function to compute a forward rate on a zero
 * coupon curve (using BOND stub).
 */

int
DrlDCurveForwardRate(
	DCurve* that,		/* (I) zero curve */
	DDate startDate,	/* (I) start date */
	DDate maturityDate,	/* (I) maturity date */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	double *yield)		/* (I) forward yield */
{
static	char	routine[] = "DrlDCurveForwardRate";
	int		status = FAILURE;
#if !defined(DRL_CLIB)

	double	discStart ,
		discEnd,
		discFact,
		annuity = 0e0,
		dcf;
	DDate	cfDate,
		cfDateNext;
	int	numPer;
	DInterval	payInterval,
			offsetInterval;


	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (DrlFreqToDInterval((long) freq, &payInterval)
			!= SUCCESS) goto done;
		offsetInterval = payInterval;

		/* Compute the annuity factor */
		cfDateNext =  maturityDate;
		numPer = 0;
		annuity = 0e0;
		do {
			numPer--;
			offsetInterval.prd = numPer * payInterval.prd;

			IF_FAILED_DONE( DrlDDateFwdAny(
				maturityDate, &offsetInterval, &cfDate));

			IF_FAILED_DONE( DrlDayCountFract(
				cfDate, cfDateNext, dayCountConv, &dcf));

			IF_FAILED_DONE( DrlDCurveDiscFact(
				that,
				cfDateNext,
				&discFact));
			if (numPer == 0) {
				discEnd = discFact;
			}

			annuity += discFact * dcf;

			/*printf("\t%4s %10s->%10s  DCF=%lf Z=%lf\n",
				DrlDIntervalPrint(NULL, &offsetInterval),
				DrlDDatePrint(NULL, cfDate),
				DrlDDatePrint(NULL, cfDateNext),
				dcf,
				discFact);*/

			cfDateNext = cfDate;

		} while (cfDate >  startDate);

		/* compute the stub */
		IF_FAILED_DONE( DrlDCurveDiscFact(
			that,
			startDate,
			&discStart));

		IF_FAILED_DONE( DrlDayCountFract(
			cfDate, startDate, dayCountConv, &dcf));
		annuity -= dcf * discStart;


		IF_FAILED_DONE( DrlDCurveDiscFact(
			that,
			maturityDate,
			&discEnd));

		*yield = (discStart - discEnd) / annuity;

		break;
	case 0:
		IF_FAILED_DONE( DrlDayCountFract(
			startDate, maturityDate, dayCountConv, &dcf));

		IF_FAILED_DONE( DrlDCurveDiscFact(
			that,
			startDate,
			&discStart));

		IF_FAILED_DONE( DrlDCurveDiscFact(
			that,
			maturityDate,
			&discEnd));

		*yield = (discStart / discEnd - 1e0) / dcf;

		break;
	default:
		DrlErrMsg("%s: bad frequency %d.\n", routine, freq);
		goto done;
	}

	/*printf("YLD: S=%s E=%s F=%2d DCC=%10s Y=%12.8f %% A=%lf\n",
		DrlDDatePrint(NULL, startDate),
		DrlDDatePrint(NULL, maturityDate),
		freq,
		DrlDDayCountPrint(NULL, dayCountConv),
		(*yield)*1e2, annuity);*/
#else
	DInterval	payInterval;

	switch (freq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (DrlFreqToDInterval((long) freq, &payInterval)
			!= SUCCESS) goto done;
		if (GtoZerosToCouponsPoint(
			that,
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
		DrlErrMsg("%s: bad frequency.\n", routine);
		goto done;
	}

	/*printf("YLD: S=%s E=%s F=%2d DCC=%10s Y=%12.8f %%\n",
		DrlDDatePrint(NULL, startDate),
		DrlDDatePrint(NULL, maturityDate),
		freq,
		DrlDDayCountPrint(NULL, dayCountConv),
		(*yield)*1e2);*/
#endif

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Convenience function to compute a forward rate on a zero
 * coupon curve.
 */

int
DrlDCurveForwardRate2(
	DCurve* that,		/* (I) zero curve */
	DDate startDate,	/* (I) reset date */
	DInterval maturity,	/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	double *yield)		/* (I) forward yield */
{
static	char	routine[] = "DrlDCurveForwardRate2";
	int	status = FAILURE;
	DDate	maturityDate;

	IF_FAILED_DONE( DrlDDateFwdAny(
		startDate, &maturity, &maturityDate));

	IF_FAILED_DONE( DrlDCurveForwardRate(
		that,
		startDate,
		maturityDate,
		freq,
		dayCountConv,
		yield));

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Convenience routine to compute a forward rate
 * on a given zero curve.
 */

int
DrlDCurveForwardRate3(
	DCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	double *yield)		/* (O) forward rate */
{
static	char	routine[] = "DrlDCurveForwardRate3";
	int	status = FAILURE;

	DInterval	payInterval;
	DDate		startDate,
			maturityDate;

	IF_FAILED_DONE( DrlDDateAdvanceYears(
		DRLTCURVE_BASEDATE(that), tExp, &startDate));

	IF_FAILED_DONE( DrlDDateAdvanceYears(
		startDate, tMat, &maturityDate));

	IF_FAILED_DONE ( DrlDCurveForwardRate(
		that,
		startDate,
		maturityDate,
		freq,
		dayCountConv,
		yield));

	status = SUCCESS;
done:
	if (status != SUCCESS)
		DrlErrMsg("%s: failed.\n", routine);
	return(status);
}

