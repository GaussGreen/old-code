/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	smdat.c
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 ****************************************************************/
#define	_pg_SRC
#include "dritkwrp.h"

#include <ctype.h>
#include <string.h>

#include "macros.h"
#include "date_sup.h"
#include "convert.h"
#include "zr2coup.h"		/* Analytics C Library */
#include "zr2simp.h"		/* Analytics C Library */
#include "zr2sv.h"		/* Analytics C Library */

#include "drlts.h"
#include "drlio.h"
#include "drlstr.h"		/* DrlCurPrint */
#include "drltime.h"		/* DrlTDatePrint */


#undef	SQR
#define	SQR(x)	((x)*(x))


static	int	forwardATMRate(
	TCurve *zcCurve,	/* (I) zero curve */
	TDate startDate,	/* (I) reset date */
	double rateMat,		/* (I) rate maturity */
	int rateFreq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	long rateDayCountConv,	/* (I) day count convention */
	double *atmRate);	/* (I) forward atm rate */

#define	__DEBUG__
#undef	__DEBUG__

/*f--------------------------------------------------------------
 * Computes the 3-points smile adjustment factor from the
 * coefficients $V_{-}$, $V_{+}$, $h$, $q$ and $\alpha$.
 * The {\tt smileType} determines the type of smile function used:
 * \begin{itemize}
 * \item if {\tt smileType} is 0, no smile function is computed.
 * \item if {\tt smileType} is 1, a 2-points smile power function is
 * computed using $q$ and $\alpha$.
 * \item if {\tt smileType} is 2, a 3-points smile power/parabolic function is
 * computed using $V_{+}$, $h$, $q$ and $\alpha$.
 * \end{itemize}
 * Returns SUCCESS/FAILURE.
 */

static	int
SmileFunc3Pts(
	double f,		/* (I) atm rate */
	double k,		/* (I) strike */
	double v_hi,		/* (I) V+ coefficient */
	double h,		/* (I) h coefficient */
	double q,		/* (I) q coefficient */
	double alpha,		/* (I) alpha coefficient (>=0) */
	int smiltype,		/* (I) 0=no smile, 1=2-pts, 2=3-points */
	double *smadj)		/* (O) interp smile adjustment factor */
{
static	char	routine[] = "SmileFunc3Pts";
	int	status = FAILURE;

	double	u_lo, u_hi,
		u, w,
		powf,
		parf;

	double v_lo;		/* (I) V- coefficient */

	switch (smiltype) {
	case 0:
		/* no smile */
		*smadj = 1e0;
		break;

	case 1:
		/* 2-points (no parabolic) */
		u = log(k/f);

		/* power function */
		if (IS_ALMOST_ZERO(alpha)) {
			powf = exp(-q * u);
		} else {
			powf = exp(-q * u / (1e0 - alpha*q*u));
		}

		/* smile function */
		*smadj = powf;
		break;

	case 2:
		/* 3-points (with parabolic) */
		u_hi =  h;
		u_lo = -h;
		u = log(k/f);

		/* power function */
		if (IS_ALMOST_ZERO(alpha)) {
			v_lo = exp(q*u_hi);
			powf = pow(k / f, -q);
		} else {
			v_lo = exp(- q*u_lo / (1e0 - q*u_lo*alpha));
			powf = exp(-q * u / (1e0 - alpha*q*u));
		}


		/* parabolic function */
		parf = 1e0
			+ (v_hi - 1e0) * u * (u+h) * 0.5e0 / SQR(h)
			+ (v_lo - 1e0) * u * (u-h) * 0.5e0 / SQR(h);

		/* weight */
		if (u <= u_lo) {
	    		w = 1e0;
		} else if (u >= u_hi) {
	    		w = 0e0;
		} else {
	    		w = ((h-u)/h) * ((h-u)/h) * ((h-u)/h) *
			(0.1875e0 * SQR(u/h) +
			0.5625e0 * u/h + 0.5e0);
		}

		/* smile function */
		*smadj = w * powf + (1e0-w) * parf;

		break;
	default:
		GtoErrMsg("%s: bad smile type %d.\n",
			routine, smiltype);
		goto done;
	}


	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*f---------------------------------------------------------------------
 * Default constructor (smile is off).
 */

TSmile3Data*
DrlTSmile3DataNewEmpty(void)
{
	TSmile3Data	*that = NULL;

	if ((that = NEW(TSmile3Data)) == NULL)
		return(NULL);

	that->fSmileType = 0;		/* Default is smile off */
	that->fQcoeff   = NULL;
	that->fVhicoeff = NULL;
	that->fHcoeff = 0e0;
	that->fAlpha = 0e0;
	that->fBackBoneQ = NULL;
	that->fAtmRate = NULL;

	return(that);
}

/*f---------------------------------------------------------------------
 * Destructor.
 */

void
DrlTSmile3DataDelete(TSmile3Data *that)
{
	if (that) {
		DrlTSmile3DataCleanup(that);
		FREE(that);
	}
}

/*f---------------------------------------------------------------------
 * Private function that frees all internally allocated
 * memory.
 */

int
DrlTSmile3DataCleanup(TSmile3Data *that)
{
	DrlTSwaptionMatrix2DFree(that->fQcoeff);
	DrlTSwaptionMatrix2DFree(that->fVhicoeff);
	DrlTSwaptionMatrix2DFree(that->fBackBoneQ);
	DrlTSwaptionMatrix2DFree(that->fAtmRate);

	that->fSmileType = 0;
	that->fQcoeff   = NULL;
	that->fVhicoeff = NULL;
	that->fHcoeff = 0e0;
	that->fAlpha = 0e0;
	that->fBackBoneQ = NULL;
	that->fAtmRate = NULL;

	return(SUCCESS);
}

#ifdef	_SKIP
/*f---------------------------------------------------------------------
 * Read class from file pointer {\tt fp}.
 */

void
TSmile3Data::FpRead(FILE *fp)
{
static	char	routine[] = "TSmile3Data::FpRead";

    try {
	/* cleanup */
	Cleanup();

	FP_ADVANCE("SMILE_TYPE:");
	IF_FAILED_DONE(DrlFScanVType(fp, DRL_INT_T, &fSmileType));

	FP_ADVANCE("Q_COEFF:");
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpRead(&fQcoeff, fp,
		TSWAPTION_MATRIX_FMT_NUM));

	FP_ADVANCE("VHI_COEFF:");
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpRead(&fVhicoeff, fp,
		TSWAPTION_MATRIX_FMT_NUM));

	FP_ADVANCE("H_COEFF:");
	IF_FAILED_DONE(DrlFScanVType(fp, DRL_DOUBLE_T, &fHcoeff));

	FP_ADVANCE("ALPHA_COEFF:");
	IF_FAILED_DONE(DrlFScanVType(fp, DRL_DOUBLE_T, &fAlpha));

	FP_ADVANCE("BACKBONE_Q_COEFF:");
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpRead(&fBackBoneQ, fp,
			TSWAPTION_MATRIX_FMT_NUM));

	FP_ADVANCE("ATM_RATE:");
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpRead(&fAtmRate, fp,
			TSWAPTION_MATRIX_FMT_NUM));

	return;
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}

#endif	/*_SKIP*/

/*f---------------------------------------------------------------------
 * Write class to file pointer {\tt fp}.
 */

int
DrlTSmile3DataFpWrite(TSmile3Data *that, FILE *fp)
{
	DrlFPrintf(fp, "SMILE_TYPE:\n");
	DrlFPrintf(fp, "%d\n", that->fSmileType);

	DrlFPrintf(fp, "Q_COEFF:\n");
	DrlTSwaptionMatrix2DFpWrite(that->fQcoeff,
		fp, TSWAPTION_MATRIX_FMT_NUM);

	DrlFPrintf(fp, "VHI_COEFF:\n");
	DrlTSwaptionMatrix2DFpWrite(that->fVhicoeff,
		fp, TSWAPTION_MATRIX_FMT_NUM);

	DrlFPrintf(fp, "H_COEFF:\n");
	DrlFPrintVType(fp, DRL_DOUBLE_T, &that->fHcoeff);
	DrlFPrintf(fp, "\n");

	DrlFPrintf(fp, "ALPHA_COEFF:\n");
	DrlFPrintVType(fp, DRL_DOUBLE_T, &that->fAlpha);
	DrlFPrintf(fp, "\n");

	DrlFPrintf(fp, "BACKBONE_Q_COEFF:\n");
	DrlTSwaptionMatrix2DFpWrite(that->fBackBoneQ,
		fp, TSWAPTION_MATRIX_FMT_NUM);

	DrlFPrintf(fp, "ATM_RATE:\n");
	DrlTSwaptionMatrix2DFpWrite(that->fAtmRate,
		fp, TSWAPTION_MATRIX_FMT_NUM);

	return (SUCCESS);
}

/*----------------------------------------------------------------------
 * 
 */

int
DrlTSmile3DataReadFromExport(
	TSmile3Data *that,		/* (B) smile data */
	char *curCode,			/* (I) currency code (USD,etc)  */
	char *dirname,			/* (I) directory (or "nil") */
        TDate baseDate)                 /* (I) reference date */
{
static	char	routine[] = "DrlTSmile3DataReadFromKapitalExport";
	int	status = FAILURE;
	int	freq = 1;		/* vol frequency */
	FILE	*fp;
	char	cur[32], fnam[1024];

#define	FP_OPEN(x)	{if ((fp = fopen((x), "r")) == NULL) {\
			GtoErrMsg("%s: can't open `%s'.\n", routine,\
			x); goto done;}}
#define	FP_CLOSE	{if (fp) fclose(fp); fp = NULL;}

#define FP_ADVANCE(token) {if (DrlFAdvanceToToken(fp, token) != SUCCESS) { \
                        GtoErrMsg("%s: can't locate token `%s'\n", token, \
			routine);  goto done;}}



	/* cleanup */
	IF_FAILED_DONE(DrlTSmile3DataCleanup(that));

	/* currency code to lower case */
	strcpy(cur, curCode);
	{ char *p; for (p=cur; *p !='\0'; p++) *p = tolower(*p);}

	/* If directory name is "nil", don't load */
	strcpy(fnam, dirname);
	{ char *p; for (p=fnam; *p !='\0'; p++) *p = tolower(*p);}
	if (!strcmp(fnam, "nil")) {
		return (SUCCESS);
	}


	/* 3-pts smile */
	that->fSmileType = 2;

	/* Read q coeffs */
	sprintf(fnam, "%s/base_smq.%s", dirname, cur);
	IF_FAILED_DONE(DrlStrSubsEnv(fnam));
	FP_OPEN(fnam);
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpReadBVFmt(
		&that->fQcoeff, fp, baseDate, freq));
	FP_CLOSE;

	/* Read v+ coeffs */
	sprintf(fnam, "%s/base_vpl.%s", dirname, cur);
	IF_FAILED_DONE(DrlStrSubsEnv(fnam));
	FP_OPEN(fnam);
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpReadBVFmt(
		&that->fVhicoeff, fp, baseDate, freq));
	FP_CLOSE;


	/* Read v+ coeffs */
	sprintf(fnam, "%s/dr_smile.%s", dirname, cur);
	FP_OPEN(fnam);
	FP_ADVANCE("H_COEFF:");
	IF_FAILED_DONE(DrlFScanVType(fp, DRL_DOUBLE_T, &that->fHcoeff));

	FP_ADVANCE("ALPHA_COEFF:");
	IF_FAILED_DONE(DrlFScanVType(fp, DRL_DOUBLE_T, &that->fAlpha));

	FP_ADVANCE("BACKBONE_Q_COEFF:");
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpRead(&that->fBackBoneQ, fp,
			TSWAPTION_MATRIX_FMT_NUM));
	FP_ADVANCE("ATM_RATE:");
	IF_FAILED_DONE(DrlTSwaptionMatrix2DFpRead(&that->fAtmRate, fp,
			TSWAPTION_MATRIX_FMT_NUM));
	FP_CLOSE;


	/* Read ATM rates */

	/* Logging */
#ifdef	__DEBUG__
	DrlFPrintf(NULL, \
		"%s: read %s data from `%s' (baseDate %s):\n", \
		routine, curCode, dirname, \
		(baseDate > 0L ? DrlTDatePrint(NULL, baseDate) : "N/A"));\
	DrlTSmile3DataFpWrite(that, NULL);
#endif	/*__DEBUG__*/

	/* Made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);

#undef	FP_OPEN
#undef	FP_CLOSE
#undef	FP_ADVANCE
}


/*f--------------------------------------------------------------
 * Computes the smile volatility adjustment for a given
 * point interpolated off the volatility grid.\\
 */

int
DrlTSmile3DataInterp(
	TCurve *zcCurve,		/* (I) zero curve */
	TSmile3Data *that,		/* (I) smile data */
	TDate baseDate,			/* (I) base date */
	double strikeRate,		/* (I) strike rate */
	int freq,			/* (I) rate freq (0,1,2,4,12) */
	long dayCountConv,		/* (I) rate day count conv */
	TDate volDate,			/* (I) vol exp date */
	double volMat,			/* (I) vol fwd maturity */
	int volFreq,			/* (I) vol frequency */
	double volRate,			/* (I) vol */
	double *volAdj)			/* (O) vol smile adjustment */
{
static	char	routine[] = "DrlTSmile3DataInterp";
	int	status = FAILURE;

	double	atmVol = volRate,
		atmRate,
		v_hi,
		q, smadj;

	/* */
	if (that == NULL) {
		GtoErrMsg("%s: no smile data available.\n", routine);
		goto done;
	}

	/* if no smile, we are done */
	if (that->fSmileType == 0) {
		*volAdj = 0e0;
		status = SUCCESS;
		goto done;
	}


	/*
	 * compute forward ATM rate
	 */
	if (forwardATMRate(
		zcCurve,
		volDate,
		volMat,
		freq,
		dayCountConv,
		&atmRate) != SUCCESS)
			goto done;


	if (DrlTSwaptionMatrix2DInterpDate(
		that->fVhicoeff, &v_hi,
		baseDate, volDate, volMat,
		FALSE) != SUCCESS)
			goto done;

	if (DrlTSwaptionMatrix2DInterpDate(
		that->fQcoeff, &q,
		baseDate, volDate, volMat,
		FALSE) != SUCCESS)
			goto done;


	/* compute smile */
	if (SmileFunc3Pts(
		atmRate, strikeRate,
		v_hi,
		that->fHcoeff,
		q,
		that->fAlpha,
		that->fSmileType,
		&smadj) != SUCCESS) goto done;

	*volAdj = (smadj - 1e0) * atmVol;


#ifdef	__DEBUG__
#endif	/*__DEBUG__*/

	DrlFPrintf(NULL,
		"SMILE: %5.2fx%5.2f ATM=%5.2f%% STK=%5.2f%% ATMVOL=%5.2f%%: "
		"Q=%6.4f V+=%6.4f ALPHA=%6.4f H=%6.4f TYPE=%d:\n"
		"SMILE: ADHJ=%6.4f VOLADJ=%5.2f%% VOLSML=%5.2f%%\n",
		(volDate - zcCurve->fBaseDate) / 365e0, volMat, 
		atmRate*1e2, strikeRate*1e2, atmVol*1e2, 
		q, v_hi, that->fAlpha, that->fHcoeff, that->fSmileType,
		smadj, *volAdj*1e2, atmVol*1e2 + *volAdj*1e2);

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*f--------------------------------------------------------------
 * Convenience routine to compute a forward ATM rate.
 */

static	int
forwardATMRate(
	TCurve *zcCurve,	/* (I) zero curve */
	TDate startDate,	/* (I) reset date */
	double rateMat,		/* (I) rate maturity */
	int rateFreq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	long rateDayCountConv,	/* (I) day count convention */
	double *atmRate)	/* (I) forward atm rate */
{
static	char	routine[] = "forwardATMRate";
	int	status = FAILURE;
	TDate		matDate;
	TDateInterval	payInterval;


	if (GtoTDateAdvanceYears(startDate, rateMat, &matDate) != SUCCESS)
		goto done;

	switch (rateFreq) {
	case 1:
	case 2:
	case 4:
	case 12:
	    if (GtoFreq2TDateInterval((long)rateFreq, &payInterval) != SUCCESS)
		goto done;

	    if (GtoZerosToCouponsPoint(zcCurve,
			GTO_LINEAR_INTERP,
			startDate,
			&payInterval,
			matDate,
			rateDayCountConv,
			GTO_STUB_BOND,
			FALSE,
			atmRate) != SUCCESS) goto done;
	    break;

	case 0:
	    if (GtoZerosToSimplePoint(zcCurve,
			GTO_LINEAR_INTERP,
			startDate,
			matDate,
			rateDayCountConv,
			atmRate) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad rate freq %d.\n", routine, rateFreq);
		goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}




/*f--------------------------------------------------------------
 * Computes the smile volatility adjustment for an array of
 * volatility point interpolated off the volatility grid.
 * if {\tt volAdj} is not NULL, then {\tt volRates} is
 * unchangedand {\tt volAdj} contains the volatility differrences
 * due to the smile adjustment.
 * If it is NULL, then {\tt volRates} is
 * modified on exit (contains the asjuetd vol).
 */

int
DrlTSmile3DataInterpVolCurve(
	TSmile3Data *that,		/* (I) smile data */
	TCurve *zcCurve,		/* (I) zero curve */

	int freq,			/* (I) rate freq (0,1,2,4,12) */
	long dayCountConv,		/* (I) day count convention */

	int numVolDates,		/* (I) # of volatility dates */
	TDate *volDates,		/* (I) array of vol dates */
	double *strikeRates,		/* (I) array of strike rate */
	double *volMat,			/* (I) array of vol fwd maturities */
	int *volFreq,			/* (I) array of vol frequency */
	double *volRates,		/* (I/B) array of vol */
	double *volAdj)			/* (O) array of vol adj (or NULL)*/
{
static	char	routine[] = "DrlTSmile3DataInverpVolCurve";
	int	status = FAILURE;

	int	n;
	double	volAdjustment;

	/* interp the cube */
	for (n=0; n<=numVolDates-1; n++) {

	    IF_FAILED_DONE(DrlTSmile3DataInterp(
		zcCurve,
		that,
		zcCurve->fBaseDate,
		strikeRates[n],
		freq,
		dayCountConv,
		volDates[n],
		volMat[n],
		volFreq[n],
		volRates[n],
		(volAdj != NULL ? &volAdj[n] : &volAdjustment)));

		if (volAdj == NULL)
			volRates[n] += volAdjustment;

	}


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


