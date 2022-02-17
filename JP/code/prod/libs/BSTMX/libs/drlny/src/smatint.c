/****************************************************************
 * Module:	DRL
 * Submodule:	SMAT
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>	
#include <float.h>	
#include <stdarg.h>	

#include "gtomat.h"
#include "date_sup.h"
#include "convert.h"
#include "tcurve.h"
#include "check.h"		/* checkFrequency */

#include "drlinter.h"		/* DrlLinInterpXXX */
#include "drlmem.h"		/* DrlDoubleVectAlloc */
#include "drltime.h"		/* DrlTDateIntervalNil() */
#include "drlts.h"

#include "drlsmat.h"		/* Prototype consistency */

#undef	SWNEXP
#undef	SWNMAT
#undef	SWFREQ
#undef	SWTEXP
#undef	SWTMAT
#undef	SWVOL

#define	SWNEXP		(that->table->matrix->numDim1)
#define	SWNMAT		(that->table->matrix->numDim2)
#define	SWFREQ		(that->swapPayFreq)
#define	SWDIAG		(that->diagonal)
#define	SWTEXP(idxExp)	(that->table->dim1Values[idxExp])
#define	SWTMAT(idxMat)	(that->table->dim2Values[idxMat])
#define	SWVOL(idxExp, idxMat) \
			(that->table->matrix->data[idxExp][idxMat])

/*
static	int	_getInterpCoeffs(TSwaptionMatrix2D *that,
			double tExp, double tMat,
			int *ie, int *im, double *w);
static	int	_getMatInterpCoeffs(TSwaptionMatrix2D *that,
			int ie, double tMat,
			int *imlo, int *imhi,
			double *wmlo, double *wmhi);
*/

/*#define	__DEBUG__*/

/*---------------------------------------------------------------
 * Computes the maturity interpolation time on
 * a matrix according to the following conventions.
 * <br>
 * <br> <i> vertical matrix and final maturity:</i>
 * maturity interpolation time $T_m$ is obtained
 * by subtracting 
 * expiration time (ACT/365F from base date) from
 * final maturity time (ACT/365F from base date)
 * $$T_m = T_e - \bigl(\mbox{mat date}\bigr)_{ACT/365F}.$$
 * <br> <i> diagonal matrix and final maturity:</i>
 * maturity interpolation time $T_m$ is obtained from final
 * maturity date as ACT/365F.
 * $$T_m = \bigl(\mbox{mat date}\bigr)_{ACT/365F}.$$
 * <br> <i> vertical matrix and constant maturity:</i>
 * maturity interpolation time $T_m$ is obtained from 
 * maturity interval as 30/360.
 * $$T_m = \bigl(\mbox{mat interval}\bigr)_{30/360}.$$
 * <br> <i> diagonal matrix and constant maturity:</i>
 * maturity interpolation time $T_m$ is obtained
 * by first adding maturity interval
 * to expiration date (ACT/365F form base date)
 * to find final maturity date
 * and finding time to final maturity (ACT/365F from base date)
 * $$T_m = T_e + \bigl(\mbox{exp date} + \mbox{mat interval}
 * \bigr)_{ACT/365F}.$$
 * <br>
 * Also computes the forward rate maturity in year fraction (30/360)
 * and the rate frequency from the swaption matrix
 * (unless it is a CMS rate in vertical matrix with maturity
 * less than the frequency, in which case the rate is assumed to
 * be simple).
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DInterpMatTime(
	TSwaptionMatrix2D *that,/* (I) swaption matrix */
	TDate baseDate,		/* (I) base date */
	double tExp,		/* (I) time to exp (ACT/365F) */
	int finalFlag,		/* (I) TRUE=final maturity, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *interpMat,	/* (O) matrix mat interp time (can be NULL) */
	double *tMat,		/* (O) fwd mat (30/360) (can be NULL) */
	int *freq)		/* (O) rate freq (0,1,2,4,12) (can be NULL) */
{
static	char	routine[] = "DrlTSwaptionMatrix2DInterpMatTime";
	int	status = FAILURE;
	TDate	expDate;
	double	fm;

	/* compute interp time */
	if (interpMat != NULL) {
	if ((finalFlag == TRUE) && (that->diagonal == FALSE)) {
	    /*
	     * Final maturity in vertical matrix
	     */
	    IF_FAILED_DONE( GtoDayCountFraction(
		baseDate,
		matDate,
		GTO_ACT_365F,
		interpMat));

	    *interpMat -= tExp;

	} else if ((finalFlag == TRUE) && (that->diagonal == TRUE)) {
	    /* final maturity in diagonal matrix */
	    if (GtoDayCountFraction(baseDate, matDate,
		GTO_ACT_365F, interpMat)
		!= SUCCESS) goto done;

	} else if ((finalFlag == FALSE) && (that->diagonal == FALSE)) {
	    /* constant maturity in vertical matrix */
	    if (GtoDateIntervalToYears(&matInt, interpMat) != SUCCESS)
		goto done;

	} else if ((finalFlag == FALSE) && (that->diagonal == TRUE)) {
	    /* constant maturity in diagonal matrix */
	    expDate = baseDate + (long) (tExp*365e0);
	    if (GtoDtFwdAny(expDate, &matInt, &matDate)
		!= SUCCESS) goto done;

	    if (GtoDayCountFraction(baseDate, matDate,
		GTO_ACT_365F, interpMat)
		!= SUCCESS) goto done;
	}
	}

        /* compute forward maturity time of swaption */
	if (tMat != NULL) {
            if (finalFlag == FALSE) {
                if (GtoDateIntervalToYears(&matInt, tMat) != SUCCESS)
                        goto done;
            } else {
		expDate = baseDate + (long) (tExp*365e0);
                if (GtoDayCountFraction(expDate, matDate,
                        GTO_B30_360, tMat) != SUCCESS)
                        goto done;
            }
	}

        /* frequency: for a CMS rate in a vertical matrix,
	 * rates with maturity shorter than 1/freq are simple
	 */
	if (freq != NULL) {
            if ((finalFlag == FALSE) && (that->diagonal == FALSE)) {
                if (GtoDateIntervalToYears(&matInt, &fm) != SUCCESS)
                        goto done;
		*freq = (fm * ((double)that->swapPayFreq) < 1e0 ?
			0 :that->swapPayFreq);
            } else {
		*freq = that->swapPayFreq;
            }
	}

#if defined(__DEBUG__)
	fprintf(stdout, "%s: tExp=%5.2f x %10s : interpMat=%lf\n",
		routine, tExp,
		(finalFlag == FALSE ? GtoFormatDateInterval(&matInt) :
			GtoFormatDate(matDate)),
		*interpMat);
#endif

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);

}

/*f--------------------------------------------------------------
 * Performs a direct/adjoint interpolation of the whole
 * matrix <i> interpMat</i>. The date <i> baseDate</i> is used for date
 * interpolation consistency.
 * <br>
 * <br> If <i> adjointFlag</i> is FALSE, linearly interpolates the values
 * of the matrix "that" at the points of the matrix "interpMat".
 * <br> If <i> adjointFlag</i> is TRUE, rebuckets the values
 * of the matrix "that" at the points of the matrix "interpMat"
 * based on a linear interpolation from "interpMat" ("that" is unchanged,
 * and the interp values are added to "interpMat").
 * <br>
 * Returns 0 iff successful.
 */


DLL_EXPORT(int)
DrlTSwaptionMatrix2DInterpMatrix(
	TSwaptionMatrix2D *that,	/* (B) swaption matrix */
	TSwaptionMatrix2D *interpMat,	/* (B) interp swaption matrix */
        TDate baseDate,			/* (I) base date */
	int adjointFlag)		/* (I) TRUE=rebucket, FALSE=interp */
{
static	char	routine[] = "DrlTSwaptionMatrix2DInterpMatrix";
	int	status = FAILURE;

	TDate		matDate = 0L;
	TDateInterval	matInt = DrlTDateIntervalNil();
	int		i, j;
	double		tExp, tMat;

#define	IMNEXP		(interpMat->table->matrix->numDim1)
#define	IMNMAT		(interpMat->table->matrix->numDim2)
#define	IMFREQ		(interpMat->swapPayFreq)
#define	IMDIAG		(interpMat->diagonal)
#define	IMTEXP(idxExp)	(interpMat->table->dim1Values[idxExp])
#define	IMTMAT(idxMat)	(interpMat->table->dim2Values[idxMat])
#define	IMVOL(idxExp, idxMat) \
			(interpMat->table->matrix->data[idxExp][idxMat])
#ifdef	__DEBUG__
	fprintf(stdout, "%s: Base date %s.\n", routine,
		GtoFormatDate(baseDate));
#endif


	if (adjointFlag != FALSE) {
	    /* adjoint interpolation */
	    /* reset to zero */
	    if (DrlTSwaptionMatrix2DOperScalar(interpMat, "=", 0e0) != SUCCESS)
		goto done;

#ifdef	__DEBUG__
	    fprintf(stdout, "%s: Adjoint Interpolation:\n", routine);
	    fprintf(stdout, "%s: input matrix (that)\n", routine);
            DrlTSwaptionMatrix2DFpWrite(that, stdout, TSWAPTION_MATRIX_FMT_STD);
	    fprintf(stdout, "%s: output matrix (interpMat)\n", routine);
            DrlTSwaptionMatrix2DFpWrite(interpMat, stdout,
			TSWAPTION_MATRIX_FMT_STD);
#endif

	    for (i=0; i<=SWNEXP-1; i++)
	    for (j=0; j<=SWNMAT-1; j++) {

		tExp = SWTEXP(i);

		/* compute maturity date/interval */
		if (SWDIAG == TRUE) {
			matDate = baseDate + (long) (365e0 * SWTMAT(j));
		} else {
			if (GtoYearsToDateInterval(SWTMAT(j), &matInt)
			    != SUCCESS) goto done;
		}

		/* compute interpolation maturity */
		if (DrlTSwaptionMatrix2DInterpMatTime(
			interpMat,
			baseDate,
			tExp,
			SWDIAG, matDate, matInt,
			&tMat, NULL, NULL) != SUCCESS)
				goto done;

		/* perform interp */
		if (DrlTSwaptionMatrix2DInterpExpMat(
			interpMat, 
			&SWVOL(i,j),
			tExp,
			tMat,
			TRUE) != SUCCESS)
				goto done;

	    }
	} else {
	    /* direct interpolation */
	    /* set to 0 */
	    if (DrlTSwaptionMatrix2DOperScalar(interpMat, "=", 0e0) != SUCCESS)
		goto done;
#ifdef	__DEBUG__
	    fprintf(stdout, "%s: Direct Interpolation:\n", routine);
	    fprintf(stdout, "%s: input matrix (that)\n", routine);
            DrlTSwaptionMatrix2DFpWrite(that, stdout, TSWAPTION_MATRIX_FMT_STD);
	    fprintf(stdout, "%s: output matrix (interpMat)\n", routine);
            DrlTSwaptionMatrix2DFpWrite(interpMat, stdout,
			TSWAPTION_MATRIX_FMT_STD);
#endif

	    for (i=0; i<=IMNEXP-1; i++)
	    for (j=0; j<=IMNMAT-1; j++) {

		tExp = IMTEXP(i);

		/* compute maturity date/interval */
		if (IMDIAG == TRUE) {
			matDate = baseDate + (long) (365e0 * IMTMAT(j));
		} else {
			if (GtoYearsToDateInterval(IMTMAT(j), &matInt)
			    != SUCCESS) goto done;
		}

		/* compute interpolation maturity */
		if (DrlTSwaptionMatrix2DInterpMatTime(
			that,
			baseDate,
			tExp,
			IMDIAG, matDate, matInt,
			&tMat, NULL, NULL) != SUCCESS)
				goto done;

		/* perform interp */
		if (DrlTSwaptionMatrix2DInterpExpMat(
			that, 
			&IMVOL(i,j),
			tExp,
			tMat,
			FALSE) != SUCCESS)
				goto done;

	    }
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
#undef	IMNEXP
#undef	IMNMAT
#undef	IMFREQ
#undef	IMDIAG
#undef	IMTEXP
#undef	IMTMAT
#undef	IMVOL
}


/*---------------------------------------------------------------
 * Computes interpolation indices and weights.
 */


DLL_EXPORT(int)
DrlGetInterpCoeffs(
	TSwaptionMatrix2D *that,/* (I) swaption matrix */
	double tExp,		/* (I) time to exp */
	double tMat,		/* (I) maturity interpolation time */
	int *ie,		/* (O) exp idx [4] */
	int *im,		/* (O) mat idx [4] */
	double *w)		/* (O) weights [4] */
{
static	char	routine[] = "DrlGetInterpCoeffs";
	int	status = FAILURE;

	int	ilo, ihi;
	double	wlo, whi;

	/* interp expiration */
	if (DrlLinearInterp1dWeights(&SWTEXP(0), SWNEXP,
		tExp,
		&ilo, &ihi,
		&wlo, &whi) != SUCCESS)
				goto done;
	ie[0] = ilo;
	ie[1] = ilo;
	ie[2] = ihi;
	ie[3] = ihi;

	w[0] = wlo;
	w[1] = wlo;
	w[2] = whi;
	w[3] = whi;


	/* interp mat on lo exp point */
	if (DrlGetMatInterpCoeffs(that, ie[0], tMat,
		&ilo, &ihi,
		&wlo, &whi) != SUCCESS) goto done;

	im[0] = ilo;
	im[1] = ihi;
	w[0] *= wlo;
	w[1] *= whi;

	/* interp mat on hi exp point */
	if (DrlGetMatInterpCoeffs(that, ie[2], tMat,
		&ilo, &ihi,
		&wlo, &whi) != SUCCESS) goto done;

	im[2] = ilo;
	im[3] = ihi;
	w[2] *= wlo;
	w[3] *= whi;

	if (((that->diagonal != FALSE) && (tExp > tMat)) ||
	    ((that->diagonal == FALSE) && (tMat < 0e0 ))) {
		w[0] = 0e0;
		w[1] = 0e0;
		w[2] = 0e0;
		w[3] = 0e0;
	}
	


#ifdef	__DEBUG__
	fprintf(stdout, "%s: input Matrix\n", routine);
        DrlTSwaptionMatrix2DFpWrite(that, stdout, TSWAPTION_MATRIX_FMT_STD);
	fprintf(stdout, "\ttExp=%8.4f x tMat=%8.4f: \n", tExp, tMat);
	for (ilo=0; ilo<=3; ilo++) {
	   fprintf(stdout, "\tie=%2d  im=%2d  w=%lf  \n",
		ie[ilo], im[ilo], w[ilo]);
	}
#endif

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);

}

/*---------------------------------------------------------------
 * Computes interpolation indices and weights
 * in the maturity dimension on the matrix <i> that</i>
 * at expiration time of index <i> ie</i>.
 */

DLL_EXPORT(int)
DrlGetMatInterpCoeffs(
	TSwaptionMatrix2D *that,/* (I) swaption matrix */
	int ie,			/* (I) exp idx */
	double tMat,		/* (I) mat to interp */
	int *imlo,		/* (O) lo mat idx */
	int *imhi,		/* (O) hi mat idx */
	double *wmlo,		/* (O) lo mat weight */
	double *wmhi)		/* (O) hi mat weight */
{
static	char	routine[] = "DrlGetMatInterpCoeffs";
	int	status = FAILURE;

	int	io;		/* offset index */

	if (that->diagonal) {
	    /* diagonal matrix: * offset maturities
	     * find the first maturity index larger or equal
	     * to the expiration.
	     */
	    for (io=0; io<=SWNMAT-1; io++) {
		if (SWTMAT(io) >= SWTEXP(ie))
		    break;
	    }

	    /* interp the maturities on the offset vector
	     * so that no weight is put on exp/final mat points
	     * that have negative forward maturities.
	     */
	    if (DrlLinearInterp1dWeights(&SWTMAT(io), SWNMAT-io,
		tMat,
		imlo, imhi,
		wmlo, wmhi) != SUCCESS)
				goto done;
	    *imlo += io;
	    *imhi += io;

	} else {
	    if (DrlLinearInterp1dWeights(&SWTMAT(0), SWNMAT,
		tMat,
		imlo, imhi,
		wmlo, wmhi) != SUCCESS)
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
 * Performs a direct/rebucketing interpolation of the point
 * of coordinates <i> tExp</i> and <i> tMat</i>.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DInterpExpMat(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	double *value,		/* (B) interpolated value */
	double tExp,		/* (I) time to expiration */
	double tMat,		/* (I) maturity interpolation time */
	int adjointFlag)	/* (I) TRUE=adjoint, FALSE=direct*/
{
static	char	routine[] = "DrlTSwaptionMatrix2DInterpExpMat";
	int	status = FAILURE;
	int	k, ie[4], im[4];
	double	w[4];

	/* get interpolation coefficients */
	if (DrlGetInterpCoeffs(that, tExp, tMat, ie, im, w) != SUCCESS)
		goto done;

	if (adjointFlag == FALSE) {
	    /* direct interpolation */
	    *value = 0e0;
	    for (k=0; k<=3; k++) {
		*value += SWVOL(ie[k], im[k]) * w[k];
	    }
	} else {
	    /* adjoint interpolation */
	    for (k=0; k<=3; k++) {
		SWVOL(ie[k], im[k]) += (*value) * w[k];
	    }
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
 * Same as <i> DrlTSwaptionMatrix2DInterpExpMat</i>,
 * but takes expiration date.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DInterpDate(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	double *value,		/* (B) interpolated value */
	TDate baseDate,		/* (I) base date */
	TDate expDate,		/* (I) expiration date */
	double tMat,		/* (I) maturity interpolation time */
	int adjointFlag)	/* (I) TRUE=adjoint, FALSE=direct*/
{
	double	tExp;

	if (GtoDayCountFraction(baseDate, expDate, GTO_B30_360,
			&tExp) != SUCCESS)
				return(FAILURE);

	if (DrlTSwaptionMatrix2DInterpExpMat(
		that,
		value,
		tExp,
		tMat,
		adjointFlag) != SUCCESS)
				return(FAILURE);

	return(SUCCESS);
}


/*f--------------------------------------------------------------
 * Interpolates a volatility point on a swaption matrix 
 * for a constant or final maturity.
 * <br>
 * <br> If "adjointFlag" is FALSE, linearly interpolates the values
 * of the matrix "that" at the point
 * ("that" is NOT changed, "\*value" IS changed).
 * <br> If "adjointFlag" is TRUE, rebuckets (adds) the values "value"
 * at the corresponding interpolation points
 * ("that" IS changed, "\*value" is NOT changed).
 * <br>
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DInterpValue(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	TDate baseDate,		/* (I) base date */
	double tExp,		/* (I) time to exp (ACT/365F) */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *value,		/* (B) interpolated volatility */
	double *fwdMat,		/* (I) fwd maturity for interp (or NULL) */
	int adjointFlag)	/* (I) TRUE=rebucket, FALSE=interp */
{
static	char	routine[] = "DrlTSwaptionMatrix2DInterpVol";
	int	status = FAILURE;

	double		tMat;

	/* get interp time */
	if (DrlTSwaptionMatrix2DInterpMatTime(
		that,
		baseDate,
		tExp,
		finalFlag,
		matDate,
		matInt,
		&tMat, fwdMat, NULL) != SUCCESS)
			goto done;

	/* do the interp */
	if (DrlTSwaptionMatrix2DInterpExpMat(
		that,
		value,
		tExp,
		tMat,
		adjointFlag)
		!= SUCCESS) goto done;

#if defined(__DEBUG__)
	fprintf(stdout, "%s: tExp=%5.2f x %10s : value=%lf\n",
		routine, tExp,
		(finalFlag == FALSE ? GtoFormatDateInterval(&matInt) :
			GtoFormatDate(matDate)),
		*value);
#endif

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);

}

/*f--------------------------------------------------------------
 * Interpolates a volatility point on a swaption matrix 
 * for a constant or final maturity.
 * <br>
 * <br> If "adjointFlag" is FALSE, linearly interpolates the values
 * of the matrix "that" at the point
 * ("that" is NOT changed, "\*value" IS changed).
 * <br> If "adjointFlag" is TRUE, rebuckets (adds) the values "value"
 * at the corresponding interpolation points
 * ("that" IS changed, "\*value" is NOT changed).
 * <br>
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DDateInterpValue(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	TDate baseDate,		/* (I) base date */
	TDate expDate,		/* (I) expiration date */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *value,		/* (B) interpolated volatility */
	int adjointFlag)	/* (I) TRUE=rebucket, FALSE=interp */
{
static	char	routine[] = "DrlTSwaptionMatrix2DDateInterpValue";
	int	status = FAILURE;

	double		tExp, tMat;

	/* get interp time */
	IF_FAILED_DONE(GtoDayCountFraction(
		baseDate,
		expDate,
		GTO_B30_360,
		&tExp));

	if (finalFlag == TRUE) {
		/* final maturity */
		if (that->diagonal == FALSE) {
			/* final maturity in vertical matrix */
	    		IF_FAILED_DONE(GtoDayCountFraction(
				expDate,
				matDate,
				GTO_B30_360,
				&tMat));
		} else {
			/* final maturity in diagonal matrix */
			IF_FAILED_DONE(GtoDayCountFraction(
				baseDate,
				matDate,
				GTO_ACT_365F,
				&tMat));
		}
	} else {
		/* constant maturity */
		/* get time to maturity (30/360) */
		IF_FAILED_DONE(GtoDateIntervalToYears(
			&matInt,
			&tMat));

		if (that->diagonal != FALSE) {
			/* constant maturity in diagonal matrix */
			tMat -= tExp;
		}
	}

	/* do the interp */
	if (DrlTSwaptionMatrix2DInterpExpMat(
		that,
		value,
		tExp,
		tMat,
		adjointFlag)
		!= SUCCESS) goto done;


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);

}


/*f--------------------------------------------------------------
 * Interpolates a volatility curve on a swaption matrix 
 * for a constant or final maturity.
 * <br>
 * <br> If "adjointFlag" is FALSE, linearly interpolates the values
 * of the matrix "that" at the point of the vol curve
 * ("that" is unchanged).
 * <br> If "adjointFlag" is TRUE, rebuckets the values
 * of the array "volRates" (or the curve "volCurve" if "valRates" is
 * NULL) on  the matrix "that" ("that" IS changed). 
 * <b> WARNING:</b> in this case, no check is made on arrays length.
 * <br>
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DInterpVolCurve(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	TDate baseDate,		/* (I) base date */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */

	double tMatMin,		/* (I) Minimum muturity 		*/
	int adjointFlag,	/* (I) TRUE=rebucket, FALSE=interp */

	int *nVols,		/* (O) # of vol points (can be NULL) */
	TDate **volDates,	/* (O) array of vol exp dates (can be NULL) */
	double **volExp,	/* (O) array of vol exp time (can be NULL) */
	double **volMat,	/* (O) array of vol fwd mat (can be NULL) */
	int **volFreq,		/* (O) array of vol freq (can be NULL) */
	double **volRates,	/* (B) array of vol (can be NULL) */

	TCurve **volCurve)	/* (B) interpolated TCurve (can be NULL) */
{
static	char	routine[] = "DrlTSwaptionMatrix2DInterpVolCurve";
	int	status = FAILURE;

	int		n, nActive = 0;
	double		tMat, vol;
	double		tMatMax = SWTMAT(SWNMAT-1),
			*vMat = NULL;

	/* allocate memory */
	if ((vMat = DrlDoubleVectAlloc(0, SWNEXP-1)) == NULL)
		goto done;


	/* Find num of  points in the interpolated curve
	 */
	if (finalFlag == TRUE) {
#ifndef	NO_MIN_MAT
		/* find time to maturity (30/360) */
		if (GtoDayCountFraction(baseDate, matDate, GTO_B30_360,
			&tMat) != SUCCESS) goto done;

		/* check maturity OK */
		if (tMat <= 0e0) {
		     GtoErrMsg("%s: too short maturity date %s (%lf).\n",
			routine, GtoFormatDate(matDate), tMat);
		     goto done;
	        }

		/* find # of active points : Use all maturities */
		nActive = SWNEXP;

#else	/*NO_MIN_MAT*/
		/* Final Maturity */
		/* find time to maturity (30/360) */
		if (GtoDayCountFraction(baseDate, matDate, GTO_B30_360,
			&tMat) != SUCCESS) goto done;

		/* check maturity OK */
		if (tMat <= 0e0) {
		     GtoErrMsg("%s: too short maturity date %s (%lf).\n",
			routine, GtoFormatDate(matDate), tMat);
		     goto done;
	        }

		/* find # of active points */
		nActive = 0;
		for (n=0; n<=SWNEXP-1;n++) {
			if (tMat - SWTEXP(n) <= 0e0)
				break;
			nActive++;
		}

		if (nActive == 0) {
		     GtoErrMsg("%s: too short maturity date %s (%lf yrs),"
			"first maturity on matrix is %lf yrs.\n",
			routine, GtoFormatDate(matDate), tMat, SWTEXP(0));
		     goto done;
	        }
#endif/*NO_MIN_MAT*/

	} else {
		/* Constant Maturity */
		nActive = SWNEXP;
		/* get time to maturity (30/360) */
		if (GtoDateIntervalToYears(&matInt, &tMat) != SUCCESS)
			goto done;
	}



	/* Compute interpolation time and cutoff */
	for (n=0; n<=nActive-1; n++) {
	    if (DrlTSwaptionMatrix2DInterpMatTime(
		that,
		baseDate,
		SWTEXP(n),
		finalFlag,
		matDate,
		matInt,
		&vMat[n],
		NULL,
		NULL) != SUCCESS)
			goto done;

	    /* Cutoff at min and max maturities.
	     * Only adjust for final maturity
	     */
	    if (finalFlag == TRUE){ 
	    	vMat[n] = MIN(vMat[n], tMatMax);
	    	vMat[n] = MAX(vMat[n], tMatMin);
	    }

	}

	/* perform interpolation/adjoint */

	if (adjointFlag == FALSE) {
	    /**
	     ** regular interpolatin
	     **/
	    if (volDates != NULL) *volDates = NULL;
	    if (volExp   != NULL) *volExp = NULL;
	    if (volMat   != NULL) *volMat = NULL;
	    if (volFreq  != NULL) *volFreq = NULL;
	    if (volRates != NULL) *volRates = NULL;
	    if (volCurve != NULL) *volCurve = NULL;


	    /* do proper allocations */
	    if (volDates != NULL) {
		if ((*volDates = NEW_ARRAY(TDate, nActive)) == NULL)
			goto done;
		for (n=0; n<=nActive-1; n++) {
		    if (GtoTDateAdvanceYears(
			baseDate,
			SWTEXP(n),
			&(*volDates)[n]) != SUCCESS)
				goto done;
		    }
	    }
	    if (volExp != NULL) {
		if ((*volExp = NEW_ARRAY(double , nActive)) == NULL)
			goto done;
		for (n=0; n<=nActive-1; n++) {
			(*volExp)[n] = SWTEXP(n);
		}
	    }
	    if (volFreq != NULL) {
		if ((*volFreq= NEW_ARRAY(int, nActive)) == NULL)
			goto done;
		for (n=0; n<=nActive-1; n++) {
			(*volFreq)[n] = SWFREQ;
		}
	    }
	    if (volMat != NULL) {
		if ((*volMat = NEW_ARRAY(double , nActive)) == NULL)
			goto done;
		for (n=0; n<=nActive-1; n++) {
			(*volMat)[n] = vMat[n];
		}
	    }
	    if (volRates != NULL) {
		if ((*volRates = NEW_ARRAY(double , nActive)) == NULL)
			goto done;
	    }


	    /* create the TCurve */
	    if (volCurve != NULL) {
		if ((*volCurve = GtoNewTCurve(
			baseDate,
			nActive,
			(long) SWFREQ,
			GTO_ACT_365F)) == NULL)
				goto done;
		for (n=0; n<=nActive-1; n++) {
			if (GtoTDateAdvanceYears(
				baseDate,
				SWTEXP(n),
				&(*volCurve)->fArray[n].fDate) != SUCCESS)
					goto done;
		}
	    }

	    if (nVols != NULL) *nVols = nActive;

	    /* Perform the interpolation */
	    for (n=0; n<=nActive-1; n++) {

		/* do the interp */
		if (DrlTSwaptionMatrix2DInterpExpMat(
			that,
			&vol,
			SWTEXP(n),
			vMat[n],
			adjointFlag)
			!= SUCCESS)
				goto done;

		if (volCurve != NULL)
			(*volCurve)->fArray[n].fRate = vol;
		if (volRates != NULL)
			(*volRates)[n] = vol;
	    }

#if defined(__DEBUG__)
	    fprintf(stdout, "%s:\n", routine);
	    if (volCurve != NULL) DrlTCurveFpWrite(*volCurve, stdout,
			DRL_TCURVE_FMT_STD);
#endif


	} else {
	    /**
	     ** adjoint interpolation
	     **/
	    for (n=0; n<=nActive-1; n++) {
		/* */
		if (volRates != NULL) {
			vol = (*volRates)[n];
		} else if (volCurve != NULL) {
			vol = (*volCurve)->fArray[n].fRate;
		} else {
			GtoErrMsg("%s: adjoint interpolation, but volRates and"
				" volCurve both NULL.\n", routine);
				goto done;
		}


		/* do the interp */
		if (DrlTSwaptionMatrix2DInterpExpMat(
			that,
			&vol,
			SWTEXP(n),
			vMat[n],
			TRUE)
			!= SUCCESS)
				goto done;
	    }
	}


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(vMat, 0, SWNEXP-1);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*f--------------------------------------------------------------
 * Performs the point by point scalar product of
 * two swaption matrix "that" and "mat2" (that MUST have
 * same number of points, same expirations and maturities):
 * if $(a_{ij})$ and $(b_{ij})$ are the matrices values,
 * the product is defined by
 * <br>
 * \sum_{i=1}^{\#exp} \sum_{j=1}^{\#mat} a_{ij}b_{ij}
 * <br>
 * The output result is put in "val".
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DL2Product(
	TSwaptionMatrix2D *that,	/* (I) matrix #1 */
	TSwaptionMatrix2D *mat2,	/* (I) matrix #2 */
	double *val)			/* (O) */
{
static	char	routine[] = "DrlTSwaptionMatrix2DL2Product";
register int	i, j;

#define	SW2NEXP		(mat2->table->matrix->numDim1)
#define	SW2NMAT		(mat2->table->matrix->numDim2)
#define	SW2FREQ		(mat2->swapPayFreq)
#define	SW2TEXP(idxExp)	(mat2->table->dim1Values[idxExp])
#define	SW2TMAT(idxMat)	(mat2->table->dim2Values[idxMat])
#define	SW2VOL(idxExp, idxMat) \
			(mat2->table->matrix->data[idxExp][idxMat])

	if (!DrlTSwaptionMatrix2DIsSameType(that, mat2)) {
		GtoErrMsg("%s: matrices of different type.\n", routine);
		return(FAILURE);
	}

	*val = 0e0;
	for (i=0; i<=SWNEXP-1; i++)
	for (j=0; j<=SWNMAT-1; j++) {
		*val += SWVOL(i,j)*SW2VOL(i,j);
	}

	return(SUCCESS);
#undef	SW2NEXP
#undef	SW2NMAT
#undef	SW2FREQ
#undef	SW2TEXP
#undef	SW2TMAT
#undef	SW2VOL
}


