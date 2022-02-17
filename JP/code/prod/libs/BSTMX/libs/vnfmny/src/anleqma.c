/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>

#include "convert.h"

#include "drltime.h"
#include "drlio.h"
#include "drlsort.h"
#include "drlinter.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"

#define	__DEBUG__
#undef	__DEBUG__


#undef	K
#define K(a,b,alpha)    (ABS((alpha)*((b)-(a))) < 1e-8 ? \
                        (exp(-(alpha)*(a))*(  \
                         (b) - (a) - 0.5*(alpha)*((b)-(a))*((b)-(a))) ) : \
                        (exp(-(alpha)*(a)) - exp(-(alpha)*(b))) / (alpha))


#undef	K1
#define K1(a,b,alpha)    (ABS((alpha)*((b)-(a))) < 1e-8 ? \
			(1e0 + (alpha)*(a))*exp(-(alpha)*(a)) * ( \
				((b)-(a))*(a) + \
				((b)-(a))*((b)-(a))*0.5) \
			: \
			(( (1e0 + (alpha)*(a))*exp(-(alpha)*(a)) \
			  -(1e0 + (alpha)*(b))*exp(-(alpha)*(b))) \
			/ ((alpha)*(alpha))))


/*f-------------------------------------------------------------
 * Calculate the A1 coefficient.
 *                                                             
 * <br><br>
 * Compute and returns coefficient $A$ (used for computation
 * of the MM rates volatility) defined by
 * <blockquote>
 * A_j(S,T) = integral_S^T r^{1-2q}(t) e^{-beta_j*(t-S)} dt.
 * </blockquote>
 * where q is the backboneq with cevPower = 1-2q
 * (0=lognormal, 0.5=normal).
 *
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).
 * Returns 0 iff OK.
 */

int
VnfmA1(
	VnfmData *that,		/* (I) model parameters */
	double S,		/* (I) forward rate start */
	double T,		/* (I) forward rate end */
	double *A,		/* (O) A_i(S,T) */
	double *A1)		/* (O) dA_i/dbeta(S,T) */
{
	int	i, j, m, n,
		nDim = NDIM;

	DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);
	DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);

	for (j=0; j<=nDim-1; j++) {

	    A[j]  = 0e0;
	    A1[j] = 0e0;

	    for (i=n+1;i<m;i++) {
		A[j]  += FN(RATE[i],that->fBackBoneQ)
				* K(ZTT[i]-S, ZTT[i+1]-S, BETA[j]);
		A1[j] += FN(RATE[i],that->fBackBoneQ)
				* K1(ZTT[i]-S, ZTT[i+1]-S, BETA[j]);
	    }

	    if (n+1<=m) {
		/* L(..) same as K(0e0, ZTT[n+1]-S, BETA[j]) */
		A[j] += FN(RATE[n],that->fBackBoneQ)
				* L((ZTT[n+1]-S), BETA[j]);
		A[j] += FN(RATE[m],that->fBackBoneQ)
				* K(ZTT[m]-S, T-S, BETA[j]);

		A1[j] += FN(RATE[n],that->fBackBoneQ)
				* K1(0e0,      ZTT[n+1]-S, BETA[j]);
		A1[j] += FN(RATE[m],that->fBackBoneQ)
				* K1(ZTT[m]-S, T-S,        BETA[j]);


	    } else {
		/* L(..) same as K(0e0, T-S, BETA[j]) */
		A[j] += FN(RATE[n],that->fBackBoneQ)
				* L((T-S), BETA[j]);
		A1[j] += FN(RATE[n],that->fBackBoneQ)
				* K1(0e0, T-S, BETA[j]);

	    }
 	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "\tVnfmA1: S=%8.4f T=%8.4f n=%2d m=%2d : ",
		S, T, m, n);
	for (j=0; j<=nDim-1; j++)
		DrlFPrintf(vnfmFpLog, "%12.8f  %12.8f", A[j], A1[j]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif

#undef	FN

	return(SUCCESS);
}




/*f-------------------------------------------------------------
 * Calculate the B1 coefficient.
 *                                                             
 * <br><br>
 * Computes and returns the V-coefficients for coupon yield volatility.
 * "S" is the time to reset of the option,
 * "tMat" the forward maturity of the bond,
 * "freq" its frequency (1,2,4,12).\\
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).
 * Returns 0 iff OK.
 * <b>
 * WARNING: the $B$ in the code is the B of the note
 * multiplied by the $1/(1-Z)$ factor. </b>
 */

int
VnfmBulletB1(
	VnfmData *that,	/* (I) model parameters */
	double S,	/* (I) expiration */
	double cpn,	/* (I) coupon or -1e0 for forward */
	double tMat,	/* (I) fwd bond maturity */
	int freq,	/* (I) bond frequency */
	double *B,	/* (O) B[i](S,tMat) */
	double *B1)	/* (O) dB/dbeta[i](S,tMat) */
{
	/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*
	 *	New Version:					*
	 *	We do simple stubs to make it consistent	*
	 *	with simple rate				*
	 *xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

	int	nDim = NDIM,
		m, n, idxC, nC, idxF;
	double	T,
		A[VNFM_NDIMMAX],  An[VNFM_NDIMMAX],	/* A coefficients */
		A1[VNFM_NDIMMAX], A1n[VNFM_NDIMMAX],	/* A1 coefficients */
		z, zn,					/* forward zero */
		b[VNFM_NDIMMAX],  b1[VNFM_NDIMMAX],	/* temporar storage */
		dur,					/* duration */
		dcf,					/* accrual */
		y,					/* pay yield */
		dt;
#define	ONE_DAY	2.73972602e-3

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmBulletB1: ");
#endif

	dt      = 1e0/freq;
	nC      = (int) floor((tMat-ONE_DAY)/dt) + 1;	/* # coupons */
	dur     = 0e0;
	for (idxF=0; idxF<=nDim-1; idxF++) {
		b[idxF]  = 0e0;
		b1[idxF] = 0e0;
	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "S=%lf cpn=%8.4f tMat=%lf freq=%d : nC=%d\n",
		S, cpn, tMat, freq, nC);
#endif


	for (idxC=1; idxC<=nC; idxC++) {
		/* Coupon payment date */
		T = (tMat - (nC - idxC) * dt) + S;

		/* compute A coefficients */
		VnfmA1(that, S, T, A, A1);

		/* */
		DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);
		DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);

		/* compute forward zero */
		z = (that->fZero[m] / that->fZero[n])
			*exp( - RATE[m] * (T - ZTT[m])
			      + RATE[n] * (S - ZTT[n]) );

		/* last coupon: add principal */
		if (idxC == nC) {
			for (idxF=0; idxF<=nDim-1; idxF++) {
				An[idxF]  = A[idxF];
				A1n[idxF] = A1[idxF];
			}
			zn  = z;
		}

		/* 1s coupon: case of a front stub */
		if ((idxC == 1) && (T - dt < S)) {
			dcf = (T - S);
		} else {
			dcf = dt;
		}


#ifdef	__DEBUG__
		DrlFPrintf(vnfmFpLog, "VnfmB(%3d): cpn %2d/%2d "
			"S=%7.4f T=%7.4f (%7.4f)  n=%2d m=%2d dcf=%lf "
			"z=%lf zrcc=%lf zr1c=%lf\n",
			__LINE__, idxC, nC, S, T, T-S,
			n, m, dcf, z, -log(z) / (T-S),
			pow(z, -1./(T-S))-1.);
#endif

		/* compute duration and beta-duration */
		dur += z * dcf;
		for (idxF=0; idxF<=nDim-1; idxF++) {
			b[idxF]   += A[idxF]  * z * dcf;
			b1[idxF]  += A1[idxF] * z * dcf;
		}
	}

	/* Compute par yield: simple stub */
	if (cpn < 0e0) {
		y = (1e0 - zn) / dur;
	} else {
		y = cpn;
	}

#ifdef	__DEBUG__
	/*DrlFPrintf(vnfmFpLog, "VnfmB1(%3d):\n", __LINE__);
	for (idxF=0; idxF<=nDim-1; idxF++)
		DrlFPrintf(vnfmFpLog, "\tB[%2d] =%12.8f ", idxF, b[idxF]);
	DrlFPrintf(vnfmFpLog, "\n");
	for (idxF=0; idxF<=nDim-1; idxF++)
		DrlFPrintf(vnfmFpLog, "\tB1[%2d]=%12.8f ", idxF, b1[idxF]);
	DrlFPrintf(vnfmFpLog, "\n");*/
#endif

	for (idxF=0; idxF<=nDim-1; idxF++) {
		B[idxF]  = (y * b[idxF]  + An[idxF]  * zn);
		B1[idxF] = (y * b1[idxF] + A1n[idxF] * zn);
	}


#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog,
		"VnfmB1(%3d): S=%lf T=%lf freq=%3d y=%lf dur=%lf acc=% | ",
		__LINE__, S, tMat, freq, y, dur, (S-T));
	for (idxF=0; idxF<=nDim-1; idxF++)
		DrlFPrintf(vnfmFpLog, " [%2d] %12.8f %12.8f (%12.8f)",
			idxF, B[idxF], B1[idxF],
			B1[idxF]/B[idxF]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif


	return(SUCCESS);
}



/*f-------------------------------------------------------------
 *
 */

int
VnfmSolveBulletEqmat(
	VnfmData *that,	/* (I) model parameters */
	double S,	/* (I) expiration */

	double *Q,	/* (I) input B[i](S,tMat) */
	double *Q1,	/* (I) input dB/dbeta[i](S,tMat) */

	double tMatMin,	/* (I) minimum equivalent maturity */
	double tMatMax,	/* (I) minimum equivalent maturity */
	int freq,	/* (I) frequency for equivalent bullet */

	double *teqMat)	/* (I) fwd bond maturity */
{
static	char	routine[] = "VnfmSolveBulletEqmat";
	int	status = FAILURE;

	int	nDim = NDIM;
	double	B[VNFM_NDIMMAX],
		B1[VNFM_NDIMMAX];	/* bullet V-coefficients */

	int	nC, nCMin, nCMax;
	double	tMat, tMat0,
		zB, zB0,
		zA;
	int	nIter = 0;

	double	tMin, tMax,
			zBMin, zBMax;


	/* Checj 1F */
	if (nDim != 1) {
		GtoErrMsg("%s: num fact (%d) != 1.\n", routine, nDim);
		goto done;
	}
	if (tMatMin > tMatMax) {
		GtoErrMsg("%s: min mat (%lf) > max mat (%lf).\n",
			routine, tMatMin, tMatMax);
		goto done;
	}

	if (tMatMin > 1e0/(double)freq) {
		GtoErrMsg("%s: min mat (%lf) > 1/freq (%lf).\n",
			routine, tMatMin, (1e0/(double)freq));
		goto done;
	}


	zA = Q1[0] / Q[0];

	nCMin = (int) floor(tMatMin * freq + 1e-3);
	nCMin = MAX(nCMin, 1);
	nCMax = (int) floor(tMatMax * freq + 1e-3);

#undef	_OLDSCHEME
#define	_OLDSCHEME
#ifdef	_OLDSCHEME
	for (nC = nCMin; nC <= nCMax; nC++) {

		tMat = ((double) nC) / ((double) freq);

		/*GtoErrMsg("nIter=%3d nCMin=%3d mCMax=%3d nC=%3d\n",
			nIter, nCMin, nCMax, nC); */
		IF_FAILED_DONE( VnfmBulletB1(
			that,
			S,
			-1e0,	/* fwd ATM */
			tMat,
			freq,
			B,
			B1));

		zB = B1[0] / B[0];

#ifdef	__DEBUG__
		GtoErrMsg("%s: nC=%d tMat=%lf zA=%lf zB=%lf\n", routine,
			nC, tMat, zA, zB);
#endif


		/* found first larger than target */
		if (zB > zA) {
			if (nC != nCMin) {
				tMat = tMat0 + (zA - zB0) / (zB - zB0)
					/ ((double) freq);
			} else {
				tMat = tMatMin;
			}
			break;
		}
		zB0 = zB;
		tMat0 = tMat;

		nIter++;
	}
	*teqMat = tMat;
#else
	/* New faster algorithm: binary search
	 */
		tMin = ((double) nCMin) / ((double) freq);
		IF_FAILED_DONE( VnfmBulletB1(
			that,
			S,
			-1e0,	/* fwd ATM */
			tMin,
			freq,
			B,
			B1));
		zBMin = B1[0] / B[0];

		tMax = ((double) nCMax) / ((double) freq);
		IF_FAILED_DONE( VnfmBulletB1(
			that,
			S,
			-1e0,	/* fwd ATM */
			tMax,
			freq,
			B,
			B1));
		zBMax = B1[0] / B[0];

		nIter++;
		nIter++;

		if (zBMax <= zA) {
			*teqMat = tMax;
		} else if (zBMin >= zA) {
			*teqMat = tMin;
		} else {
		    while (nCMin+1 < nCMax) {
			nC = (int) floor((tMin + tMax)*0.5 * freq + 1e-3);

GtoErrMsg("nIter=%3d nCMin=%3d mCMax=%3d nC=%3d\n", nIter, nCMin, nCMax, nC);

			tMat = ((double) nC) / ((double) freq);
			IF_FAILED_DONE( VnfmBulletB1(
				that,
				S,
				-1e0,	/* fwd ATM */
				tMat,
				freq,
				B,
				B1));
			zB = B1[0] / B[0];

			if (zB >= zA) {
				nCMax = nC;
				tMax  = tMat;
				zBMax = zB;
			} else {
				nCMin = nC;
				tMin  = tMat;
				zBMin = zB;
			}
			nIter++;
		    }

		    ASSERT_OR_DONE(nCMin+1 == nCMax);
		    *teqMat = tMin + (zA - zBMin) / (zBMax - zBMin)
					/ ((double) freq);

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



/*f-------------------------------------------------------------
 * Calculate the B1 coefficients for given cash-flows.
 *                                                             
 * <br><br>
 * Computes and returns the price V-coefficients for
 * a general set of fixed cash-flows.
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling VnfmComputeCoeff).
 * Returns 0 iff OK.
 */

int
VnfmCashFlowsB1(
	VnfmData *that,		/* (I) model parameters */
	TDate resetDate,	/* (I) reset effective date */
	long numCf,		/* (I) number of cash flows */
	TDate *cfDates,		/* (I) cash flow dates */
	double *cfPaymts,	/* (I) cash flow amounts */

	double *B,		/* (O) B[i](S,tMat) */
	double *B1)		/* (O) dB/dbeta[i](S,tMat) */
{
static	char	routine[] = "VnfmBondEqmat";
	int	status = FAILURE;
	int	nDim = NDIM;


	int	idxC, idxF, n, m;
	double	S, T,
		A[VNFM_NDIMMAX], 		/* A coefficients */
		A1[VNFM_NDIMMAX], 		/* A1 coefficients */
		z;				/* forward zero */


	/* Clear  */
	for (idxF=0; idxF<=nDim-1; idxF++) {
		B[idxF]  = 0e0;
		B1[idxF] = 0e0;
	}


	/* COmpute time to expiration */
	IF_FAILED_DONE( GtoDayCountFraction(
		REFDATE,
		resetDate,
		GTO_ACT_365F,
		&S));

	for (idxC=0; idxC<numCf; idxC++) {

		/* Coupon cashflow time */
		IF_FAILED_DONE( GtoDayCountFraction(
			REFDATE,
			cfDates[idxC],
			GTO_ACT_365F,
			&T));

		/* compute A coefficients */
		IF_FAILED_DONE( VnfmA1(
			that, S, T, A, A1));

		/* compute forward zero */
		DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);
		DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);
		z = (that->fZero[m] / that->fZero[n])
			*exp( - RATE[m] * (T - ZTT[m])
			      + RATE[n] * (S - ZTT[n]) );


		/* Add */
		for (idxF=0; idxF<=nDim-1; idxF++) {
			B[idxF]   += A[idxF]  * z * cfPaymts[idxC];
			B1[idxF]  += A1[idxF] * z * cfPaymts[idxC];
		}



#ifdef	__DEBUG__
		DrlFPrintf(vnfmFpLog, "VnfmCFB1(%3d): cpn %2d/%2d "
			"S=%7.4f T=%7.4f (%7.4f)  n=%2d m=%2d dcf=%lf "
			"z=%lf zrcc=%lf zr1c=%lf\n",
			__LINE__, idxC+1, numCf, S, T, T-S,
			n, m, cfPaymts[idxC], z, -log(z) / (T-S),
			pow(z, -1./(T-S))-1.);
#endif





	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog,
		"VnfmCFB1(%3d): S=%lf ", __LINE__, S);
	for (idxF=0; idxF<=nDim-1; idxF++)
		DrlFPrintf(vnfmFpLog, " [%2d] %12.8f %12.8f (%12.8f)",
			idxF, B[idxF], B1[idxF],
			B1[idxF]/B[idxF]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}










/*--------------------------------------------------------------
 * Calculate the B1 coefficients for coupon and principal cash-flows.
 *                                                             
 * <br><br>
 * Computes and returns the price V-coefficients for
 * a general set of fixed coupon and principal cash-flows.
 * The coupon flows are being stubbed whereas the principal are not.
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling VnfmComputeCoeff).
 * Returns 0 iff OK.
 */

int
VnfmFlowsB1(
	VnfmData *that,		/* (I) model parameters */
	TDate resetDate,	/* (I) reset effective date */

	int numCoupDates,	/* (I) */
	TDate *coupDates,	/* (I) */
	double *coupPaymts,	/* (I) */

	int numPrinDates,	/* (I) */
	TDate *prinDates,	/* (I) */
	double *prinPaymts,	/* (I) */

	long numExer,		/* (I) number of exercise dates */
	TDate *exerDates,	/* (I) */
	double *exerStrikes,	/* (I) */
	int exerAmer,		/* (I) FALSE=euro */

	double *B,		/* (O) B[i](S,tMat) */
	double *B1)		/* (O) dB/dbeta[i](S,tMat) */
{
static	char	routine[] = "VnfmBondEqmat";
	int	status = FAILURE;
	int	nDim = NDIM;


	int	idxC, idxF,
		numCf;
	double	dcf;
	TDate	effDate,   *cfDates = NULL;
	double	effStrike, *cfPaymts = NULL;

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "\n%s:============ %10s ===================="
		"=======================================================\n\n",
		routine, DrlTDatePrint(NULL, resetDate));
#endif

	/* Clear  */
	for (idxF=0; idxF<=nDim-1; idxF++) {
		B[idxF]  = 0e0;
		B1[idxF] = 0e0;
	}


	/**
	 ** (1) If any, find first exercise date after reset effective date
	 **/
	if (numExer != 0) {
	    if (exerAmer != FALSE) {
		/* American: peform linear interp on strike */
		effDate = MAX(exerDates[0], resetDate);
		IF_FAILED_DONE( DrlTDateLinearInterp1d(
			exerDates,
			exerStrikes,
			numExer,
			effDate,
			&effStrike));
	    } else {
		GtoErrMsg("%s: ME not supported.\n", routine);
		goto done;
	    }
	} else {
	    effDate = resetDate;
	    effStrike = 0e0;
	}


	/**
	 ** (2) set up cash flows.
	 **/

	cfDates = NEW_ARRAY(TDate, numCoupDates + numPrinDates + 2);
	ASSERT_OR_DONE(cfDates != NULL);
	cfPaymts = NEW_ARRAY(double, numCoupDates + numPrinDates + 2);
	ASSERT_OR_DONE(cfPaymts);
	numCf = 0;

	/* Add coupon flows: may need to stub (SIMPLE ACT/ACT)
	 */
	for (idxC=0; idxC<numCoupDates; idxC++) {
	    /* Don't add past or null coupons
	     */
	    if (IS_ALMOST_ZERO(coupPaymts[idxC])) continue;
	    if (coupDates[idxC] < effDate) continue;

	    if (idxC != 0) {
		/* Not first coupon: stub ACT/ACT when applicable
		 */
     		dcf = MIN((double)(coupDates[idxC] - effDate)
               		/(double)(coupDates[idxC] - coupDates[idxC-1]), 1e0);
		cfDates[numCf] = coupDates[idxC];
		cfPaymts[numCf] = coupPaymts[idxC] * dcf;
		numCf++;

	    } else {
		/* If in future, 1st cpn date must be zero,
		 * otherwise bad stub calculaton
		 */
		if (!IS_ALMOST_ZERO(coupPaymts[idxC])) {
		    GtoErrMsg("%s: first coupon date (%s) is non zero (%lf)"
			" and comes after first exer date (%s).\n",
			routine,
			DrlTDatePrint(NULL, coupDates[idxC]),
			coupPaymts[idxC],
			DrlTDatePrint(NULL, effDate));
		    goto done;
		}
	    }
	}


	/* Add principal flows (no stub)
	 */
	for (idxC=0; idxC<numCoupDates; idxC++) {
	    /* Don't add past or null prin
	     */
	    if (IS_ALMOST_ZERO(prinPaymts[idxC])) continue;
	    if (prinDates[idxC] < effDate) continue;

	    cfDates[numCf] = prinDates[idxC];
	    cfPaymts[numCf] = prinPaymts[idxC];
	    numCf++;

	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "%s: Calculated flows at %s:\n",
		routine, DrlTDatePrint(NULL, resetDate));
	for (idxC=0; idxC<numCf; idxC++) {
		DrlFPrintf(vnfmFpLog, "\t[%3d]    %10s    %22.8f\n",
			idxC, DrlTDatePrint(NULL, cfDates[idxC]),
			cfPaymts[idxC]);
	}
#endif


	/**
	 ** (3) Compute V-coefficients
	 **/
	IF_FAILED_DONE( VnfmCashFlowsB1(
		that,
		resetDate,
		numCf,
		cfDates,
		cfPaymts,
		B,
		B1));


	/* made it through OK */
	status = SUCCESS;
done:
	FREE(cfDates);
	FREE(cfPaymts);
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*--------------------------------------------------------------
 *
 */


DLL_EXPORT(long)
VnfmBondEqmat(
	TCurve *zcCurve,	/* (I) zero coupon curve */
	long numCoupDates,	/* (I) */
	TDate *coupDates,	/* (I) */
	double *coupPaymts,	/* (I) */
	long  numPrinDates,	/* (I) */
	TDate *prinDates,	/* (I) */
	double *prinPaymts,	/* (I) */

	double tMatMin,		/* (I) minimum equivalent maturity */
	double tMatMax,		/* (I) minimum equivalent maturity */
	int freq,		/* (I) frequency for equivalent bullet */

	long numResetDates,	/* (I) number of reset dates */
	TDate *resetDates,	/* (I) reset dates */
	double *resetEqmat)	/* (O) reset equivalent maturities */

{
static	char	routine[] = "VnfmBondEqmat";
	int	status = FAILURE;


















	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



