/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "convert.h"
#include "cerror.h"
#include "cgeneral.h"
#include "macros.h"
#include "ldate.h"			/* GTO_ACT_365F */

#include "drlio.h"			/* DrlFPrintf */

#define	_vnfm_SOURCE
#include "vnfmanly.h"

# undef __DEBUG__
#if defined(_WINDLL) 
# undef __DEBUG__
#endif

/*--------------------------------------------------------------
 * Solves in $(X,Y)$ the system of quadratic equations
 * <br>
 * a_1 X^2 + 2 b_1 XY * c_1 Y^2 &=& v_1,\\
 * a_2 X^2 + 2 b_2 XY * c_2 Y^2 &=& v_2.\\
 * <br>
 *
 */

static	int
vnfmBiQuadSolve(
	double a1,		/* (I) */
	double b1,		/* (I) */
	double c1,		/* (I) */
	double v1,		/* (I) */
	double a2,		/* (I) */
	double b2,		/* (I) */
	double c2,		/* (I) */
	double v2, 		/* (I) */
	int *nSol,		/* (O) # of solutions */
	double *Xsol1,		/* (O) 1st set of solutions */
	double *Xsol2)		/* (O) 2nd set of solutions */
{
static	char	routine[] = "BiQuadSolve";
	int	status = FAILURE;
	double	av, bv, cv, ba, ca,
		det, z, z1, z2, ysq1, ysq2, X, Y;

	/* */
#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s:\n", routine);\
	DrlFPrintf(vnfmFpLog, "a1=%lf\tb2=%lf\tc1=%lf\tv1=%lf\n",\
		a1, b1, c1, v1);\
	DrlFPrintf(vnfmFpLog, "a2=%lf\tb2=%lf\tc2=%lf\tv2=%lf\n",\
		a2, b2, c2, v2));
#endif

	*nSol = 0;
	av = v2*a1 - v1*a2;
	bv = v2*b1 - v1*b2;
	cv = v2*c1 - v1*c2;
	det = bv*bv - av * cv;
	if (det < 0e0) {
	    GtoErrMsg("%s: det=%g < 0\n", routine, det);
	    goto done;
	} else if (fabs(det) < 1e-12) {
	    z1 = - bv / av;
	    if (z1 <= 0e0) {
		GtoErrMsg("%s: det=0(zero), z1=%g < 0\n", routine, z1);
		goto done;
	    }
	} else {
	    /* det > 0 */

	    z1 = (- bv + sqrt(det)) / av;
	    z2 = (- bv - sqrt(det)) / av;

	    /* set z1 to be the largest */
	    if (z1 < z2) { z = z1; z1 = z2; z2 = z; }

	    if ((z1 <= 0e0) && (z2 <= 0e0)) {
		/* no >0  root */
		GtoErrMsg("%s: det=%g, but z1=%g, z2=%g < 0\n",
			routine, det, z1, z2);
		goto done;
	    } else if ((z1 > 0e0) && (z2 > 0e0)) {

		ba = b1*a2 - b2 * a1;
		ca = c1*a2 - c2 * a1;
		ysq1 = -av / (2e0*ba*z1 + ca);
		ysq2 = -av / (2e0*ba*z2 + ca);
		if ((ysq1 < 0e0) && (ysq2 < 0e0)) {
		    GtoErrMsg("%s: det=%g, z1=%g, z2=%g > 0\n",
			routine, det, z1, z2);
		    GtoErrMsg("%s: but ysq1=%g < 0, ysq2=%g < 0\n",
			routine, ysq1, ysq2);
		    GtoErrMsg("%s: no solution\n", routine);
		    goto done;
		} else if ((ysq1 >= 0e0) && (ysq2 >= 0e0)) {
		    /*GtoErrMsg("%s: det=%g, z1=%g, z2=%g > 0\n",
			routine, det, z1, z2);
		    GtoErrMsg("%s: but ysq1=%g > 0, ysq2=%g > 0\n",
			routine, ysq1, ysq2);
		    GtoErrMsg("%s: too many solution\n", routine);*/

		    *nSol = 2;
		    Xsol1[1] = sqrt(ysq1);
		    Xsol1[0] = z1 * Xsol1[1];

		    Xsol2[1] = sqrt(ysq2);
		    Xsol2[0] = z2 * Xsol2[1];

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	X = Xsol1[0]; Y = Xsol1[1];
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, " two solutions:\n");\
	DrlFPrintf(vnfmFpLog, " 1st solution: X=%g Y=%g\n", X, Y);\
	DrlFPrintf(vnfmFpLog, "eq1 = %g\n", a1*X*X + 2e0*b1*X*Y + c1*Y*Y - v1);\
	DrlFPrintf(vnfmFpLog, "eq2 = %g\n", a2*X*X + 2e0*b2*X*Y + c2*Y*Y - v2));
	X = Xsol2[0]; Y = Xsol2[1];
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, " 2nd solution: X=%g Y=%g\n", X, Y);\
	DrlFPrintf(vnfmFpLog, "eq1 = %g\n", a1*X*X + 2e0*b1*X*Y + c1*Y*Y - v1);\
	DrlFPrintf(vnfmFpLog, "eq2 = %g\n", a2*X*X + 2e0*b2*X*Y + c2*Y*Y - v2));
#endif
		    status = SUCCESS;
		    goto done;
		} else {
		   /* one root */
		    if (ysq2 > ysq1) {
			z1 = z2;
		    }
		}
	    }
	}

	/* z1 is the good root */
	ba = b1*a2 - b2 * a1;
	ca = c1*a2 - c2 * a1;
	ysq1 = -av / (2e0*ba*z1 + ca);
	if (ysq1 <= DBL_EPSILON) {
		GtoErrMsg("%s: det=%g z1=%g but ysq=%g < 0\n",
			routine, det, z1, ysq1);
		goto done;
	}
	*nSol = 1;
	Xsol1[1] = sqrt(ysq1);
	Xsol1[0] = z1 * Xsol1[1];
	status = SUCCESS;


#if !defined(NO_LOGGING) && defined(__DEBUG__)
	X = Xsol1[0]; Y = Xsol1[1];
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, " solution: X=%g Y=%g\n", X, Y);\
	DrlFPrintf(vnfmFpLog, "eq1 = %g\n", a1*X*X + 2e0*b1*X*Y + c1*Y*Y - v1);\
	DrlFPrintf(vnfmFpLog, "eq2 = %g\n", a2*X*X + 2e0*b2*X*Y + c2*Y*Y - v2));
#endif


	/*made it through OK */
done:
#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "a1=%g\tb1=%g\tc1=%g\tv1=%g\n", a1,b1, c1, v1);\
	DrlFPrintf(vnfmFpLog, "a2=%g\tb2=%g\tc2=%g\tv2=%g\n", a2,b2, c2, v2));
#endif
	if (status != SUCCESS) {
	    GtoErrMsg("a1=%lf\tb2=%lf\tc1=%lf\tv1=%lf\n", a1, b1, c1, v1);
	    GtoErrMsg("a2=%lf\tb2=%lf\tc2=%lf\tv2=%lf\n", a2, b2, c2, v2);
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*f-------------------------------------------------------------
 * Bootstrap 2D spot volatility (obsolete).
 *                                                             
 * <br><br>
 * Low-level routine for calibration of a two time-dependent
 * spot volatilities in an 2-factor model.
 * On entry, "that" contains the model parameters.
 * The routine takes as input two arrays of swaption volatilities
 * of (increasing) expirations corresponding exactly to the timeline
 * of the structure "that" (argument "fDate" of the structure "VnfmData").
 * The calibration starts at timeline index "idxStart" and
 * end at "idxEnd" (both included).
 * The input arrays "tMat1", "freq1", "vol1", "tMat2", "freq2" and "vol2"
 * should have length
 * of at least "idxEnd" and contain respectively the swaptions
 * maturities, frequencies (0,1,2,4,12) and volatilities
 * (the elements with index between 0 and "idxStart"$-1$ are
 * disregarded).
 * The routine returns 0 if the bootstrapping can be done.
 * If a failure occurs (because of a too low input volatility),
 * the bootstrapping stops and the routine exits with an nonzero
 * error code.
 * An error message is written to the error log only
 * if the flag "errMsgFlag" is set to TRUE.
 */

DLL_EXPORT(int)
VnfmVolCalib2VArbitrary(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	TDate *dReset1,		/* (I) 1st array of swapt reset [0..idxEnd] */
	double *tMat1,		/* (I) 1st array of swapt mat [0..idxEnd] */
	int *freq1,		/* (I) 1st array of swapt freq [0..idxEnd] */
	double *vol1,		/* (I) 1st array of swapt vol [0..idxEnd] */
	TDate *dReset2,		/* (I) 2st array of swapt reset [0..idxEnd] */
	double *tMat2,		/* (I) 2nd array of swapt mat [0..idxEnd] */
	int *freq2,		/* (I) 2nd array of swapt freq [0..idxEnd] */
	double *vol2,		/* (I) 2nd array of swaption vol [0..idxEnd] */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue)	/* (O) <0 if calib OK, >0 otherwise */
{
static	char		routine[] = "VnfmVolCalib2VArbitrary";
	int		errCode, status = FAILURE;
	int		nExp, n, nSol, idxSRefRatio,
			idxF, nDim = NDIM;
	double		x1[VNFM_NDIMMAX], x2[VNFM_NDIMMAX],
			tReset1, tReset2,
			a1, b1, c1, v1, a2, b2, c2, v2,
			spv1[2], spv2[2],
			r1, r2, targetRatio = 0.33e0,
			S, dt;

	double		fwdRate;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	double		vm1, vm2;		/* used in debuge mode */
#endif

#define	J00(idx)	(that->fJ[0][idx])
#define	J01(idx)	(that->fJ[1][idx])
#define	J11(idx)	(that->fJ[2][idx])

#define	L00		(that->fBeta[0]+that->fBeta[0])
#define	L01		(that->fBeta[0]+that->fBeta[1])
#define	L11		(that->fBeta[1]+that->fBeta[1])


	/*
	 *
	 */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: start.\n", routine);\
	DrlFPrintf(vnfmFpLog, "\tidxStart=%d (S=%lf)\tidxEnd=%d (S=%lf)\n",\
		idxStart, TT[idxStart], idxEnd, TT[idxEnd]));
#endif

	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, that->fNDates-2);

	/* check if anything to do */
	if (idxStart > idxEnd) {
	    /*GtoErrMsg("%s: idxStart=%d > idxEnd=%d\n",
		routine, idxStart, idxEnd);*/
	    VnfmComputeCoeff(that);
	    status = SUCCESS;
	    goto done;
	}



	/* check two-factor */
	if (that->fNf != 2) {
	    GtoErrMsg("%s: must be 2-factor (got %d)\n", routine, that->fNf);
	    goto done;
	}

	/* main calibration loop */
	for (nExp=idxStart; nExp <= idxEnd; nExp++) {
	    n = nExp-1;

#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: calibrating idx=%2d.\n", routine, nExp));
#endif

	    /* nExp = index of base vol point to be fitted */
	    S = that->fTime[nExp];
	    dt = S - that->fTime[n];

	    /* compute time to rates reset */
	    if (GtoDayCountFraction(REFDATE, dReset1[nExp],
		GTO_ACT_365F, &tReset1) != SUCCESS)
		    goto done;

	    /* compute time to rates reset */
	    if (GtoDayCountFraction(REFDATE, dReset2[nExp],
		GTO_ACT_365F, &tReset2) != SUCCESS)
		    goto done;


	    /* check maturity */
	    if (tMat1[nExp] < 1.92307692e-2) {
		GtoErrMsg("%s: (%s) calibrated rate 1 "
		    "has too low maturity (%lf yrs)\n",
		    routine, GtoFormatDate(that->fDate[nExp]), tMat1[nExp]);
		goto done;
	    }

	    if (tMat2[nExp] < 1.92307692e-2) {
		GtoErrMsg("%s: (%s) calibrated rate 2 "
		    "has too low maturity (%lf yrs)\n",
		    routine, GtoFormatDate(that->fDate[nExp]), tMat2[nExp]);
		goto done;
	    }

	    if (tReset1 < S) {
		GtoErrMsg("%s: reset rate 1 (%s, %lf) < exp rate 1 (%lf)\n",
		    routine, GtoFormatDate(dReset1[nExp]), tReset1, S);
		goto done;
	    }
	    if (tReset2 < S) {
		GtoErrMsg("%s: reset rate 2 (%s, %lf) < exp rate 2 (%lf)\n",
		    routine, GtoFormatDate(dReset2[nExp]), tReset2, S);
		goto done;
	    }



	    /* 
	     * compute the Q coefficients 
	     */
	    if (freq1[nExp] > 0) {
		VnfmB(that, tReset1, tMat1[nExp], freq1[nExp], &fwdRate, x1);
	    } else {
		VnfmQ(that, tReset1, tReset1+tMat1[nExp], &fwdRate, x1);
	    }

	    /* 
	     * Convert to the old percentage B/Q coefficients
	     */
	    for (idxF=0; idxF<=nDim-1;idxF++)
		    x1[idxF] /= fwdRate;


	    if (!IS_ALMOST_ZERO(tReset1-S)) {
		for (idxF=0; idxF<=nDim-1;idxF++)
		    x1[idxF] *= exp(-(tReset1-S)*BETA[idxF]);
	    }

	    if (freq2[nExp] > 0) {
		VnfmB(that, tReset2, tMat2[nExp], freq2[nExp], &fwdRate, x2);
	    } else {
		VnfmQ(that, tReset2, tReset2+tMat2[nExp], &fwdRate, x2);
	    }

	    /* 
	     * Convert to the old percentage B/Q coefficients
	     */
	    for (idxF=0; idxF<=nDim-1;idxF++)
		    x2[idxF] /= fwdRate;

	    
	    if (!IS_ALMOST_ZERO(tReset2-S)) {
		for (idxF=0; idxF<=nDim-1;idxF++)
		    x2[idxF] *= exp(-(tReset2-S)*BETA[idxF]);
	    }

#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog,"\tnExp=%2d  n=%2d  S=%7.4f dt=%7.4f\n",\
		nExp, n, S, dt);\
	    DrlFPrintf(vnfmFpLog,"\tS=%s\tdres1=%s\tdres2=%s\n",\
		GtoFormatDate(that->fDate[nExp]),\
		GtoFormatDate(dReset1[nExp]),\
		GtoFormatDate(dReset2[nExp]));\
	    DrlFPrintf(vnfmFpLog, \
		"\tres1=%lf\tres2=%lf\n", tReset1, tReset2);\
	    DrlFPrintf(vnfmFpLog, \
		"\tmat1=%lf\tmat2=%lf\n", tMat1[nExp], tMat2[nExp]);\
	    DrlFPrintf(vnfmFpLog, \
		"\tvol1=%lf\tvol2=%lf\n", vol1[nExp], vol2[nExp]);\
	    DrlFPrintf(vnfmFpLog, "\tx1[0]=%lf\tx1[1]=%lf\n", x1[0], x1[1]);\
	    DrlFPrintf(vnfmFpLog, "\tx2[0]=%lf\tx2[1]=%lf\n", x2[0], x2[1]));
#endif

	    /* compute coefficients */
	    v1 = SQR(vol1[nExp])*S - (
		SQR(x1[0]) * exp(-dt*2e0*BETA[0]) * J00(n) +
		SQR(x1[1]) * exp(-dt*2e0*BETA[1]) * J11(n) +
		2e0 * x1[0] * x1[1] * exp(-dt*(BETA[0]+BETA[1])) * J01(n));

	    v2 = SQR(vol2[nExp])*S - (
		SQR(x2[0]) * exp(-dt*2e0*BETA[0]) * J00(n) +
		SQR(x2[1]) * exp(-dt*2e0*BETA[1]) * J11(n) +
		2e0 * x2[0] * x2[1] * exp(-dt*(BETA[0]+BETA[1])) * J01(n));


	    a1 = x1[0] * x1[0] * SQR(ALPHA[0]) * L(dt, 2e0*BETA[0]);
	    b1 = x1[0] * x1[1] * ALPHA[0] * ALPHA[1] * RHO[0][n] *
			L(dt, BETA[0]+BETA[1]);
	    c1 = x1[1] * x1[1] * SQR(ALPHA[1]) * L(dt, 2e0*BETA[1]);

	    a2 = x2[0] * x2[0] * SQR(ALPHA[0]) * L(dt, 2e0*BETA[0]);
	    b2 = x2[0] * x2[1] * ALPHA[0] * ALPHA[1] * RHO[0][n] *
			L(dt, BETA[0]+BETA[1]);
	    c2 = x2[1] * x2[1] * SQR(ALPHA[1]) * L(dt, 2e0*BETA[1]);


	    /* solve */
	    errCode = vnfmBiQuadSolve(a1, b1, c1, v1, 
				   a2, b2, c2, v2, 
				   &nSol, spv1, spv2);
	    if ((errMsgFlag != FALSE) && (errCode != 0)) {
		/* abort calibration */
		GtoErrMsg("%s: date %s (%lf) bootstrapping failed\n",
		    routine, GtoFormatDate(that->fDate[nExp]),
		    that->fTime[nExp]);
		GtoErrMsg("\tvol1=%lf vol2=%lf\n",
		    vol1[nExp], vol2[nExp]);
		goto done;
	    }

	    switch (nSol) {
	    case 1:
	    	SIGMA[0][n] = spv1[0];
	    	SIGMA[1][n] = spv1[1];
		break;
	    case 2:
		/* two solutions:
		 * If not 1st step, take closest to previous ratio.
		 * Otherwise, apply the following quite arbitrary rule...
		 */
		if (nExp <= 1) {
		  /* 1st step */

		  if ((ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1])*
		      (ALPHA[0]*spv2[0] + ALPHA[1]*RHO[0][n]*spv2[1]) < 0e0 ) {
		    /*
		     * One inverted factor
		     */
#ifndef	NO_LOGGING
		    GTO_IF_LOGGING(\
		    DrlFPrintf(vnfmFpLog, " One inverted factor.\n");\
		    DrlFPrintf(vnfmFpLog, "s1=%lf\ts2=%lf\ts1+r s2=%lf\n",\
			spv1[0], spv1[1],\
			ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1]);\
		    DrlFPrintf(vnfmFpLog, "s1=%lf\ts2=%lf\ts1+r s2=%lf\n",\
			spv2[0], spv2[1],\
			ALPHA[0]*spv2[0] + ALPHA[1]*RHO[0][n]*spv2[1]));
#endif
		    if ((ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1]) > 0e0) {
			SIGMA[0][n] = spv1[0];
			SIGMA[1][n] = spv1[1];
		    } else {
			SIGMA[0][n] = spv2[0];
			SIGMA[1][n] = spv2[1];
		    }
		  } else if (
		    (ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1] > 0e0) &&
		    (ALPHA[0]*spv2[0] + ALPHA[1]*RHO[0][n]*spv2[1] > 0e0)) {
		    /*
		     * Two inverted factors.
		     * Check ratio closest to target ratio 0.33
		     */
		    /* compute ratios of factors (at t=0) */
		    r1 = (ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1])
			/ (sqrt(1e0-RHO[0][n]*RHO[0][n])*spv1[1]);
		    r2 = (ALPHA[0]*spv2[0] + ALPHA[1]*RHO[0][n]*spv2[1])
			/ (sqrt(1e0-RHO[0][n]*RHO[0][n])*spv1[1]);
		    r1 = fabs(log(r1) - log(targetRatio));
		    r1 = fabs(log(r2) - log(targetRatio));

#ifndef	NO_LOGGING
		    GTO_IF_LOGGING(\
		    DrlFPrintf(vnfmFpLog, " Two inverted factors."\
		    	" Check ratio closest to target ratio 0.33.\n");\
		    DrlFPrintf(vnfmFpLog, "s1=%lf\ts2=%lf\ts1+r s2=%lf\n",\
			spv1[0], spv1[1],\
			ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1]);\
		    DrlFPrintf(vnfmFpLog, "s1=%lf\ts2=%lf\ts1+r s2=%lf\n",\
			spv2[0], spv2[1],\
			ALPHA[0]*spv2[0] + ALPHA[1]*RHO[0][n]*spv2[1]);\
		    DrlFPrintf(vnfmFpLog, \
			"log-d ratio target=%lf\tr1=%lf\tr2=%lf\n",\
			targetRatio, r1, r2));
#endif

		    if (r1 <= r2) {
			SIGMA[0][n] = spv1[0];
			SIGMA[1][n] = spv1[1];
		    } else {
			SIGMA[0][n] = spv2[0];
			SIGMA[1][n] = spv2[1];
		    }

		  } else {
		    GtoErrMsg("%s: date %s (%lf) bootstrapping failed\n",
			routine, GtoFormatDate(that->fDate[nExp]),
			that->fTime[nExp]);
		    GtoErrMsg("  s1=%lf\ts2=%lf\ta1*s1+a2*r*s2=%lf\n",
			spv1[0], spv1[1],
			ALPHA[0]*spv1[0] + ALPHA[1]*RHO[0][n]*spv1[1]);
		    GtoErrMsg("  s1=%lf\ts2=%lf\ta1*s1+a2*r*s2=%lf\n",
			spv2[0], spv2[1],
			ALPHA[0]*spv2[0] + ALPHA[1]*RHO[0][n]*spv2[1]);
		    GtoErrMsg("%s: too many solutions\n", routine);
		    goto done;
		  }
		} else {
		  /* not 1st step */
		  /*idxSRefRatio = n-1;*/
		  idxSRefRatio = 0;

		  r1 = fabs(log(spv1[1]/spv1[0]) - 
			log(SIGMA[1][idxSRefRatio]/SIGMA[0][idxSRefRatio]));
		  r2 = fabs(log(spv2[1]/spv2[0]) - 
			log(SIGMA[1][idxSRefRatio]/SIGMA[0][idxSRefRatio]));

#ifndef	NO_LOGGING
		  GTO_IF_LOGGING(\
		  DrlFPrintf(vnfmFpLog, "closest ratio \tr1=%lf\tr2=%lf\n",\
			targetRatio, r1, r2));
#endif

		  if (r1 <= r2) {
			SIGMA[0][n] = spv1[0];
			SIGMA[1][n] = spv1[1];
		  } else {
			SIGMA[0][n] = spv2[0];
			SIGMA[1][n] = spv2[1];
		  }

		} /* if (nExp ... */
		break;
	    default:
	    	SIGMA[0][n] = 0e0;
	    	SIGMA[1][n] = 0e0;
		break;
	    }
#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, " sel solution: X=%g Y=%g\n",\
			SIGMA[0][n], SIGMA[1][n]));
#endif


	    /* OK ! */

	    /* update coefficients J */
	    n++;
	    dt = that->fTime[n] - that->fTime[n-1];

	    J00(n) =
		exp(-dt*2e0*BETA[0]) * J00(n-1) +
		SQR(SIGMA[0][n-1]) * SQR(ALPHA[0]) * L(dt, 2e0*BETA[0]);
	    J11(n) =
		exp(-dt*2e0*BETA[1]) * J11(n-1) +
		SQR(SIGMA[1][n-1]) * SQR(ALPHA[1]) * L(dt, 2e0*BETA[1]);
	    J01(n) =
		exp(-dt*(BETA[0]+BETA[1])) * J01(n-1) +
		SIGMA[0][n-1] * ALPHA[0] * SIGMA[1][n-1] * ALPHA[1] *
			RHO[0][n-1] * L(dt, BETA[0]+BETA[1]);

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	    GTO_IF_LOGGING( \
	    DrlFPrintf(vnfmFpLog, "\tmod1=%lf\tmod2=%lf\n", \
		sqrt((SQR(x1[0])*JJ[0][n]+2e0*x1[0]*x1[1]*JJ[1][n] + \
			SQR(x1[1])*JJ[2][n]) / S), \
		sqrt((SQR(x2[0])*JJ[0][n]+2e0*x2[0]*x2[1]*JJ[1][n] + \
			SQR(x2[1])*JJ[2][n]) / S)); \
	    VnfmComputeCoeff(that); \
	    VnfmAvgQBVol(that, 0e0, S, S, tMat1[nExp], freq1[nExp], LOGVOL, &vm1); \
	    VnfmAvgQBVol(that, 0e0, S, S, tMat2[nExp], freq2[nExp], LOGVOL, &vm2); \
	    DrlFPrintf(vnfmFpLog, "\td1=%lf\td2=%lf\n",\
			(vm1 - vol1[nExp])/vm1, \
			(vm2 - vol2[nExp])/vm2) \
	    );
#endif

	}


	/* fill the remaining spot vol */
	for (; n<=that->fNDates-1; n++) {
		SIGMA[0][n] = SIGMA[0][n-1];
		SIGMA[1][n] = SIGMA[1][n-1];
	}
	if (resValue != NULL)
		*resValue = 0e0;


#ifndef	NO_LOGGING
	/* test vol values */
	VnfmComputeCoeff(that);
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: calibrated parameters.\n",\
		routine));
	GTO_IF_LOGGING(VnfmFpWrite(that, vnfmFpLog));
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: input/calibrated vol testing.\n", routine);\
	DrlFPrintf(vnfmFpLog, "\t\tvIn1   vOut1  diff1       vIn2   vOUt2  diff2\n");\
	for (nExp=idxStart; nExp <= idxEnd; nExp++) {\
	    /* compute time to rates reset */\
	    if (GtoDayCountFraction(REFDATE, dReset1[nExp],\
		GTO_ACT_365F, &tReset1) != SUCCESS)\
		    goto done;\
\
	    /* compute time to rates reset */\
	    if (GtoDayCountFraction(REFDATE, dReset2[nExp],\
		GTO_ACT_365F, &tReset2) != SUCCESS)\
		    goto done;\
\
	    if (VnfmAvgQBVol(that, 0, TT[nExp],\
		tReset1, tMat1[nExp], freq1[nExp], LOGVOL, &v1) != SUCCESS)\
			goto done;\
\
	    if (VnfmAvgQBVol(that, 0, TT[nExp],\
		tReset2, tMat2[nExp], freq2[nExp], LOGVOL, &v2) != SUCCESS)\
			goto done;\
\
	    DrlFPrintf(vnfmFpLog, "%3d\t%6.4f %6.4f %10.8f  %6.4f %6.4f %10.8f.\n",\
		nExp,\
		vol1[nExp], v1, (vol1[nExp]- v1), \
		vol2[nExp], v2, (vol2[nExp]- v2)\
		);\
	});
#endif



	/* made it through OK */
	status = SUCCESS;
done:
	/* tesing of accuracy */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: done (status = %s).\n",\
		routine, (status == SUCCESS ? "SUCCESS" : "FAILURE")));
#endif
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


