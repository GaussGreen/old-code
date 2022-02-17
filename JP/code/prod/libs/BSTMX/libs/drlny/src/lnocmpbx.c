/************************************************************************
 * Module:	DRL
 * Submodule:	LNO - Linear and Nonlinear Optimization
 * File:	
 * Function:	Optimisation Routines
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>

#include "drlio.h"
#include "drlmem.h"

#include "drllno.h"		/* prototype consistency */



static	double	uniform_ran1(long *idum, double xlow, double xhigh);

#define	__DEBUG__
#undef	__DEBUG__
#undef	GTO_IF_LOGGING
#if defined(__DEBUG__)
#define GTO_IF_LOGGING(printStatement)    \
{                                         \
    printStatement;                       \
}
#else
#define	GTO_IF_LOGGING(printStatement)
#endif




/*f---------------------------------------------------------------------
 * An implementation of the Complex box method
 * linear and non-linear constraints.\\
 */

static	int
_drlNLinProgComplexBox(
	int NVAR,		/* (I) number of dimensions */
	double *blsc,		/* (I) lbounds state variables */
	double *busc,		/* (I) ubounds state variables */
	int NCNLN,		/* (I) number of nlin constr */
	double *blcnln,		/* (I) lbounds for nlin constr */
	double *bucnln,		/* (I) ubounds for nlin constr */
	TLNOObjConFunc objCon,	/* (I) objective and constraints */
	int ITER,		/* (I) maximum number of iterations */
	double *C,		/* (B) init/final  values of nlin constr */
	double *OBJF,		/* (B) initial value of obj */
	double *X,		/* (B) starting/optimal point */
	void *userData,		/* (I) user data */
	TLNOParams *params)	/* (I) optimization parameters */
{
static	char	routine[] = "_drlNLinProgComplexBox";
	int	status = FAILURE;

	int	n,		/* number of dimensions */
		m,		/* number of constraints */
		N,		/* number of points in the simplex */
		i,
		j,
		k,		/* major iterations counter */
		mu,		/* index for simplex */
		nu,
		needo,		/* indicates if needs objective function */
		*needc = NULL,	/* indicates if constraint needed */
		iterA,		/* iteration for constraints satisfied */
		iterAMax = 100,
		iterC,		/* iteration for constraint satisf */
		iter7,		/* iteration for improvement */
		w,		/* worse vertex index */
		b,		/* bset vertex index */
		retCode;	/* returnd code of oobjective function */


	double	**x = NULL,	/* simplex */
		xinew,		/* temprary variable */
		f,		/* Current objective function */
		*Fx = NULL,	/* objective values at simplex */
		Fw,		/* worse objective value over simplex */
		Fb,		/* best objective value over simplex */
		*xb = NULL,	/* opposite barycenter */
		*xp = NULL,	/* over reflected */
		Fxp,		/* objective values at reflected point */
		*g = NULL,	/* Current value of constraints */
		gj,		/* value of constraint */
		smpstep,	/* initial simplex size */
		xtol,		/* tolerance on variable */
		ftol,		/* tolerance on objective function */
		btol,		/* tolerance on variable bounds */
		alpha;		/* over reflection factor */
	int	restart = 0;

	long	idum = 0L;


	/*----------------------------------------------*
	 * Constants and memory allocation 
	 *----------------------------------------------*/
	alpha = 1.3e0;
	smpstep = 1e-2;
	ftol = params->ftol;
	xtol = params->xtol;
	btol = params->btol;


	n = NVAR;
	N = 2 * NVAR + 1;
	m = NCNLN;

	/* memory allocation */
	x = DrlDoubleMatrAlloc(0, N-1, 0, n-1);
	needc = DrlIntVectAlloc(0, n-1);
	Fx = DrlDoubleVectAlloc(0, N-1);
	g = DrlDoubleVectAlloc(0, m);
	xb = DrlDoubleVectAlloc(0, n-1);
	xp = DrlDoubleVectAlloc(0, n-1);

	ASSERT_OR_DONE(x != NULL);
	ASSERT_OR_DONE(needc != NULL);
	ASSERT_OR_DONE(Fx != NULL);
	ASSERT_OR_DONE(g != NULL);
	ASSERT_OR_DONE(xb != NULL);
	ASSERT_OR_DONE(xp != NULL);

	/*$$$ always compute everything */
	for (i=0; i<n; i++) {
		needc[i] = 1;
	}
	needo = 1;

	GTO_IF_LOGGING(GtoErrMsg("%s:\n", routine);\
	GtoErrMsg("Variables (blsc,busc,init):\n");\
	DrlFPrintDoubleVect(NULL, "\t%g", blsc, NVAR);\
	DrlFPrintDoubleVect(NULL, "\t%g", busc, NVAR);\
	DrlFPrintDoubleVect(NULL, "\t%g", X, NVAR);\
	GtoErrMsg("Constraints (bl,bu):\n");\
	DrlFPrintDoubleVect(NULL, "\t%g", blcnln, NCNLN);\
	DrlFPrintDoubleVect(NULL, "\t%g", bucnln, NCNLN));


	/*----------------------------------------------*
	 * Initialization
	 *----------------------------------------------*/
step0:
	GTO_IF_LOGGING(GtoErrMsg("Step0: restart %d\n", restart));

	/* Set up initial simplex */
	for (i=0; i<n; i++) {
		x[0][i] = X[i];
	}

	for (nu=1; nu<N; nu++) {
	    for (i=0; i<n; i++) {
		/*x[nu][i] = x[0][i] + smpstep * (i == (nu-1));*/
		x[nu][i] = uniform_ran1(
			&idum,
			MAX(x[0][i] - smpstep, blsc[i]+btol),
			MIN(x[0][i] + smpstep, busc[i]-btol));
	    }
	    GTO_IF_LOGGING(GtoErrMsg("Step0: initial nu=%d\n", nu); \
	    DrlFPrintDoubleVect(NULL, "\t%g", x[nu], n));

	    /* Check all constraints satisfied on simplex */
	    for (j=0, iterA=0;
		 j < m;
		 iterA++)
	    {
		/* compute constraint  gj = constraint(j, x[nu]) */
		retCode = (*objCon)(n, x[nu], needo, &f, m, needc, g, userData);
		if (retCode != 0) {
			GtoErrMsg("%s: user defined objective function "
				"returned code %d.\n", routine, retCode);
			goto done;
		}
		gj = g[j];

		GTO_IF_LOGGING( GtoErrMsg("Constraint %d: %12.8f (%lf,%lf)\n",\
			j, gj, blcnln[j], bucnln[j]));

		/* Check constraint satisfied */
		if ((gj < blcnln[j]) || (gj > bucnln[j])) {
		    /* neg constraint: adjust */
	    	    for (i=0; i<n; i++) {
			xinew = 0e0;
			for (mu=0; mu<nu; mu++) {
				xinew += x[mu][i];
			}
			xinew /= ((double) nu);
			xinew += x[nu][i];
			xinew *= 0.5e0;
			x[nu][i] = xinew;
		    }

		} else {
		    /* constr OK: next j */
		    j++;
		}
		if (iterA > iterAMax) {
		    GtoErrMsg("%s: max iter %d to satisfy constraint idx %d "
			"for simplex idx %d.\n", routine,
			iterAMax, j, nu);
		    goto done;
		}
	    }
	}

	/* Compute values at simplex */
	for (nu=0; nu<N; nu++) {
		/* Compute Fx[nu] = objective(x[nu]) */
		retCode = (*objCon)(n, x[nu], needo, &f, m, needc, g, userData);
		if (retCode != 0) {
			GtoErrMsg("%s: user defined objective function "
				"returned code %d.\n", routine, retCode);
			goto done;
		}
		Fx[nu] = f;
	}

	/* set iter count */
	k = 0;
	iter7 = 0;

	GTO_IF_LOGGING(GtoErrMsg("simplex:\n"); \
	DrlFPrintDoubleMatr(NULL, "\t%g", x, N, n));


	/*----------------------------------------------*
	 * Reflection
	 *----------------------------------------------*/
step1:

	GTO_IF_LOGGING(\
	GtoErrMsg("Step1: ------- k = %4d -------\n", k); \
	GtoErrMsg("simplex:\n"); \
	DrlFPrintDoubleMatr(NULL, "\t%g", x, N, n); \
	GtoErrMsg("value:\n"); \
	DrlFPrintDoubleVect(NULL, "\t%g", Fx, N));

	if (k > ITER) {
		GtoErrMsg("%s: too many iterations %d.\n", routine, ITER);
		goto done;
	}

	/* determine worse point on simplex */
	Fw = -DBL_MAX;
	for (nu=0; nu<N; nu++) {
		if (Fx[nu] > Fw) {
			w = nu;
			Fw = Fx[nu];
		}
	}
	GTO_IF_LOGGING( GtoErrMsg("\t worse = %d\n", w));

	/* construct barycenter xb over opposite face */
	for (i=0; i<n; i++) {
		xb[i] = 0e0;
		for (nu=0; nu<N; nu++) {
			if (nu != w) xb[i] += x[nu][i];
		}
		xb[i] /= ((double) N - 1);
	}

	GTO_IF_LOGGING(GtoErrMsg("barycentre:\n"); \
	DrlFPrintDoubleVect(NULL, "\t%g", xb, n));

	/* construct over-reflected xp */
	for (i=0; i<n; i++) {
		xp[i] = xb[i] + alpha * (xb[i] - x[w][i]);
	}

	GTO_IF_LOGGING( GtoErrMsg("reflected:\n"); \
	DrlFPrintDoubleVect(NULL, "\t%g", xp, n));

	iterC = 0;

	/*----------------------------------------------*
	 * Check for constraints
	 *----------------------------------------------*/
step2:
	GTO_IF_LOGGING( GtoErrMsg("Step2:\n"));

	/* Check bounds on variables */
	for (i=0; i<n; i++) {
		xp[i] = MAX(xp[i], blsc[i]+btol);
		xp[i] = MIN(xp[i], busc[i]-btol);
	}


	/* Check implicit constraints */
	j = 0;
	if (m <= 0) goto step6;
	goto step5;


	/*----------------------------------------------*
	 * Check implicit constraints
	 *----------------------------------------------*/
step5:
	GTO_IF_LOGGING( GtoErrMsg("Step5:\n"));

	/* Compute constraint gj = constraint(j, xp) */
	retCode = (*objCon)(n, xp, needo, &f, m, needc, g, userData);
	if (retCode != 0) {
		GtoErrMsg("%s: user defined objective function "
			"returned code %d.\n", routine, retCode);
		goto done;
	}
	gj = g[j];

	GTO_IF_LOGGING( GtoErrMsg("Constraint %d: %12.8f (%lf,%lf)\n", \
		j, gj, blcnln[j], bucnln[j]));


	if ((gj >= blcnln[j]) && (gj <= bucnln[j])) {
		/* Constraint OK */
		goto step6;
	} else {
		if (++iterC < 6) {
			goto step8;
		} else {
			if ((gj >= blcnln[j]-sqrt(DBL_EPSILON)) &&
			    (gj <= bucnln[j]+sqrt(DBL_EPSILON))) {
				goto step9;
			} else {
				/* Too many iteration */
				GtoErrMsg("%s: more that 5 search to satisfy "
					" constraints.\n", routine);
				goto done;
			}
		}
	}


	/*----------------------------------------------*
	 * Implicit constraints loop
	 *----------------------------------------------*/
step6:
	GTO_IF_LOGGING( GtoErrMsg("Step6:\n"));

	if (j < m-1) {
		j += 1;
		goto step5;
	}


	/*----------------------------------------------*
	 * Check for improvement
	 *----------------------------------------------*/
step7:
	GTO_IF_LOGGING( GtoErrMsg("Step7:\n"));

	/* Compute Fxp = objective(xp) */
	retCode = (*objCon)(n, xp, needo, &f, m, needc, g, userData);
	if (retCode != 0) {
		GtoErrMsg("%s: user defined objective function "
			"returned code %d.\n", routine, retCode);
		goto done;
	}
	Fxp = f;

	/* determine if improvement */
	for (nu=0; nu<N; nu++) {
	    if (nu != w) {
		if (Fx[nu] - Fxp > ftol) {
			break;
		}
	    }
	}
	if (nu < N) {
		/* improved: replace worse point by new and next iteration */
		for (i=0; i<n; i++) {
			x[w][i] = xp[i];
			Fx[w] = Fxp;
		}
		iter7 = 0;
		k++;
		goto step1;
	} else {
		/* Failure to improve: if too many tries, done. */
		iter7++;
		if (iter7 > 5) {
			goto step9;
		}
		/* Otherwise try again contract wrt barycenter */
		goto step8;
	}


	/*----------------------------------------------*
	 * Contraction
	 *----------------------------------------------*/
step8:
	GTO_IF_LOGGING( GtoErrMsg("Step8:\n"));

	for (i=0; i<n; i++) {
		xp[i] = (xb[i] + xp[i]) * 0.5e0;
	}
	goto step2;



	/*----------------------------------------------*
	 * Termination
	 *----------------------------------------------*/
step9:
	GTO_IF_LOGGING( GtoErrMsg("Step9: k=%d\n", k));

	/* determine best index in simplex  */
	Fb = DBL_MAX;
	for (nu=0; nu<N; nu++) {
		if (Fx[nu] < Fb) {
			b = nu;
			Fb = Fx[nu];
		}
	}

	/* Returns values */
	for (i=0; i<n; i++) {
		X[i] = x[b][i];
	}
	*OBJF = Fb;

	/* Restart ? */
	if (--restart >= 0)
		goto step0;



	GTO_IF_LOGGING(GtoErrMsg("Optimal vector:\n"); \
	DrlFPrintDoubleVect(NULL, "\t%g", X, n); \
	GtoErrMsg("Optimal value: %12.8f\n", *OBJF)); 



	/*
	 * End
	 */
	status = SUCCESS;
done:
	/* Free memory */
	DrlDoubleMatrFree(x, 0, N-1, 0, n-1);
	DrlIntVectFree(needc, 0, n-1);
	DrlDoubleVectFree(Fx, 0, N-1);
	DrlDoubleVectFree(g, 0, m);
	DrlDoubleVectFree(xb, 0, n-1);
	DrlDoubleVectFree(xp, 0, n-1);


	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*----------------------------------------------------------------------
 *
 */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static	double	uniform_ran1(long *idum, double xlow, double xhigh)
{
	int	j;
	long	k;
	static long	iy = 0;
	static long	iv[NTAB];
	double	temp, retVal;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) 
			*idum = 1;
		else 
			*idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0) 
				*idum += IM;
			if (j < NTAB) 
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0) 
		*idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;

	if ((temp = AM * iy) > RNMX) 
		retVal = (xlow + (xhigh-xlow)*RNMX);
	else 
		retVal = (xlow + (xhigh-xlow)*temp);

	/*GtoErrMsg("uniform_ran1:%ld  %lf.\n", *idum, retVal);*/
	return(retVal);
}


#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*----------------------------------------------------------------------
 *
 */

#define	MAX_CON	32

typedef	struct	{
	void		*userData;
	int		ncnln;
	int		needc[MAX_CON];
	double		cnln[MAX_CON];
	double		*blcnln;
	double		*bucnln;
	TLNOObjConFunc	objConFunc;
} TFeasObjData;

static	int
feasObjFunc(
	int NVAR,		/* (I) number of variables */
	double *X,		/* (I) point where obective valued [NVAR] */
	int NEEDO,		/* (I) indicates if objective  needed */
	double *OBJF,		/* (O) value of objective */
	int NCNLN,		/* (I) number of nlin constr */
	int *NEEDC,		/* (I) array indicate if const needed [NCNLN] */
	double *CNLN,		/* (O) array of values of constr [NCNLN] */
	void *userData) 	/* (I) user data */
{
	TFeasObjData *feasObj = (TFeasObjData*)userData;
	int	j;
	double	dummy;

	if (feasObj->objConFunc(
		NVAR,
		X,
		FALSE,
		&dummy,
		feasObj->ncnln,
		feasObj->needc,
		feasObj->cnln,
		feasObj->userData) != SUCCESS)
			return(FAILURE);


	*OBJF = 0e0;
	for (j=0; j<feasObj->ncnln; j++) {
		*OBJF +=  (feasObj->cnln[j] < feasObj->blcnln[j] ? 
			  feasObj->blcnln[j] - feasObj->cnln[j] : 0e0);
		*OBJF +=  (feasObj->cnln[j] > feasObj->bucnln[j] ? 
			- feasObj->bucnln[j] + feasObj->cnln[j] : 0e0);
	}
	GTO_IF_LOGGING( GtoErrMsg("feasObjFunc: %12.8f\n", *OBJF)); 
	return(SUCCESS);
}



/*f---------------------------------------------------------------------
 * An implementation of the Complex box method
 * linear and non-linear constraints.\\
 */

DLL_EXPORT(int)
DrlNLinProgComplexBox(
	int NVAR,		/* (I) number of dimensions */
	double *blsc,		/* (I) lbounds state variables */
	double *busc,		/* (I) ubounds state variables */
	int NCNLN,		/* (I) number of nlin constr */
	double *blcnln,		/* (I) lbounds for nlin constr */
	double *bucnln,		/* (I) ubounds for nlin constr */
	TLNOObjConFunc objCon,	/* (I) objective and constraints */
	int ITER,		/* (I) maximum number of iterations */
	double *C,		/* (B) initial/final values of nlin constr */
	double *OBJF,		/* (B) initial/final value of obj */
	double *X,		/* (B) starting/optimal  point */
	void *userData,		/* (I) user data */
	TLNOParams *params)	/* (I) optimization parameters */
{
static	char	routine[] = "DrlNLinProgComplexBox";
	int	status = FAILURE;
	TFeasObjData	feasObjData;
	int	j;


	/*
	 * Find feasible point
	 */
	if (NCNLN > MAX_CON) {
		GtoErrMsg("%s: NCNLN %d > MAX_CON %d.\n", routine,
			NCNLN, MAX_CON);
		goto done;
	}

	GTO_IF_LOGGING(GtoErrMsg("%s: FEASIBLE POINT\n", routine));
	feasObjData.userData = userData;
	feasObjData.objConFunc = objCon;
	feasObjData.ncnln = NCNLN;
	feasObjData.blcnln = blcnln;
	feasObjData.bucnln = bucnln;
	for (j=0; j<feasObjData.ncnln; j++) {
		feasObjData.needc[j] = 1;
	}

	if (_drlNLinProgComplexBox(
		NVAR,
		blsc,
		busc,
		0,
		NULL,
		NULL,
		feasObjFunc,
		ITER,
		C,
		OBJF,
		X,
		(void*) &feasObjData,
		params) != SUCCESS) {
			GtoErrMsg("%s: no feasible point found.\n", routine);
			goto done;
	}
	if (*OBJF > 2*DBL_EPSILON) {
		GtoErrMsg("%s: no feasible point found (best %e).\n",
			routine, *OBJF);
		goto done;
	}

	/*
	 * Optimize
	 */
	GTO_IF_LOGGING(GtoErrMsg("%s: OPTIMIZATION\n", routine));
	if (_drlNLinProgComplexBox(
		NVAR,
		blsc,
		busc,
		NCNLN,
		blcnln,
		bucnln,
		objCon,
		ITER,
		C,
		OBJF,
		X,
		userData,
		params) != SUCCESS) {
			goto done;
	}


	/*
	 * End
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}
























