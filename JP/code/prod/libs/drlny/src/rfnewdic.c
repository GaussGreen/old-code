/************************************************************************
 * Module:	DRL
 * Submodule:	ROOT
 * File:	
 * Function:	Root Finding
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>

#include "drlmem.h"

#include "drlroot.h"		/* prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__

#define	X	(*x)		/* current point where f is computed */
#define	FX	(fx)		/* current f(X) */
#define	IND	(*ind)		/* IND = 1 iff root found */
#define	COMPT	(*count) 	/* number of iterations */
#define	TOLX	(tolx)		/* required accuracy for the zero */

#define	FXOLD	c[ 1]
#define	XOLD	c[ 2]
#define	DX	c[ 3]		/* current dx for f' estimation */
#define	XMIN	c[ 4]		/* current low estimate */
#define	XMAX	c[ 5]		/* current high estimate */
#define	FXMIN	c[ 6]		/* current f(XMIN) */
#define	FXMAX	c[ 7]		/* current f(XMAX) */
#define	RSTATE	(*((int*)&(c[15])))	/* current status */

#define	NEW_ITER_MAX	6	/* # of iterations to switch to dichotomy */



#define	RSTAT_INIT  ((int) 0)	/* Newton */
#define	RSTAT_HOUT  ((int) 1)	/* New->Dicho, upbnd now computed */
#define	RSTAT_HBRKT ((int) 2)	/* Dicho, upbnd done, lobnd now comp.*/
#define	RSTAT_LOUT  ((int) 3)	/* New->Dicho, lobnd now computed */
#define	RSTAT_LBRKT ((int) 4)	/* Dicho, lobnd done, upbnd now comp.*/
#define	RSTAT_BRKN  ((int) 5) 	/* Newton w/brkt */
#define	RSTAT_BRKD  ((int) 6)	/* Dicho w/brkt */

static	char	statName[10][10] = {
		"INIT", "HOUT", "HBRKT", "LOUT", "LBRKT",
		"BRKN", "BRKD"};
	



/*f---------------------------------------------------------------------
 * Root finding : Newton/dichotomy 1D.
 *                                                         
 * <br><br>
 * the RF starts with a Newton algorithm. If $f'(x)$ initial is equal 
 * to 0, then $f'(x)$ is estimated on finite differences.
 * If the estimate of the root falls outside the
 * prescribed bounds, or if <i> NEW\_ITER\_MAX</i> steps are
 * done without convergence, the routine switches to a dichotomy
 * algorithm.
 * The algorithm is thus guaranteed to converge if there
 * is a root between the prescribed bounds (dichotomy).
 * In most cases however, one obtains a better convergence speed than
 * with simple dichotomy.
 * <br>
 * <br>[x] On entry initial guess. On exit, best approximation of the root.
 * <br>[xmin] Lower bound for bracketing.
 * <br>[xmax] Upper bound for bracketing.
 * <br>[fx] current value of $f(x)$.
 * <br>[count] counter: should be set to 0 on first call, and unchanged after.
 * <br>[ind] flag that indicates if rootfinder is done,
 *      (1 iff found, 0 if not). Set to 0 on first call. 
 *      Should be checked at each iteration.
 * <br>[tolx] tolerance $|f(x)|&lt;= \mbox<i> tolx</i>$.
 * <br>[c] Working array of size at least 20.
 * <br>
 * <b> WARNING:</b> NOT FULLY DEBUGGED.
 */

DLL_EXPORT(int)
DrlRootFinderNewtonDicho(
	double *x,	/* (I) current estimate of the root */
	double fx,	/* (I) value of function at current estimate */
	double xmin,	/* (I) lower bound */
	double xmax,	/* (I) upper bound */
	int *count,	/* (I) iteration counter */
	int *ind,	/* (I) TRUE if rootfinder is done */
	/* DRL_ROOT_TOLX, double tolx,		// tolerance 
	 * DRL_ROOT_MEM, double c[20],		// working memory 
	 * DRL_ROOT_FPLOG, FILE *fpLog,		// log file 
	 * DRL_ROOT_NEWDIC_DFX, double dfx,	// f'(x) if available
	 * 0) (last arg MUST be 0) 
	 */
	...)
{
static	char	routine[] = "DrlRootFinderNewtonDicho";
	int	status = FAILURE;

	int	what;
	va_list	ap;

	double	*c;
static	double	cstat[32];
	double	tolx = 1e-9, /*DBL_EPSILON;*/
		dfx = 0e0,
		XNEW;
	FILE	*fpLog = NULL;


	/*
	 *
	 */

	c = &cstat[0];

	va_start(ap, ind);
	while ((what = va_arg(ap, int)) != 0) {
	    switch (what) {
	    case DRL_ROOT_NEWDIC_DFX:
		dfx = va_arg(ap, double);
		break;
	    case DRL_ROOT_TOLX:
		tolx = va_arg(ap, double);
		break;
	    case DRL_ROOT_MEM:
		c = va_arg(ap, double*);
		break;
	    case DRL_ROOT_FPLOG:
		fpLog = (FILE*) va_arg(ap, FILE*);
		break;
	    default:
		GtoErrMsg("%s: bad option.\n", routine);
		goto done;
	    }
	}
	va_end(ap);


	if (fpLog != (FILE*)NULL) {
	    fprintf(fpLog, "%s: S [%3d,%s]"
	    	" XMIN=%8.4f FXMIN=%8.4f"
	    	" XMAX=%8.4f FXMAX=%8.4f"
	    	" X=%14.8f FX=%14.8f\n",
		routine, COMPT, statName[RSTATE],
		XMIN, FXMIN, XMAX, FXMAX,
		X, FX);
	}
#if defined(__DEBUG__)
	fprintf(stdout, "%s: Start\n", routine);
	fprintf(stdout, "\tTOLX         %.10e\n", TOLX);
	fprintf(stdout, "\tX/FX:        %12.8f  %12.8f (%.10e)\n", X, FX, FX);
	fprintf(stdout, "\tXOLD/FXOLD:  %12.8f  %12.8f\n", XOLD, FXOLD);
	fprintf(stdout, "\tXMIN/FXMIN:  %12.8f  %12.8f\n", XMIN, FXMIN);
	fprintf(stdout, "\tXMAX/FXMAX:  %12.8f  %12.8f\n", XMAX, FXMAX);
	fprintf(stdout, "\tRSTATE :     %2d %s\n", RSTATE, statName[RSTATE]);
#endif

 	/* If sufficient precision is achieved,
	 * root found stop search
	 */
	if (fabs(FX) <= TOLX) {
		IND = 1;
		status = SUCCESS;
		goto done;
	} 

	/*
	 * First step
	 */
	if (COMPT == 0) {
	    /* Check input values OK */
	    if (xmin > xmax) {
		GtoErrMsg("%s: xmin (%g) > xmax (%g).\n",
			routine, xmin, xmax);
		goto done;
	    }
	    if (X > xmax) {
		GtoErrMsg("%s: X (%g) > xmax (%g).\n",
			routine, X, xmax);
		goto done;
	    }

	    if (X < xmin) {
		GtoErrMsg("%s: X (%g) < xmin (%g).\n",
			routine, X, xmin);
		goto done;
	    }

	    XMIN = xmin;
	    XMAX = xmax;
	    DX = XMAX - XMIN;
	    FXMIN = 0e0;
	    FXMAX = 0e0;


	    /* try first a Newton */
	    XOLD = X;
	    FXOLD = FX;
	    if (IS_ALMOST_ZERO(dfx)) {
		/* If derivative is NA, use a step
		 * of (XMAX - XMIN)*1e-4
		 * to estimate numerically the derivative.
		 */
		X = X + (XMAX - XMIN) * 1e-4;
	    } else {
		X = X - FX / dfx;
	    }
	    RSTATE = RSTAT_INIT;

	    /*
	     * Check if estimate X is inside bounds [XMIN, XMAX].
	     */
	    if (X <= XMIN) {
		X = XMIN ;
		RSTATE = RSTAT_LOUT;
	    }
	    if (X >= XMAX) {
		X = XMAX ;
		RSTATE = RSTAT_HOUT;
	    }
	    COMPT++;

	    /* continue */
	    status = SUCCESS;
	    goto done;
	}


	/*
	 * these steps are performed at 2nd iter or more.
	 */
	switch (RSTATE) {
	case RSTAT_INIT:
		/* Pure Newton step:
		 * Root HAS NOT been bracketed.
		 */
		if (IS_ALMOST_ZERO(dfx)) {
			XNEW = X - FX / ((FX - FXOLD) / (X - XOLD));
		} else {
			XNEW = X - FX / dfx;
		}
		XOLD = X;
		FXOLD = FX;
		X = XNEW;

		/*
		 * If next Newton estimate does not fall within bounds,
		 * or if too many iterations, switch to a dichotomy.
		 * Bracket the root.
		 */
		if ((X >= XMAX) || (COMPT >= NEW_ITER_MAX)) {
			X = XMAX;
			RSTATE = RSTAT_HOUT;
		}
		if (X <= XMIN) {
			X = XMIN;
			RSTATE = RSTAT_LOUT;
		}
		break ;

	case RSTAT_BRKN:
	case RSTAT_BRKD:
		/* 
		 * Newton/Dicho step:
		 * Root HAS been bracketed (FXMIN*FXMAX < 0).
		 */
		if (FXMIN * FX > 0) {
			XMIN	= X;
			FXMIN	= FX; 
			DX	= XMAX - XMIN;
		} else if (FXMAX * FX > 0) {
			XMAX	= X;
			FXMAX	= FX; 
			DX	= XMAX - XMIN;
		}
		XOLD = X;
		FXOLD = FX;

		/* if derivative NA, use secant */
		if (IS_ALMOST_ZERO(dfx)) {
			X = X - FX / ((FX - FXOLD) / (X - XOLD));
		} else {
			X = X - FX / dfx;
		}

		/*
		 * If Newton sends estimate out of bounds or not interesting
		 */
		if ((X >= XMAX) ||
		    (X <= XMIN) ||
		    (!IS_ALMOST_ZERO(dfx) && (fabs(2.0*FX)>fabs(DX*dfx))))
		{
			X = (XMIN + XMAX)*0.5;
			RSTATE = RSTAT_BRKD;
		} else {
			RSTATE = RSTAT_BRKN;
		}
		break;

	case RSTAT_HOUT:
		/* Partial (up) bracketing step:
		 * f(XMAX) has been computed at X=XMAX.
		 * root HAS NOT been bracketed from below.
		 * X = XMAX.
		 */
		FXMAX = FX;
		if (FXMAX * FXOLD < 0e0) {
			/* root is in interval [XOLD, XMAX]
			 * XMIN <- XOLD
			 * XMAX <- X
			 * X    <- (XMAX+XMIN)/2
			 * XOLD unchanged
			 */
			XMIN = XOLD;
			FXMIN = FXOLD;
			XMAX = X;
			FXMAX = FX;
			X = (XOLD + X) * 0.5;

			RSTATE = RSTAT_BRKN;
		} else {
			/* root is NOT in interval [XOLD, XMAX]
			 * try bracket [XMIN,XOLD=XMAX]
			 * XMAX <- XOLD
			 * X    <- XMIN
			 * XOLD unchanged
			 */

			XMAX = XOLD;
			FXMAX = FXOLD;
			X = XMIN;

			RSTATE = RSTAT_LBRKT;
		}
		break;

	case RSTAT_HBRKT:
		/* Partial (up) bracketing step:
		 * f(XMAX) has been computed at X=XMAX.
		 * HBRKT:	root HAS been bracketed from below and
		 *		XMIN, FXMIN contain low bracket estimates.
		 */
		FXMAX = FX;

		if (FXMAX * FXMIN > 0e0) {
			goto bracket_failed;
		}

		/* root is in interval [XOLD, XMAX].
		 * Do Dicho step
		 */
		XMIN = XOLD;
		FXMIN = FXOLD;

		X = (XOLD + X) * 0.5;

		DX = XMAX - XMIN;
		RSTATE = RSTAT_BRKN;

		break;

	case RSTAT_LOUT :
		/* Partial (down) bracketing step:
		 * f(XMIN) has been computed at X=XMIN.
		 * root HAS NOT been upward bracketed.
		 */
		FXMIN = FX;
		if (FXMIN * FXOLD < 0e0) {
			/* root is in interval [XMIN, XOLD] */
			XMAX = XOLD;
			FXMAX = FXOLD;
			XMIN = X;
			FXMIN = FX;
			X = (XOLD + X) * 0.5 ;

			RSTATE = RSTAT_BRKN;
		} else {
		 	/* Do not have full bracketing of the root:
			 * root is NOT in [XMIN,XOLD]
			 * XOLD becomes XMIN (low bracket).
		 	 * The high bracketing needs to be done.
			 * Need to check f(XMAX)
			 */

			XMIN = XOLD;
			FXMIN = FXOLD;
			X = XMAX;

			RSTATE = RSTAT_HBRKT;
			/* XOLD UNCHANGED */
		}
		break ;


	case RSTAT_LBRKT :
		/* Partial (down) bracketing step:
		 * f(XMIN) has been computed at X=XMIN.
		 * root HAS been bracketed from above and
		 * XMAX, FXAX contain high bracket estimates.
		 */
		FXMIN = FX;
		if (FXMIN * FXMAX > 0e0) {
			goto bracket_failed;
		}


		/* root is in interval [XMIN, XOLD].
		 */
		XMAX = XOLD;
		FXMAX = FXOLD;
		X = (XOLD + X) * 0.5 ;

		DX = XMAX - XMIN ;
		RSTATE = RSTAT_BRKN;

		break ;

	default:
		GtoErrMsg("%s: bad state.\n");
		goto done;
	}



	/*
	 * Check if the root lies between the prescribed bounds.
	 * FMIN and FMAX are zero in Newton mode
	 */
bracket_failed:
	if (FXMAX * FXMIN > 0e0) {
		GtoErrMsg("%s: bounds do not bracket root "
			"fxmin (%g) * fxmax(%g)  > 0\n",
			routine, FXMIN, FXMAX);
		goto done;
	}

	COMPT++;


	/* made it through */
	status = SUCCESS;
done:
	if (fpLog != (FILE*)NULL) {
	    fprintf(fpLog, "%s: E [%3d,%s]"
	    	" XMIN=%8.4f FXMIN=%8.4f"
	    	" XMAX=%8.4f FXMAX=%8.4f"
	    	" X=%14.8f FX=%14.8f\n",
		routine, COMPT, statName[RSTATE],
		XMIN, FXMIN, XMAX, FXMAX,
		X, FX);
	}
#ifdef	__DEBUG__
	fprintf(stdout, "%s: End\n", routine);
	fprintf(stdout, "\tX/FX:        %12.8f  %12.8f (%.10e)\n", X, FX, FX);
	fprintf(stdout, "\tXOLD/FXOLD:  %12.8f  %12.8f\n", XOLD, FXOLD);
	fprintf(stdout, "\tXMIN/FXMIN:  %12.8f  %12.8f\n", XMIN, FXMIN);
	fprintf(stdout, "\tXMAX/FXMAX:  %12.8f  %12.8f\n", XMAX, FXMAX);
	fprintf(stdout, "\tRSTATE :     %2d %s\n",
		RSTATE, statName[RSTATE]);
#endif

	return(status);
}



