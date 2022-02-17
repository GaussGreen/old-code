/************************************************************************
 * Module:	DRL
 * Submodule:	ROOT
 * File:	
 * Function:	Root Finding
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <float.h>
#include <math.h>

#include "drlmem.h"

#include "drlroot.h"


#define EPS 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/*----------------------------------------------------------------------
 */

DLL_EXPORT(int)
DrlRootFinderBrent(
	double (*func)(double),	/* (I) function to solve */
	double x1,		/* (I) low bound */
	double x2,		/* (I) high bound */
	double tol,		/* (I) tolerance */
	int itmax,		/* (I) max iterations */
	double *xSol)		/* (I) solution */
{
static	char	routine[] = "DrlRootFinderBrent";
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		GtoErrMsg("%s: low  bound f(%6.4e) = %6.4e\n", routine, a, fa);
		GtoErrMsg("%s: high bound f(%6.4e) = %6.4e\n", routine, b, fb);
		GtoErrMsg("%s: root must be bracketed\n", routine);
		return(FAILURE);
	}
	fc=fb;
	for (iter=1;iter<=itmax;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) {
			*xSol = b;
			return(SUCCESS);
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b);
	}

	GtoErrMsg("%s: maximum number of iterations exceeded\n",
		routine);
	return(FAILURE);
}
#undef EPS
#undef	SIGN


/*
 *
 */




#define EPS 3.0e-8

#define	fa	C[1]
#define	fb	C[2]
#define	fc	C[3]
#define	TOL1	C[4]
#define	xm	C[6]
#define	min1	C[7]
#define	min2	C[8]
#define	a	C[10]
#define	b	C[11]
#define	c	C[12]
#define	d	C[13]
#define	e	C[14]
#define	p	C[15]
#define	q	C[16]
#define	r	C[17]
#define	s	C[18]
#define	ITER	(*iter)


/*f---------------------------------------------------------------------
 * Root finding : 1D Brent's in reverse call.
 *                                                         
 * <br><br>
 * Brent's root finding algorithm implemented in reverse procedure.
 * \par\noindent
 * <b> Example:</b>
 * <tt>
 * double X, Y, XMIN=0e0, XMAX=1e0,
 *        TOL = 1e-9, tmp[20];
 * int    retCode,
 *        iter=0, indic;
 * ...
 * X = XMIN;
 * do {
 *        Y =f(X);
 *        retCode = DrlRootFinderBrentRP(&X, Y, XMIN, XMAX,
 *                                       &iter, &indic, 0);
 * } while (retCode == 0);
 * if (indic) {
 *         printf("root found %lf\n", X);
 * } else {
 *         printf("failed: error code %d\n", retCode);
 * }
 * </tt>
 */

DLL_EXPORT(int)
DrlRootFinderBrentRP(
	double *X,	/* (I) root estimate */
	double Y,	/* (I) function value: WARNING: ON FIRST CALL,
			 *     MUST SET TO f(XMIN) */
	double XMIN,	/* (I) low vbound */
	double XMAX,	/* (I) hghh bound */
	int *iter,	/* (I) iteration counter */
	int *indic,	/* (O) TRUE if root found */
	/* DRL_ROOT_TOLX, double tolx,		// tolerance on dx
	 * DRL_ROOT_TOLF, double tolf,		// tolerance on f(x)
	 * DRL_ROOT_MEM, double c[20],		// working memory 
	 * DRL_ROOT_MAXITER, int itermax,	// max # iterations
	 * DRL_ROOT_FPLOG, FILE *fpLog,		// log
	 * DRL_ROOT_NOERRMSG,			// supresses all err msg
	 * 0) (last arg MUST be 0) 
	 */
	...)
{
static	char	routine[] = "DrlRootFinderBrentRP";
	int	status = FAILURE;

	va_list	ap;
static	double	tmpC[32];

	int	what;
	double	*C = tmpC,
		TOL = 1e-8,
		TOLF = DBL_EPSILON;
	int	itmax = 50;
	FILE	*fpLog = NULL;
	int	errMsgFlag = TRUE;

#define	IF_ERRMSG(statement)	if (errMsgFlag) { statement;}


	/* parse extra options */
	C = tmpC;

	va_start(ap, indic);
	while ((what = va_arg(ap, int)) != 0) {
	    switch (what) {
	    case DRL_ROOT_TOLX:
		TOL = va_arg(ap, double);
		break;
	    case DRL_ROOT_TOLF:
		TOLF = va_arg(ap, double);
		break;
	    case DRL_ROOT_MEM:
		C = va_arg(ap, double*);
		break;
	    case DRL_ROOT_MAXITER:
		itmax = va_arg(ap, int); 
		break;
	    case DRL_ROOT_FPLOG:
		fpLog = (FILE*) va_arg(ap, FILE*);
		break;
	    case DRL_ROOT_NOERRMSG:
		errMsgFlag = FALSE;
		break;
	    default:
		GtoErrMsg("%s: bad option.\n", routine);
		goto done;
	    }
	}
	va_end(ap);



	/* start here */

	ITER++;
	*indic = FALSE;

	switch (ITER) {
	case 1:
		a = XMIN;
		b = XMAX;
		fa = Y;
		*X = XMAX;

		if (fabs(fa) <= TOLF) {
			/* root found */
			*X = a;
			*indic = TRUE;
			status = SUCCESS;
			goto done;
		}
		break;

	case 2:
		fb = Y;

		if (fabs(fb) <= TOLF) {
			/* root found */
			*X = b;
			*indic = TRUE;
			status = SUCCESS;
			goto done;
		}

		/*
		 * root not bracketed
		 */
		if (fa * fb > 0e0) {
		    IF_ERRMSG(GtoErrMsg("%s: failed bracketing: "
			"f(%6.4e)=%6.4e  f(%6.4e)=%6.4e.\n",
			routine, a, fa, b, fb));
		    goto done;
		}

		fc=fb;
	default:
		/* Y contains the value of the function */
		fb = Y ;

		if (fb*fc > 0.0) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		TOL1=2.0*EPS*fabs(b)+0.5*TOL;
		xm=0.5*(c-b);
		if (fabs(xm) <= TOL1 || fabs(fb) <= TOLF) {
			/* root found */
			*X = b;
			*indic = TRUE;
			status = SUCCESS;
			goto done;
		}
		if (fabs(e) >= TOL1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(TOL1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > TOL1)
			b += d;
		else
			b += (xm > 0.0 ? fabs(TOL1) : -fabs(TOL1));

		*X = b ;
		/* fb=(*func)(b); */

		/* too many iterations */
		if (ITER > itmax) {
		    IF_ERRMSG(GtoErrMsg("%s: too many itertions "
			"(iter=%d > itmax=%d)\n", routine, ITER, itmax));
		    goto done;
		}
	}


	status = SUCCESS;
done:
	if (fpLog != (FILE*)NULL) {
	    fprintf(fpLog, "%s: it=%3d XMIN=%8.4f XMAX=%8.4f "
		"X=%14.10f Y=%14.10f %s\n",
		routine, ITER, XMIN, XMAX, *X, Y,
		(*indic ? "FOUND" : ""));
	}

	if (status != SUCCESS) {
		IF_ERRMSG(GtoErrMsg("%s: failed.\n", routine));
	}
	return(status || *indic);
}

#undef EPS
