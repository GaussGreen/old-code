/************************************************************************
 * Module:	DRL
 * Submodule:	ROOT
 * File:	
 * Function:	Root Finding
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>

#include "drlmem.h"
#include "drlroot.h"

#if defined(USE_NRC2)
#include "nrc2.h"
#endif

#if defined(UNIX)

#define EPS 1.0e-7
#define	SQR(x)	((x)*(x))
static	double	maxarg1, maxarg2;
#define	FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

static	void	fprintfErrLog(FILE *fp, char *fmt,...);


/*
 * function returning double
 */
typedef	int (*TMDFunc)(
	int n,			/* (I) # dim */
	double *v,		/* (I) point where evaluated [0..n-1] */
	double *f);		/* (O) vector values [0..n-1] */

typedef	double	(*TMDFuncUsrPtr)(
	int n,			/* (I) # dim */
	double *x,		/* (I) point where evaluated */
	double *fvec,		/* (O) on exit, contains values */
	void *usrPtr);		/* (I) usr defined pointer */



/*----------------------------------------------------------------------
 *
 */

#define ALF 1.0e-4
#define TOLX 1.0e-7

static	int
LineSearch(
	int n,
	double xold[],
	double fold,
	double g[],
	double p[],
	double x[],
	double *f,
	double stpmax,
	int *check,
	double *fvec,		/* (I) used to stor vect values */
	TMDFuncUsrPtr func,
	void *usrPtr)		/* (I) usr defined pointer */
{
static	char	routine[] = "LineSearch";

	int	i;
	double	a, alam, alam2, alamin,
	    b, disc,
	    f2, fold2, rhs1, rhs2, slope, sum, temp, 
	    test, tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
	    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
	    slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
	    temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
	    if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;

	for (;;) {
	    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

	    *f = (*func)(n, x, fvec, usrPtr); /*$$$*/

	    if (alam < alamin) {
		for (i=1;i<=n;i++) x[i]=xold[i];
		*check=1;
		return(SUCCESS);
	    } else if (*f <= fold+ALF*alam*slope) return(SUCCESS);
	    else {
		if (alam == 1.0)
		    tmplam = -slope/(2.0*(*f-fold-slope));
		else {
		    rhs1 = *f-fold-alam*slope;
		    rhs2=f2-fold2-alam2*slope;
		    a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
		    b=(-alam2*rhs1/(alam*alam)+
			alam*rhs2/(alam2*alam2))/(alam-alam2);
		    if (a == 0.0) tmplam = -slope/(2.0*b);
		    else {
		    	disc=b*b-3.0*a*slope;
		    	if (disc<0.0) {
			    DrlErrMsg("%s: Roundoff problem\n.", routine);
			    return(FAILURE);
			}
		    	else tmplam=(-b+sqrt(disc))/(3.0*a);
		    }
		    if (tmplam>0.5*alam)
		    	tmplam=0.5*alam;
		}
	    }
	    alam2=alam;
	    f2 = *f;
	    fold2=fold;
	    alam=FMAX(tmplam,0.1*alam);
	}
	return(SUCCESS);
#undef ALF
#undef TOLX
}

/*----------------------------------------------------------------------
 * Evaluates a finite difference jacobian of a 
 * vector function $f:R^n\to R^n$.
 */


static	int
FwdDiffJacobian(
	int n,			/* (I) # dim */
	double *x,		/* (I) point where valued [1..n]  */
	double *fvec,		/* (I) function value [1..n]*/
	double **df,		/* (O) finite diff jacobian [1..n][1..n] */
	TMDFunc vecfunc)	/* (I) */
{
	int	status = FAILURE;
	int	i,j;
	double	h, temp, *f;
#define DEPS 1.0e-4

	if ((f = DrlDoubleVectAlloc(1, n)) == NULL)
	    goto done;

	for (j=1; j<=n; j++) {
	    temp = x[j];
	    h = DEPS*fabs(temp);
	    if (h == 0.0) h = DEPS;
	    x[j] = temp+h;
	    h = x[j]-temp;
#ifndef	_NRC2_ARRAY_CONVENTIONS
	    (*vecfunc)(n,x+1,f+1);
#else
	    (*vecfunc)(n,x,f);
#endif	/*_NRC2_ARRAY_CONVENTIONS*/
	    x[j] = temp;
	    for (i=1; i<=n; i++)
	    	df[i][j] = (f[i]-fvec[i])/h;
	}

	status = SUCCESS;
done:
	DrlDoubleVectFree(f, 1, n);
	return(status);
#undef	DEPS
}


/*----------------------------------------------------------------------
 * Given a vector function $F:R^n\to R^n$,
 * evaluates ${1/ 2} F(x).F(x)$ at the point $x\in R^n$.
 *
 */

static	double
FuncMin(	
	int n,			/* (I) # dim */
	double *x,		/* (I) point where evaluated */
	double *fvec,		/* (O) on exit, contains values */
	void *usrPtr)
{
	int	i;
	double	sum;
	TMDFunc	*vecfunc = (TMDFunc*) usrPtr;

#ifndef	_NRC2_ARRAY_CONVENTIONS
	(*vecfunc)(n, x+1, fvec+1);
#else
	(*vecfunc)(n, x, fvec);
#endif	/*_NRC2_ARRAY_CONVENTIONS*/

	for (sum=0.0,i=1; i<=n; i++)
	    sum += SQR(fvec[i]);
	return(0.5*sum);
}


#define TOLFDEFAULT	1.0e-4
#define TOLX		EPS
#define STPMX		100.0
#define TOLMIN		1.0e-6

#endif	/*defined(UNIX)*/

/*f---------------------------------------------------------------------
 * Root finding : multidimensional using Broyden's method.
 *                                                         
 * <br><br>
 * Multidimensional root finder using Broyden's method. 
 * Performs a search of the zero of a $n$-dimensional
 * vector function <i> vecfunc</i>. On entry, <i> x</i> contains the starting
 * Point of the alogorithm. On successful exit, the zero
 * of the function. <br>
 *
 * The type definition for a n-dimensional function <i>vecfunc</i>
 * of an n-dimensional vector x=(x_0,...,x_{n-1})
 * <blockquote>
 * x = (x_0,...,x_{n-1}) to
 *     f(x)=(f_0(x),...,f_{n-1}(x)).
 * </blockquote>
 * <br><br>
 * Extra options (the last argument of the routine call should always
 * be 0):
 * <br>
 * <br>[DRL\_ROOT\_MAXITER] set the maximum number of
 * itrations to <i> iterMax</i>,
 * <br>[DRL\_ROOT\_FPLOG] prints messages on progress
 * of algorithm on <i> fpLog</i>,
 * <br>[DRL\_ROOT\_TOLF] sets the tolarance on the functional
 * to be <i> tolf</i>, i.e. the routine exits as soon as 
 * $\max_i |f(x_i)| &lt;= \mbox<i> tolf</i>$.
 * <br>[DRL\_ROOT\_BROYDEN\_JACOBIAN] returns the jacobian
 * of the function at the solution (<i> jac</i> should be allocated),
 * <br>
 */

int
DrlRootFindBroyd(
	double *x,				/* (I/O) start/end point */
	int n,					/* (I) number of dim */
	int *check,				/* (I) */
	int (*vecfunc) (int,double*,double*),	/* (I) vector function */
	/*			       Optional arguments:
	 * DRL_ROOT_BROYDEN_JACOBIAN, double **jac,
	 * DRL_ROOT_MAXITER, int iterMax,
	 * DRL_ROOT_FPLOG, FILE *fpLog,
	 * DRL_ROOT_TOLF, double tolf,
	 *			       Last argument must always be 0:
	 * 0,
	 */
	...)
{
static	char	routine[] = "DrlRootFindBroyden";
	int	status = FAILURE;

#if !defined(USE_NRC2)
	DrlErrMsg("%s: routine not available.\n", routine);
	return(FAILURE);
#else

	int	i, its=0, j, k, restrt, sing, skip,
		maxIter = 200;
	double	den, f, fold, stpmax, sum, temp, test,
		tolf = TOLFDEFAULT,	/* tolerance on functional */
		*c = NULL,
		*d = NULL,
		*fvcold = NULL,
		*g = NULL,
		*p = NULL,
		**qt = NULL,
		**r = NULL,		/* jacobian */
		*s = NULL,
		*t = NULL,
		*w = NULL,
		*xold = NULL;
	double	*fvec = NULL;

	va_list	ap;
	int	what, idx1, idx2;
	double	**jacobian = NULL;
	FILE	*fpLog = NULL;

#define FREERETURN {status = SUCCESS; goto done;}

#if !defined(UNIX)
	DrlErrMsg("%s: routine not available.\n", routine);
	*check = 1;
	return(FAILURE);
#endif

#ifndef	_NRC2_ARRAY_CONVENTIONS
	x--;
#endif	/*_NRC2_ARRAY_CONVENTIONS*/


	/*
	 * Parse options
	 */

	va_start(ap, vecfunc);
	while ((what = va_arg(ap, int)) != 0) {
	    switch (what) {
	    case DRL_ROOT_BROYDEN_JACOBIAN:
		jacobian = (double**) va_arg(ap, double**);
		break;
	    case DRL_ROOT_MAXITER:
		maxIter = (int) va_arg(ap, int);
		break;
	    case DRL_ROOT_FPLOG:
		fpLog = (FILE*) va_arg(ap, FILE*);
		break;
	    case DRL_ROOT_TOLF:
		tolf = (double) va_arg(ap, double);
		break;
	    case 0:
		break;
	    default:
		fprintfErrLog(fpLog, "%s: unknown option (%d).\n",
			routine, what);
		goto done;
	    }
	}
	va_end(ap);


	/*
	 *
	 */
	if (fpLog) fprintf(fpLog, "%s: start.\n", routine);

	if ((c = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((d = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((fvcold = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((g = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((p = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((qt = DrlDoubleMatrAlloc(1,n,1,n)) == NULL) goto done;
	if ((r = DrlDoubleMatrAlloc(1,n,1,n)) == NULL) goto done;
	if ((s = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((t = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((w = DrlDoubleVectAlloc(1,n)) == NULL) goto done;
	if ((xold = DrlDoubleVectAlloc(1,n)) == NULL) goto done;

	if ((fvec = DrlDoubleVectAlloc(1,n)) == NULL) goto done;

	/*
	 *
	 */

	f = FuncMin(n, x, fvec, (void*) &vecfunc);


	test=0.0;
	for (i=1;i<=n;i++) {
	    if (fabs(fvec[i]) > test) {
		test=fabs(fvec[i]);
	    }
	}
	if (test<0.01*tolf) FREERETURN
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	restrt=1;


	for (its=1; its<=maxIter; its++) {
	    if (restrt) {
		/* compute jacobian */
		if (fpLog) fprintf(fpLog, "%s: iter # %d calling jac.\n",
			routine, its);
	    	if (FwdDiffJacobian(n, x, fvec, r, vecfunc) != SUCCESS)
			goto done;

		/* copy jacobian if needed */
		if (jacobian != NULL) {
		    for (idx1=0; idx1<=n-1; idx1++)
		    for (idx2=0; idx2<=n-1; idx2++)
			jacobian[idx1][idx2] = r[idx1+1][idx2+1];
		}

#ifdef	UNIX
	    	qrdcmp(r,n,c,d,&sing);
#endif
	    	if (sing) {
		    DrlErrMsg("%s: singular Jacobian:\n", routine);
	    	    DrlErrMsg("NL:%d\nNC:%d\n", n, n);
	    	    for(idx1=0; idx1<=n-1; idx1++) {
			for(idx2=0; idx2<=n-1; idx2++) {
	    	    		DrlErrMsg("\t%lf", r[idx1+1][idx2+1]);
			}
	    	        DrlErrMsg("\n");
	    	    }
		    goto done;
		}
	    	for (i=1;i<=n;i++) {
	    		for (j=1;j<=n;j++) qt[i][j]=0.0;
	    		qt[i][i]=1.0;
	    	}
	    	for (k=1;k<n;k++) {
	    		if (c[k]) {
	    			for (j=1;j<=n;j++) {
	    				sum=0.0;
	    				for (i=k;i<=n;i++)
	    					sum += r[i][k]*qt[i][j];
	    				sum /= c[k];
	    				for (i=k;i<=n;i++)
	    					qt[i][j] -= sum*r[i][k];
	    			}
	    		}
	    	}
	    	for (i=1;i<=n;i++) {
	    		r[i][i]=d[i];
	    		for (j=1;j<i;j++) r[i][j]=0.0;
	    	}
	    } else {
	    	for (i=1;i<=n;i++) s[i]=x[i]-xold[i];
	    	for (i=1;i<=n;i++) {
	    		for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
	    		t[i]=sum;
	    	}
	    	skip=1;
	    	for (i=1;i<=n;i++) {
	    		for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
	    		w[i]=fvec[i]-fvcold[i]-sum;
	    		if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i])))
				skip=0;
	    		else w[i]=0.0;
	    	}
	    	if (!skip) {
	    		for (i=1;i<=n;i++) {
	    			for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
	    			t[i]=sum;
	    		}
	    		for (den=0.0,i=1;i<=n;i++) den += SQR(s[i]);
	    		for (i=1;i<=n;i++) s[i] /= den;
#ifdef	UNIX
	    		qrupdt(r,qt,n,t,s);
#endif
	    		for (i=1;i<=n;i++) {
	    			if (r[i][i] == 0.0) {
				    fprintfErrLog(fpLog, "%s: r singular.\n",
					routine);
				    goto done;
				}
	    			d[i]=r[i][i];
	    		}
	    	}
	    }
	    for (i=1;i<=n;i++) {
	    	for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
	    	g[i]=sum;
	    }
	    for (i=n;i>=1;i--) {
	    	for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
	    	g[i]=sum;
	    }
	    for (i=1;i<=n;i++) {
	    	xold[i]=x[i];
	    	fvcold[i]=fvec[i];
	    }
	    fold=f;
	    for (i=1;i<=n;i++) {
	    	for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
	    	p[i] = -sum;
	    }
#ifdef	UNIX
	    rsolv(r,n,d,p);
#endif

	    if (LineSearch(n,xold,fold,g,p,x,&f,stpmax,check,
			fvec, (TMDFuncUsrPtr) &FuncMin,
			(void*) &vecfunc)
		!= SUCCESS) goto done;

	    test=0.0;
	    for (i=1;i<=n;i++)
	    	if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	    if (test < tolf) {
	    	*check=0;
	    	FREERETURN
	    }
	    if (*check) {
	    	if (restrt) {
		    /* failure and already reinit jacobian => failed */
		    fprintfErrLog(fpLog, "%s: line search failed.\n", routine);
		    goto done;
		} else {
	    	    test=0.0;
	    	    den=FMAX(f,0.5*n);
	    	    for (i=1;i<=n;i++) {
	    		temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	    		if (temp > test) test=temp;
	    	    }
	    	    if (test < TOLMIN) {
		        fprintfErrLog(fpLog, "%s: line search failed.\n",
				routine);
			goto done;
		    } else {
			restrt=1;
		    }
	    	}
	    } else {
	    	restrt=0;
	    	test=0.0;
	    	for (i=1;i<=n;i++) {
	    		temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
	    		if (temp > test) test=temp;
	    	}
	    	if (test < TOLX) FREERETURN
	    }
	}

	DrlErrMsg("%s: maximum # of iterations exceeded (%s).\n",
		routine, maxIter);
	goto done;



	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(fvec,1,n);
	DrlDoubleVectFree(xold,1,n);
	DrlDoubleVectFree(w,1,n);
	DrlDoubleVectFree(t,1,n);
	DrlDoubleVectFree(s,1,n);
	DrlDoubleMatrFree(r,1,n,1,n);
	DrlDoubleMatrFree(qt,1,n,1,n);
	DrlDoubleVectFree(p,1,n);
	DrlDoubleVectFree(g,1,n);
	DrlDoubleVectFree(fvcold,1,n);
	DrlDoubleVectFree(d,1,n);
	DrlDoubleVectFree(c,1,n);


	if (status != SUCCESS) {
	    fprintfErrLog(fpLog, "%s: failed in %d iter.\n",
		routine, its);
	} else {
	    if (fpLog) fprintf(fpLog, "%s: succeeded in %d iter "
		"(check = %d)\n",
		routine, its, *check);
	}
	return(status);
}
#undef EPS
#undef TOLFDEFAULT
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI





static	void
fprintfErrLog(FILE *fp, char *fmt,...)
{
	va_list ap;
	char    buf[256];

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	DrlErrMsg(buf);
	if (fp != (FILE*)NULL) fputs(buf, fp);
	return;

#endif	/*defined(UNIX)*/
}
