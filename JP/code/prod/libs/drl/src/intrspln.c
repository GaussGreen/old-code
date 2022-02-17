/************************************************************************
 * Module:	DRL - INTER
 * Function:	Linear interpolation
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <float.h>

#include "drlinter.h"		/* prototype consistency */

#undef	HUGE
#define	HUGE	(0.99e11)


/*f---------------------------------------------------------------------
 * Interpolation 1D spline on doubles: initialization.
 * 
 * <br><br>
 * Initialize a spline interpolation of an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$.
 * The values of ${dy/ dx}$ at $x_0$ and $X_{n-1}$
 * should be given (pass 1e12 for natural splines, i.e. ${d^2y/ dx^2}=0$).
 * On exit, <i> y2</i> contains the coefficients
 * to be used in the function <i> SplineInterp</i>
 * to interpolate.
 * If $x&lt;= x^a_0$ computes $y = y^a_0$
 * (and if $x&gt;= x^a_{m-1}$, $y = y^a_{n-1}$).
 * Returns 0 iff successful.
 */

int
DrlSplineInterp1dInit(
	double *x,		/* (I) array of X(i) (i=0,..,n-1) */
	double *y,		/* (I) array of Y(i) (i=0,..,n-1) */
	int n,			/* (I) # of points */
	double yp1,		/* (I) dx/dy at X(0) */
	double ypn,		/* (I) dx/dy at X(n-1) */
	double *y2)		/* (O) should be allocated on entry */
{
static	char	routine[] = "DrlSplineInterp1dInit";
	int	status = FAILURE;
	int	i,k;
	double	p,qn,sig,un,
		*u = NULL;

	u = NEW_ARRAY(double, n);
	if (u == NULL) goto done;

	/* nrc conventions */
	x--; y--; y2--; u--;


	if (yp1 > HUGE)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > HUGE)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];

	/* free tmp memory */
	FREE_ARRAY(++u);

	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f---------------------------------------------------------------------
 * Interpolation 1D spline on doubles: interpolation.
 * 
 * <br><br>
 * Performs a spline interpolation of an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$ at a point $x$.
 * Should be called AFTER <i> SplineInterp1dInit</i> without changing
 * the values of <i> y2</i>.
 * Returns 0 iff successful.
 */

int
DrlSplineInterp1dInterp(
	double *xa,		/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,		/* (I) array of X(i) (i=0,..,n-1) */
	double *y2a,		/* (I) from SplineInit */
	int n,			/* (I) # of points */
	double x,		/* (I) point to intepolate */
	double *y)		/* (O) interpolated value */
{
static	char	routine[] = "DrlSplineInterp1dInterp";
	int klo,khi,k;
	double h,b,a;

	/* nrc2 convetions */
	xa--; ya--; y2a--;

	if (x <= xa[1]) {
		*y = ya[1];
		return(SUCCESS);
	} else if (x >= xa[n]) {
		*y = ya[n];
		return(SUCCESS);
	}


	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) {
		DrlErrMsg("%s: bad xa input.\n", routine);
		return(FAILURE);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+
		(b*b*b-b)*y2a[khi])*(h*h)/6.0;

	return(SUCCESS);
}



/*f---------------------------------------------------------------------
 * Interpolation 1D spline on doubles.
 * 
 * <br><br>
 * Performs a spline interpolation of an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$ at a point $x$.
 * Calls successively <i> SplineInterp1dInit</i> and
 * <i> SplineInterp1dInterp</i>.
 * Returns 0 iff successful.
 */

int
DrlSplineInterp1d(
	double *xa,		/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,		/* (I) array of X(i) (i=0,..,n-1) */
	int n,			/* (I) # of points */
	double yp1,		/* (I) dx/dy at X(0) */
	double ypn,		/* (I) dx/dy at X(n-1) */
	double x,		/* (I) point to intepolate */
	double *y)		/* (O) interpolated value */
{
static	char	routine[] = "DrlSplineInterp1d";
#undef	NMAX
#define	NMAX	32
	double	y2[NMAX];

	if (n >= NMAX) {
		DrlErrMsg("%s: array too large "
			"(max %d, can be recompiled with larger static).\n",
			routine, NMAX);
		return(FAILURE);
	}

	if (DrlSplineInterp1dInit(xa, ya, n, yp1, ypn, y2) != SUCCESS)
		return(FAILURE);

	if (DrlSplineInterp1dInterp(xa, ya, y2, n, x, y) != SUCCESS)
		return(FAILURE);

	return(SUCCESS);
#undef	NMAX
}
