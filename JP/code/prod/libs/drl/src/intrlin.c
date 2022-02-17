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

#include "drlsort.h"		/* DrlDoubleArrayFloorIdx() */
#include "drltime.h"		/* DrlDDate routines */

#include "drlinter.h"		/* prototype consistency */

/*f---------------------------------------------------------------------
 * Interpolation 1D linear on doubles.
 * 
 * <br><br>
 * Linearly interpolates an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$ at a point $x$.
 * If $x&lt;= x^a_0$ computes $y = y^a_0$
 * (and if $x&gt;= x^a_{m-1}$, $y = y^a_{n-1}$).
 * Returns 0 iff successful.
 */

int
DrlLinearInterp1d(
	double *xa,	/* (I) array of X(i) (i=0,..,m-1) */
	double *ya,	/* (I) array of Y(i) (i=0,..,m-1) */
	int m,		/* (I) arrays size */
	double x,	/* (I) point to intepolate */
	double *y)	/* (O) interpolated value */
{
static	char	routine[] = "DrlLinearInterp1d";

	int	i;
	double	u;

	if (m < 1) {
		DrlErrMsg("%s: not enough points (%d).\n",
			routine, m);
		return(1);
	}

	if (m == 1) {
		*y = ya[0];
		return(0);
	}

	if (x <= xa[0]) {
		*y = ya[0];
		return(0);
	}
	if (x >= xa[m-1]) {
		*y = ya[m-1];
		return(0);
	}


	for (i=0; (xa[i] <= x) && (i<=m-1) ; i++);
	i--;
	if (i > m-1) i--;

	u = (x - xa[i]) / (xa[i+1] - xa[i]) ;

	*y =	(1.0-u)	* ya[i  ] +
		u	* ya[i+1] ;

	return(0) ;
}

/*f---------------------------------------------------------------------
 * Interpolation 1D linear on DDate.
 * 
 * <br><br>
 * Linearly interpolates an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$ at a point $x$.
 * If $x&lt;= x^a_0$ computes $y = y^a_0$
 * (and if $x&gt;= x^a_{m-1}$, $y = y^a_{n-1}$).
 * Returns 0 iff successful.
 */

int
DrlDDateLinearInterp1d(
	DDate *xa,	/* (I) array of X(i) (i=0,..,m-1) */
	double *ya,	/* (I) array of Y(i) (i=0,..,m-1) */
	int m,		/* (I) arrays size */
	DDate x,	/* (I) point to intepolate */
	double *y)	/* (O) interpolated value */
{
static	char	routine[] = "DrlLinearInterp1d";

	int	i;
	double	u;

	if (m < 1) {
		DrlErrMsg("%s: not enough points (%d).\n",
			routine, m);
		return(1);
	}

	if (m == 1) {
		*y = ya[0];
		return(0);
	}
#if	!defined(DRL_TDATE_MDY) || defined(DRL_CLIB)
	if (x <= xa[0]) {
#else
	if (DrlDDateDaysDiff(xa[0], x) <= 0) {
#endif
		*y = ya[0];
		return(0);
	}
#if	!defined(DRL_TDATE_MDY) || defined(DRL_CLIB)
	if (x >= xa[m-1]) {
#else
	if (DrlDDateDaysDiff(xa[m-1], x) >= 0) {
#endif
		*y = ya[m-1];
		return(0);
	}


#if	!defined(DRL_TDATE_MDY) || defined(DRL_CLIB)
	for (i=0; (xa[i] <= x) && (i<=m-1) ; i++);
#else
	for (i=0; (DrlDDateDaysDiff(xa[i], x) >= 0) && (i<=m-1) ; i++);
#endif

	i--;
	if (i > m-1) i--;

#ifdef	DRL_CLIB
	u = ((double)x - (double)xa[i]) / ((double)xa[i+1] - (double)xa[i]) ;
#else
	u = ((double) DrlDDateDaysDiff(xa[i], x)) / 
		((double) DrlDDateDaysDiff(xa[i], xa[i+1]));
#endif

	*y =	(1.0-u)	* ya[i  ] +
		u	* ya[i+1] ;

	return(SUCCESS);
}


/*f---------------------------------------------------------------------
 * Interpolation 1D linear on doubles, weights returned.
 * 
 * <br><br>
 * Computes the linear interpolation coefficients:
 * for an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,m-1}$ and a point $x$, returns the
 * indices $j_{\pm}$ and weights $w_{\pm}$ of linear interpolation 
 * defined by 
 * <br>
 * <br> $0 &lt;= j_{+}  - j_{-} &lt;= 1$
 * <br> $x^a_{j_{-}} &lt;= x &lt;= x^a_{j_{+}}$
 * <br> $w_{+} = 1-w_{-} = {
 *       x - x^a_{j_{-}} / x^a_{j_{+}} - x^a_{j_{-}}}$
 * <br>
 * Returns 0 iff successful.
 */

int
DrlLinearInterp1dWeights(
	double *xa,	/* (I) array of X(i) (i=0,..,m-1) */
	int m,		/* (I) arrays size */
	double x,	/* (I) point to intepolate */
	int *jlo,	/* (O) lo index */
	int *jhi,	/* (O) hi index */
	double *wlo,	/* (O) lo weight */
	double *whi)	/* (O) hi weight */
{
static	char	routine[] = "DrlLinearInterp1dWeights";

	if (m < 1) {
		DrlErrMsg("%s: not enough points (%d).\n",
			routine, m);
		return(FAILURE);
	} else if (m == 1) {
		*jlo = *jhi = 0;
		*wlo = 1e0;
		*whi = 0e0;
		return(SUCCESS);
	} else if (x <= xa[0]) {
		*jlo = *jhi = 0;
		*wlo = 1e0;
		*whi = 0e0;
		return(SUCCESS);
	} else if (x >= xa[m-1]) {
		*jlo = *jhi = m-1;
		*wlo = 1e0;
		*whi = 0e0;
		return(SUCCESS);
	} else {
		DrlDoubleArrayFloorIdx(xa, m, x, jlo);
		*jhi = *jlo+1;
		*whi = (x - xa[*jlo]) / (xa[*jhi] - xa[*jlo]);
		*wlo = 1e0 - *whi;
		return(SUCCESS);
	}
}


/*f---------------------------------------------------------------------
 * Interpolation 2D linear on doubles.
 * 
 * <br><br>
 * Linearly interpolates a 2-D discrete function
 * $(y^a_{ij})$ given at $(x^{1a}_i, x^{2a}_j)$ where
 * $i=0,\dots,n_1-1$, $j=0,\dots,n_2-1$
 * at some point $(x^1,x^2)$.
 * On exit, {tt y} contains the interpolated value.
 * Returns 0 iff success.
 */

int
DrlLinearInterp2d(
	double *x1a,		/* (I) array 1st dim */
	double *x2a,		/* (I) array 2nd dim */
	double **ya,		/* (I) array ov values */
	int n1,			/* (I) len of array 1st dim */
	int n2,			/* (I) len of array 2nd dim */
	double x1,		/* (I) desired 1st dim value */
	double x2,		/* (I) desired 1st dim value */
	double *y)		/* (O) interpolated value */
{
static	char	routine[] = "DrlLinearInterp2d";
	int	i, ip1, j, jp1;
	double	u1, u2;

	/* x1-coord */
	if (n1 < 1) {
		DrlErrMsg("%s: not enough points in dim 1 (%d).\n",
			routine, n1);
		return(FAILURE);
	} else if (n1 == 1) {
		i = ip1 = 0;
		u1 = 0;
	} else {
		if (x1 <= x1a[0]) {
			i = ip1 = 0;
			u1 = 0;
		} else if (x1 >= x1a[n1-1]) {
			i = ip1 = n1-1;
			u1 = 0;
		} else {
			for (i=0; (x1a[i] <= x1) && (i<=n1-1) ; i++);
			i--;
			if (i > n1-1) i--;
			ip1 = i+1;
			if (ip1 > n1-1) ip1--;

			if (i != ip1) 
				u1 = (x1 - x1a[i]) / (x1a[i+1] - x1a[i]);
			else
				u1 = 0;
		}
	}



	/* x2-coord */
	if (n2 < 1) {
		DrlErrMsg("%s: not enough points in dim 2 (%d).\n",
			routine, n2);
		return(FAILURE);
	} else if (n2 == 1) {
		j = jp1 = 0;
		u2 = 0;
	} else {
		if (x2 <= x2a[0]) {
			j = jp1 = 0;
			u2 = 0;
		} else if (x2 >= x2a[n2-1]) {
			j = jp1 = n2-1;
			u2 = 0;
		} else {
			for (j=0; (x2a[j] <= x2) && (j<=n2-1) ; j++);
			j--;
			if (j > n2-1) j--;
			jp1 = j+1;
			if (jp1 > n2-1) jp1--;

			if (j != jp1) 
				u2 = (x2 - x2a[j]) / (x2a[j+1] - x2a[j]);
			else
				u2 = 0;
		}
	}





	*y = (1.-u1)*(1.-u2) * ya[i  ][j  ] +
	         u1 *(1.-u2) * ya[ip1][j  ] +
	     (1.-u1)*    u2  * ya[i  ][jp1] +
	         u1 *    u2  * ya[ip1][jp1];

	return(SUCCESS);
}


