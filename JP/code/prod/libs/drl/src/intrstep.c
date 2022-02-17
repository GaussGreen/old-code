/************************************************************************
 * Module:	DRL - INTER
 * Function:	Linear interpolation
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include "drlsort.h"		/* DrlDDateArrayFloorIdx */
#include "drltime.h"		/* DrlDDate routines */

#include "drlinter.h"		/* prototype consistency */

/*f---------------------------------------------------------------------
 * Interpolation 1D step on doubles.
 * 
 * <br><br>
 * Step interpolation an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,n-1}$ at a point $x$,
 * where $x^1$ and $x$ are dates: computes
 * $j = \sup\{i \mbox{s.t.} x^a_i <= x\}$ and, on exit,
 * $y$ contains $y^a_j$.
 * If $x<x^a_0$, computes $y^a_0$ and returns error code 1.
 */


int
DrlStepInterp1d(
	double *xa,		/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,		/* (I) array of Y(i) (i=0,..,n-1) */
	int n,			/* (I) # of points */
	double x,		/* (I) point to interpolate */
	double *y)		/* (O) interpolated value */
{
	int	j;

	if (x < xa[0]) {
		*y = ya[0];
		return(1);
	}
	DrlDoubleArrayFloorIdx(xa, n, x, &j);
	*y = ya[j];
	return(SUCCESS);
}



/*f---------------------------------------------------------------------
 * Interpolation 1D step on DDate.
 * 
 * <br><br>
 * Step interpolation an array of points
 * $(x^a_i,y^a_i)_{i=0,\dots,n-1}$ at a point $x$,
 * where $x^1$ and $x$ are dates: computes
 * $j = \sup\{i \mbox{s.t.} x^a_i <= x\}$ and, on exit,
 * $y$ contains $y^a_j$.
 * If $x<x^a_0$, computes $y^a_0$ and returns error code 1.
 */


int
DrlDateStepInterp1d(
	DDate *xa,		/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,		/* (I) array of Y(i) (i=0,..,n-1) */
	int n,			/* (I) # of points */
	DDate x,		/* (I) point to interpolate */
	double *y)		/* (O) interpolated value */
{
	int	j;

#ifdef	DRL_CLIB
	if (x < xa[0]) {
#else
	if (DrlDDateDaysDiff(xa[0], x) < 0) {
#endif
		*y = ya[0];
		return(1);
	}
	DrlDDateArrayFloorIdx(xa, n, x, &j);
	*y = ya[j];
	return(SUCCESS);
}






