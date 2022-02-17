/************************************************************************
 * Module:	DRL
 * Submodule:	SORT
 * File:	
 * Function:	Root Finding
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include "drlsort.h"		/* Prototype consistency */
#include "drltime.h"


/*f----------------------------------------------------------------------
 * Sorting : index search in double array.
 * 
 * <br><br>
 * Given reals y_0<...<y_{n-1} and x, computes
 * <blockquote>
 * j = sup {i such that y_i <= x}.
 * </blockquote>
 * If x < y_0, returns j=0 with error code 1 (0 otherwise).
 */

int
DrlDoubleArrayFloorIdx(double *y, int n, double x, int *j)
{
	register int jl, ju, jm;
	jl = (-1);
	ju = n;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= y[jm])
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
	if (*j == -1) {
		(*j)++;
		return(1);
	}
	return(0);
}


/*f----------------------------------------------------------------------
 * Sorting : index search in DDate array.
 * 
 * <br><br>
 * Given dates y_0<...<y_{n-1} and x, computes
 * <blockquote>
 * j = sup {i such that y_i <= x}.
 * </blockquote>
 * If x < y_0, returns j=0 with error code 1 (0 otherwise).
 */


int
DrlDDateArrayFloorIdx(DDate *y, int n, DDate x, int *j)
{
	register int jl, ju, jm;
	jl = (-1);
	ju = n;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
#ifdef	DRL_CLIB
		if (x >= y[jm])
#else
		if (DrlDDateDaysDiff(y[jm], x) >= 0)
#endif
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
	if (*j == -1) {
		(*j)++;
		return(1);
	}
	return(0);
}


/*f----------------------------------------------------------------------
 * Sorting : index search in double array.
 * 
 * <br><br>
 * Given reals y_0<...<y_{n-1} and x, computes
 * <blockquote>
 * j = inf {i such that y_i &gt;= x}.
 * </blockquote>
 * If x < y_0, returns j=0 with error code 1 (0 otherwise).
 */

int
DrlDoubleArrayCeilIdx(double *y, int n, double x, int *j)
{
	register int jl, ju, jm;
	jl = (-1);
	ju = n;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > y[jm])
			jl=jm;
		else
			ju=jm;
	}
	*j=ju;

	if (*j == n) {
		(*j)--;
		return(1);
	}
	return(0);
}

/*f----------------------------------------------------------------------
 * Sorting : index search in DDate array.
 * 
 * <br><br>
 * Given dates y_0<...<y_{n-1} and x, computes
 * <blockquote>
 * j = inf {i such that y_i &gt;= x}.
 * </blockquote>
 * If x < y_0, returns j=0 with error code 1 (0 otherwise).
 */

int
DrlDDateArrayCeilIdx(DDate *y, int n, DDate x, int *j)
{
	register int jl, ju, jm;
	jl = (-1);
	ju = n;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
#ifdef	DRL_CLIB
		if (x > y[jm])
#else
		if (DrlDDateDaysDiff(y[jm], x) > 0)
#endif
			jl=jm;
		else
			ju=jm;
	}
	*j=ju;

	if (*j == n) {
		(*j)--;
		return(1);
	}
	return(0);
}



/*f----------------------------------------------------------------------
 * Sorting : closest index search in double array.
 * 
 * <br><br>
 * Given reals y_0<...<y_{n-1} and x, computes
 * the index $j$ such that y_j is closest to x, i.e.
 * <blockquote>
 *  |x - y_j| &lt;= inf_k |x - x_k|.$$
 * </blockquote>
 */

int
DrlDoubleArrayClosestIdx(double *y, int n, double x, int *j)
{
	register int jl, ju, jm;
	jl = (-1);
	ju = n;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= y[jm])
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
	if (*j == -1) {
		(*j)++;
		return(0);
	}
	if (*j == n-1)
		return(0);

	if ((x - y[*j]) > 0.5*(y[*j+1] - y[*j]))
		(*j)++;

	return(0);
}




