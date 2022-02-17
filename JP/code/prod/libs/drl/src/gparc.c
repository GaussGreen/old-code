/************************************************************************
 * Module:	DRL
 * Submodule:	FPARS - Formula Parsing
 * File:	
 * Function:	Formula Parsing
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "drlgpars.h"			/* prototype consistency */


/*f---------------------------------------------------------------------
 * Dynamic Formula Parser.
 * <br>
 * <br>[n] the number of variables in the formula.
 * <br>[x] array of length "n" of the values to be substituted.
 * <br>[formula] a formula of $n$ variables labelled
 *   "x0", "x1",\dots,"xN-1" that contains
 *    <br>
 *    <br>  numerical constants, such as <i> 1.0</i> or <i> 4.56e-2</i>,
 *    <br>  basic arithmetic operations $+$, $-$, $*$ and $/$
 *    <br>  exponent <i> EXP</i>, logarithm <i> LOG</i>
 *    <br>  <i> MIN(x,y)</i> and <i> MAX(x,y)</i>
 *    <br>  the <i> IF(cond,x,y)</i> statement
 *    <br>
 * <br>[retVal] on successful exit, contains the value obtained
 *    by evaluating the formula by replacing the formal arguments
 *    "x0", "x1",\dots,"xN-1" by the values given in the input array
 *   "x".
 * <br>
 * <b> Available on UNIX only</b>.\\
 */

int
DrlGParEval(
	int n,			/* (I) # of args */
	double *x,		/* (I) array of args [0..n-1] */
	char *formula,		/* (I) string formula */
	double *retVal)		/* (O) return value */
{
#if defined(UNIX) || defined(BISON)
extern	int	_gParEval(int, double*, char*, double*);

	return _gParEval(n, x, formula, retVal);
#else
static	char	routine[] = "DrlGParEval";
	DrlErrMsg("%s: routine N/A.\n", routine);
	return(FAILURE);
#endif
}

