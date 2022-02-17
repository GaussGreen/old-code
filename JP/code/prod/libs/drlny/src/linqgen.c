/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Linear Equations and Matrix Algebra
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>
#include <string.h>

#include "drllineq.h"		/* prototype consistency */


/*f---------------------------------------------------------------------
 * Vector scalar operations.
 *
 * <br><br>
 * Performs a scalar operation on a vector <br>
 * <br> (v<sub>0</sub>, ..., v<sub>n-1</sub>).<br>
 * The argument <i> operation</i> can be one of the following:<br>
 * <br>
 * <i> "="</i>: sets all ements of {tt v} to <i> scalarArg</i>.<br>
 * <i> "+"</i>: adds <i> scalarArg</i> to all ements of {tt v}.<br>
 * <i> "-"</i>: subtracts <i> scalarArg</i> to all ements of {tt v}.<br>
 * <i> "*"</i>: mutiplies all ements of <i> v</i> by <i> scalarArg</i>.<br>
 * <i> "/"</i>: divieds all ements of <i> v</i> by <i> scalarArg</i>.<br>
 * <br>
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DrlVectScalarOper(
	double *v,		/* (B) vector [0..n-1] */
	int n,			/* (I) number of elements */
	char *operation,	/* (I) operation type (=,+,-,*,/) */
	double scalarArg)	/* (I) argument */
{
static	char	routine[] = "DrlVectScalarOper";
register int	idx;

#undef	DOLOOP
#define	DOLOOP(statement)	{for (idx=0; idx<=n-1; idx++) {statement;}}

	switch (operation[0]) {
	case '=':
		DOLOOP(v[idx]  = scalarArg);
		return(SUCCESS);
	case '+':
		DOLOOP(v[idx] += scalarArg);
		return(SUCCESS);
	case '-':
		DOLOOP(v[idx] -= scalarArg);
		return(SUCCESS);
	case '*':
		DOLOOP(v[idx] *= scalarArg);
		return(SUCCESS);
	case '/':
		DOLOOP(v[idx] /= scalarArg);
		return(SUCCESS);
	default:
		GtoErrMsg("%s: bad operation `%c'.\n", routine);
		return(FAILURE);
	}
}

/*f---------------------------------------------------------------------
 * Vector unary operations.
 *
 * <br><br>
 * Performs a unary vector operation on a vector 
 * <br> (v<sub>0</sub>, ..., v<sub>n-1</sub>).<br>
 * with <i> vectArg</i> as argument (assumed to have same length as <i> v</i>).
 * The argument <i> operation</i> can be one of the follwing: <br><br>
 * <i> "="</i>: copies <i> vectarg</i> to <i> v</i>. <br>
 * <i> "+"</i>: adds <i> vectarg</i> to <i> v</i>. <br>
 * <i> "-"</i>: subtracts <i> vectarg</i> form <i> v</i>. <br>
 * <i> "*"</i>: mutiplies <i> v</i> by <i> vectArg</i> (element by element). <br>
 * <i> "/"</i>: divides <i> v</i> by <i> vectArg</i> (element by element). <br>
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DrlVectUnaryOper(
	double *v,		/* (B) vector [0..n-1] */
	int n,			/* (I) number of elements */
	char *operation,	/* (I) operation type (=,+,-,*,/) */
	double *vectArg)	/* (I) vector argument [0..n-1] */
{
static	char	routine[] = "DrlVectUnaryOper";
register int	idx;

#undef	DOLOOP
#define	DOLOOP(statement)	{for (idx=0; idx<=n-1; idx++) {statement;}}

	switch (operation[0]) {
	case '=':
		DOLOOP(v[idx]  = vectArg[idx]);
		return(SUCCESS);
	case '+':
		DOLOOP(v[idx] += vectArg[idx]);
		return(SUCCESS);
	case '-':
		DOLOOP(v[idx] -= vectArg[idx]);
		return(SUCCESS);
	case '*':
		DOLOOP(v[idx] *= vectArg[idx]);
		return(SUCCESS);
	case '/':
		DOLOOP(v[idx] /= vectArg[idx]);
		return(SUCCESS);
	default:
		GtoErrMsg("%s: bad operation `%c'.\n", routine);
		return(FAILURE);
	}
}

/*f---------------------------------------------------------------------
 * Vector norm.
 *
 * <br><br>
 * Computes and returns the norm of a vector
 * <br> (v<sub>0</sub>, ..., v<sub>n-1</sub>).<br>
 * The string argument <i> normType</i> can be one of the follwing:<br><br>
 * <i> "L1"</i>: computes the $L^1$ norm $\sum |x_i|$. <br>
 * <i> "L2"</i>: computes the $L^2$ norm $\sqrt{\sum |x_i|^2}$. <br>
 * <i> "LINF"</i>: computes the $L^\infty$ norm $\max |x_i|$. <br>
 * <i> "LMIN"</i>: computes smallest element $\min x_i$. <br>
 * <i> "LMAX"</i>: computes largest element $\min x_i$. <br>
 * <br>
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DrlVectNorm(
	double *v,		/* (I) vector [0..n-1] */
	int n,			/* (I) number of elements */
	char *normType,		/* (I) norm type (L1,L2,LINF)*/
	double *retVal)		/* (O) vector norm */
{
static	char	routine[] = "DrlVectNorm";
register int	i;
	double	norm = 0e0;

#define	ISNORMTYPE(s)	(!strcmp(normType, s))


	if (ISNORMTYPE("L1")) {
		for (i=0; i<=n-1; i++) {
		    norm += fabs(v[i]);
		}
	} else if (ISNORMTYPE("L2")) {
		for (i=0; i<=n-1; i++) {
		    norm += v[i]*v[i];
		}
		norm = sqrt(norm);
	} else if (ISNORMTYPE("LINF")) {
		for (i=0; i<=n-1; i++) {
		    norm = MAX(fabs(v[i]), norm);
		}
	} else if (ISNORMTYPE("LMIN")) {
		norm = 1e32;
		for (i=0; i<=n-1; i++) {
		    norm = MIN(v[i], norm);
		}
	} else if (ISNORMTYPE("LMAX")) {
		norm = -1e32;
		for (i=0; i<=n-1; i++) {
		    norm = MAX(v[i], norm);
		}
	} else {
		GtoErrMsg("%s: bad norm type `%s'.\n", routine, normType);
		return(FAILURE);
	}

	*retVal = norm;
	return(SUCCESS);
#undef	ISNORMTYPE
}



