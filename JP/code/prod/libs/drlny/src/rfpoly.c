/************************************************************************
 * Module:      DRL
 * Submodule:   ROOT
 * File:
 * Function:    Root Finding
 * Author:      Taken from Lionnel Pradier
 ************************************************************************/
#include "drlstd.h"             /* platform compatibility */
 
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>
 
 
#define __DEBUG__
#undef  __DEBUG__


/*--------------------------------------------------------------
 * Root finding : 1D polynomials.
 *                                                         
 * <br><br>
 * Newton-Raphson for polynomials.
 * Taken from Lionnel
 */
DLL_EXPORT(int)
DrlNRPoly(
	double	guess,		/* (I) First guess for the root       */
	double	*a,		/* (I) Coefficients of the polynomial */
	int n,			/* (I) Degree of the polynomial       */
	double *retVal)		/* (O) Root (when successful)         */
{
static	char	routine[] = "DrlNRPoly";

const	int	MAXITER = 20;
const	double	MAXXERR = 1e-10;

    double  x;    /* Working value of the root */
    double  dx;   /* Variation of x from one iteration to the other */
    double  P;    /* Value of the polynomial */
    double  dP;   /* Value of the derivative of the polynomial */

    int     i;
    double  j;    /* Iteration index */


    x  = guess;

    for (j = 0; j < MAXITER; j++) 
    {
        P  = a[n];
        dP = 0.;

        /* Evaluate P and dP at the current value x */
        for (i = n - 1; i >= 0; i--)
        {
            dP = P + dP * x;
            P  = a[i] + P * x;
        }

        if (IS_ALMOST_ZERO(dP)) {
		GtoErrMsg("%s: derivative should not vanish.\n",
			routine);
	}

        dx = P/dP;
        x -= dx;    /* Next guess in the Newton-Raphson */

        if (fabs (dx) < MAXXERR) {
		*retVal = x;
		return(SUCCESS);
	}

    }  /* for j */


	GtoErrMsg("%s: maximum number of iterations exceeded.\n",
		routine);
	GtoErrMsg("%s: failed.\n", routine);
	return(FAILURE);
}

