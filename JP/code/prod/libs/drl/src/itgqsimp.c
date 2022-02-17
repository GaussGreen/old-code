/************************************************************************
 * Module:	DRL
 * Submodule:	INTEG
 * File:	
 * Function:	Compute 1-D intergral 
 * Author:	From NRC2
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <float.h>
#include <stdio.h>

#include "drlinteg.h"		/* prototype consistency */

#define EPS 1.0e-6
#define JMAX 20

static	double
trapzd(double (*func)(double), double a, double b, int n);


/*f---------------------------------------------------------------------
 * Numerical integration: 1D Simpson.
 *
 * <br><br>
 * Given a function f(x) on a interval [a,b],
 * computes and returns the integral * integral_a^b f(x) dx,
 * using Simpson's rule.
 * See \cite{NRC2}.
 */

int
DrlQSimp(
	double (*func)(double),		/* (I) function to integrate */
	double a,			/* (I) lower bound */
	double b,			/* (I) upper bound */
	double *retVal)			/* (I) integral */
{
static	char	routine[] = "DrlQSimp";

	int j;
	double s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) {
			*retVal = s;
			return(SUCCESS);
		}
		os=s;
		ost=st;
	}

	DrlErrMsg("%s: too many iterations (%d).\n",
		routine, (int) JMAX);
	return(FAILURE);
}



#define FUNC(x) ((*func)(x))

static	double
trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}


