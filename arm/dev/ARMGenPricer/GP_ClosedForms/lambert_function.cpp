/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file lambert_function.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date september 2005
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>
#include <complex>

#include "gpclosedforms/gamma.h"

#include "gpbase/numericconstant.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000



/////////////////////////////////////////////////////////////////////////
///
///    solution in w of  w exp(w) == x
///
/////////////////////////////////////////////////////////////////////////


double lambertfunc (double z)
{
	int n=0;
	double wnext=(z>500)? log(z) : z;
	double w;
	double diff;
	double expw;
	do
	{	
		w=wnext;
		expw=exp(w);
		wnext=w-(w*expw-z)/(expw*(w+1.)-(w*expw-z)*(w+2.)/(2.*w+2.));
		diff=wnext-w;
	}
	while
		((fabs(diff)> ARM_CF_EPS)&&(n<ARM_CF_MAXIT));
	
	if (n > ARM_CF_MAXIT) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"too many iterations in lambertfunc" );
	
	return wnext;
	
}

 
#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/