/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file asian_lognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>
#include <complex>

#include "gpclosedforms/gaussian_integrals.h"

#include "gpclosedforms/inverse.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/asian_lognormal.h"

#include "gpbase/numericconstant.h"

CC_USING_NS(std,complex)
CC_USING_NS(std,sqrt)
CC_USING_NS(std,exp)
CC_USING_NS(std,pow)
CC_USING_NS(std,real)

CC_BEGIN_NAMESPACE(ARM)
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions : implements the geman-yor formula
///
//////////////////////////////////////////////////////////////////////////////////////////


complex<double> GemanYorLEcn(double v,double q,complex<double> l)
{
	complex<double> deux(2.,0);
	complex<double> un(1.,0);
	complex<double> q_comp(q,0);
	complex<double> v_comp(v,0);
	complex<double> m=sqrt(deux*l+v*v);
	complex<double> denom(-2*v-2,0);
	denom+=l;
	denom*=l;
	denom*=Gamma(un+m);
	complex<double> numer=	pow((deux*q_comp),(deux-m+v_comp)/deux)*
							Gamma(deux+(m+v_comp)/deux)*
							Hypergeometric1F1(-un+(m-v_comp)/deux,un+m,-un/(deux*q_comp));
	return numer/denom;

}


double GemanYorAsianVanilla_Integrand(double h,double q,double nu,double alpha,double x)
{
	complex<double> hx(0,h*x);
	complex<double> exphx;
	exphx=exp(hx);
	complex<double> alphaplusIx(alpha,x);
	complex<double> alphamoinsIx(alpha,-x);
	return real(exphx*GemanYorLEcn(nu,q,alphaplusIx)+GemanYorLEcn(nu,q,alphamoinsIx)/exphx);
}

double GemanYorAsianVanillaCall(double S,double k,double T,
								   double r,double v,double alpha,int n)
{
	GaussLaguerre_Coefficients c(0,n);
	double Sum=0;
	double x;
	int i;
	double h=T*v*v/4.;
	double q=k*T*v*v/(4.*S);
	double nu=2*r/(v*v)-1.;
	double generalfactor=4*exp(h*alpha-r*T)*S/(T*v*v*2.);
	for( i=0;i<n;i++){
		{
			x=c.get_point(i);
			Sum+= GemanYorAsianVanilla_Integrand(h,q,nu,alpha,x)*exp(x)*c.get_weight(i);
			
		}
	}
	return Sum*generalfactor/ARM_NumericConstants::ARM_PI;
}

CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/