/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *
 *	\file whittaker.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/whittaker.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/hypergeometric.h"

#include "gpbase/numericconstant.h"

CC_USING_NS(std,log)
CC_USING_NS(std,real)
CC_USING_NS(std,sin)
CC_USING_NS(std,pow)
CC_USING_NS(std,exp)


using std::complex;


#define ARM_CF_EPS 1.0e-12


///
/// Implementation of the whittaker functions (generalization of the bessel functions, an avatar of the confluent hypergeometric functions) 
/// 

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////// Whittaker M  /////////////////////////////////////////////////////////////////
complex<double>  Hypergeometric_Whittaker_M (complex<double>  a, complex<double>  b, complex<double>  z)
{
	
	complex<double>  un(1,0);
	complex<double>  demi(0.5,0);
	complex<double>  deux(2.,0);
//complex<double> toto =  Hypergeometric1F1(demi+b-a,un+deux*b,z)*exp(-z*demi)*exp((demi+b)*log(z));
	return Hypergeometric1F1(demi+b-a,un+deux*b,z)*exp(-z*demi)*exp((demi+b)*log(z));
}


double Hypergeometric_Whittaker_M_regular (double  a, double  b, double  z)
{
	
	complex<double>  un(1,0);
	complex<double>  demi(0.5,0);
	complex<double>  deux(2.,0);
	complex<double> a1(a,0);
	complex<double> b1(b,0);
	complex<double> z1(z,0);
	return real(Hypergeometric1F1(demi+b1-a1,un+deux*b1,z1)*exp(-z1*demi)*pow(z1,demi+b1));
}


// for large k
long double  Hypergeometric_Whittaker_M_asymptotic_a (long double  k, long double  m, long double  z)
{
	long double ap=acos(sqrt(z/(4*k)));
	long double tanap=tan(ap);
	long double tanap2=tanap*tanap;
	long double tanap4=tanap2*tanap2;
	long double q1=(tanap2*6.+5.+tanap4*(m*m*16.-1.)*3.)/(-tanap2*tanap*96.);
	long double q2=tanap2*924.+385.+tanap4*(-m*m*112.+121.)*6.-tanap4*tanap2*(m*m*48.-7.)*36.+tanap4*tanap4*(m*m*16.-1.)*(m*m*16.-9.)*9.;
	q2/=-tanap4*tanap2*18432.;
	long double sin2ap=sin(ap*2.);
	long double cosu=cos(k*(ARM_NumericConstants::ARM_PI-ap*2.+sin2ap)-ARM_NumericConstants::ARM_PI*m+ARM_NumericConstants::ARM_PI_BY_4);
	long double sinu=sin(k*(ARM_NumericConstants::ARM_PI-ap*2.+sin2ap)-ARM_NumericConstants::ARM_PI*m+ARM_NumericConstants::ARM_PI_BY_4);
	long_double r=gammaLD(1+2.*m);
	long_double r2(-ARM_NumericConstants::ARM_HALF_LOG_PI-log(k)*m+m*(m*m*4.-1.)/(24.*k*k),1);
	return (-q1*cosu/k+(q2/(k*k)+1.)*sinu)*ARM_NumericConstants::ARM_SQRT_2*(r*r2).todouble()/sqrt(tanap);

	
}

// for large k
long_double  Hypergeometric_Whittaker_M_asymptotic_a_L (long double  k, long double  m, long double  z)
{
	long double ap=acos(sqrt(z/(4*k)));
	long double tanap=tan(ap);
	long double tanap2=tanap*tanap;
	long double tanap4=tanap2*tanap2;
	long double q1=(tanap2*6.+5.+tanap4*(m*m*16.-1.)*3.)/(-tanap2*tanap*96.);
	long double q2=tanap2*924.+385.+tanap4*(-m*m*112.+121.)*6.-tanap4*tanap2*(m*m*48.-7.)*36.+tanap4*tanap4*(m*m*16.-1.)*(m*m*16.-9.)*9.;
	q2/=-tanap4*tanap2*18432.;
	long double sin2ap=sin(ap*2.);
	long double cosu=cos(k*(ARM_NumericConstants::ARM_PI-ap*2.+sin2ap)-ARM_NumericConstants::ARM_PI*m+ARM_NumericConstants::ARM_PI_BY_4);
	long double sinu=sin(k*(ARM_NumericConstants::ARM_PI-ap*2.+sin2ap)-ARM_NumericConstants::ARM_PI*m+ARM_NumericConstants::ARM_PI_BY_4);
	long double resu1=(-q1*cosu/k+(q2/(k*k)+1.)*sinu)*ARM_NumericConstants::ARM_SQRT_2/sqrt(tanap);
	long_double r=gammaLD(1+2.*m);
	long_double r2(-ARM_NumericConstants::ARM_HALF_LOG_PI-log(k)*m+m*(m*m*4.-1.)/(24.*k*k),1);
	long_double r3=r*r2;
	double resu2=r3.logarithm;
	double exponent=resu2+log(fabs(resu1));
	int signe=resu1/fabs(resu1);
	return long_double(exponent,signe);
	
}

double  Hypergeometric_Whittaker_M (double  a,double  b, double  z)
{

	if(a<80.)
	{
		return Hypergeometric_Whittaker_M_regular(a,b,z);
	}
	else
	{
		return Hypergeometric_Whittaker_M_asymptotic_a(a,b,z);
	}
}

long_double  Hypergeometric_Whittaker_M_L (double  a,double  b, double  z)
{

	if(a<80.)
	{
		return long_double(Hypergeometric_Whittaker_M_regular(a,b,z));
	}
	else
	{
		return Hypergeometric_Whittaker_M_asymptotic_a_L(a,b,z);
	}
}

/////////////////////////// Whittaker W ///////////////////////////////////////////////////////////////////////

/// valid sauf pour b=-0.5
complex<double>  Hypergeometric_Whittaker_W_regular ( complex<double>  a, complex<double>  b, complex<double>  z)
{
	
	complex<double>  un(1,0);
	complex<double>  demi(0.5,0);
	complex<double>  deux(2.,0);

	return Gamma(-deux*b)/Gamma(demi-b-a)*Hypergeometric_Whittaker_M(a,b,z)+Gamma(deux*b)/Gamma(demi+b-a)*Hypergeometric_Whittaker_M(a,-b,z);
}


long double  Hypergeometric_Whittaker_W_regular ( long double  a, long double  b, long double  z)
{
	complex<double> ac(a,0);
	complex<double> bc(b,0);
	complex<double> zc(z,0);
	return real(Hypergeometric_Whittaker_W_regular(ac,bc,zc));
}



// for |z|<1


// for large k
long double  Hypergeometric_Whittaker_W_asymptotic_a ( double  k, double  m, double  z)
{
	long double ap=acos(sqrt(z/(4*k)));
	long double tanap=tan(ap);
	long double tanap2=tanap*tanap;
	long double tanap4=tanap2*tanap2;
	long double q1=(tanap2*6.+5.+tanap4*(m*m*16.-1.)*3.)/(-tanap2*tanap*96.);
	long double q2=tanap2*924.+385.+tanap4*(-m*m*112.+121.)*6.-tanap4*tanap2*(m*m*48.-7.)*36.+tanap4*tanap4*(m*m*16.-1.)*(m*m*16.-9.)*9.;
	q2/=-tanap4*tanap2*18432.;
	long double sin2ap=sin(ap*2.);
	long double cosu=cos(k*(ap*2.-sin2ap)+ARM_NumericConstants::ARM_PI_BY_4);
	long double sinu=sin(k*(ap*2.-sin2ap)+ARM_NumericConstants::ARM_PI_BY_4);
	long double resu1=k*log(k)-k+(m*m*4.-ARM_NumericConstants::ARM_ONE_THIRD)/(k*8.);
	long double resu2=(q1*cosu/k+(q2/(k*k)+1.)*sinu)*2./sqrt(tanap);
	return resu2*exp(resu1);
	
}

// for large k
long_double  Hypergeometric_Whittaker_W_asymptotic_a_L ( double  k, double  m, double  z)
{
	long double ap=acos(sqrt(z/(4*k)));
	long double tanap=tan(ap);
	long double tanap2=tanap*tanap;
	long double tanap4=tanap2*tanap2;
	long double q1=(tanap2*6.+5.+tanap4*(m*m*16.-1.)*3.)/(-tanap2*tanap*96.);
	long double q2=tanap2*924.+385.+tanap4*(-m*m*112.+121.)*6.-tanap4*tanap2*(m*m*48.-7.)*36.+tanap4*tanap4*(m*m*16.-1.)*(m*m*16.-9.)*9.;
	q2/=-tanap4*tanap2*18432.;
	long double sin2ap=sin(ap*2.);
	long double cosu=cos(k*(ap*2.-sin2ap)+ARM_NumericConstants::ARM_PI_BY_4);
	long double sinu=sin(k*(ap*2.-sin2ap)+ARM_NumericConstants::ARM_PI_BY_4);
	long double resu1=k*log(k)-k+(m*m*4.-ARM_NumericConstants::ARM_ONE_THIRD)/(k*8.);
	long double resu2=(q1*cosu/k+(q2/(k*k)+1.)*sinu)*2./sqrt(tanap);
	double exponent=resu1+log(fabs(resu2));
	int signe=resu2/fabs(resu2);
	return long_double(exponent,signe);

	
}


complex<double>  Hypergeometric_Whittaker_W (complex<double>  a, complex<double>  b, complex<double>  z)
{
	return Hypergeometric_Whittaker_W_regular(a,b,z);
}

complex<double>  Hypergeometric_Whittaker_W_3 ( complex<double>  a, complex<double>  b, complex<double>  z);

/// in this implementation it coud be a small jump at k=80 of less than 50bp

double  Hypergeometric_Whittaker_W (double  a,double  b,double  z)
{
	if (b==0.5) 
	{
		complex<double> a1(a,0);
		complex<double> b1(b,0);
		complex<double> z1(z,0);
		return real(Hypergeometric_Whittaker_W_3(a1,b1,z1));
	}
	if(a<80.)
	{
		return Hypergeometric_Whittaker_W_regular(a,b,z);
	}
	else
	{
		return Hypergeometric_Whittaker_W_asymptotic_a(a,b,z);
	}
}
double  Hypergeometric_Whittaker_W_3 (double a, double b, double z);

long_double  Hypergeometric_Whittaker_W_L (double  a,double  b,double  z)
{
	if(a<80.)		/// attention n'est pas valable pour b=0.5 !!!
	{
		return long_double(Hypergeometric_Whittaker_W_regular(a,b,z));
	}
	else
	{
		return Hypergeometric_Whittaker_W_asymptotic_a_L(a,b,z);
	}
}

///////////////////////// Alternatives, may be useful one day ////////////////////

complex<double>  Hypergeometric_Whittaker_W_2 ( complex<double>  a, complex<double>  b, complex<double>  z)
{
	
	complex<double>  un(1,0);
	complex<double>  demi(0.5,0);
	complex<double>  deux(2.,0);
	return HypergeometricU(demi+b-a,un+deux*b,z)*exp(-z*demi)*pow(z,demi+b);
}


double  Hypergeometric_Whittaker_W_3 (double a, double b, double z)
{	
	double Hwm1=Hypergeometric_Whittaker_M(a,b,z);
	double G1=gamma(1.+a+b);
	double G2=gamma(b);
	double Hwm2=Hypergeometric_Whittaker_M(a+1.-b,2.-b,z);
	double G3=gamma(a);
	double G4=gamma(2.-b);
	double Fac=pow(z,1.-b);
	double r1=ARM_NumericConstants::ARM_PI/sin(ARM_NumericConstants::ARM_PI*b);
	double r2=Hwm1/(G1*G2);
	double r3=Fac*Hwm2/(G3*G4);
	double r4=r1*r2-r3;
	return ARM_NumericConstants::ARM_PI/sin(ARM_NumericConstants::ARM_PI*b)*Hwm1/(G1*G2)-Fac*Hwm2/(G3*G4);
}


complex<double>  Hypergeometric_Whittaker_W_3 ( complex<double>  a, complex<double>  b, complex<double>  z)
{	
	complex<double>  un(1,0);
	complex<double>  demi(0.5,0);
	complex<double>  deux(2.,0);
	complex<double> pi(ARM_NumericConstants::ARM_PI,0);
	complex<double>  Hwm1=Hypergeometric_Whittaker_M(a,b,z);
	complex<double> G1=Gamma(un+a-b);
	complex<double> G2=Gamma(b);
	complex<double>  Hwm2=Hypergeometric_Whittaker_M(a+un-b,deux-b,z);
	complex<double> G3=Gamma(a);
	complex<double> G4=Gamma(deux-b);
	complex<double> Fac=pow(z,un-b);
	return (pi/sin(pi*b))*(Hwm1/(G1*G2)-Fac*Hwm2/(G3*G4));
}

/// after Abramowitz pp257 formula 6.1.40 
/// valid only for Real(z)>0

complex<double> GammaLog_2(complex<double> z)
{	
	complex<double> demi(0.5,0);
	complex<double> pi(ARM_NumericConstants::ARM_PI,0);
	complex<double> zn(z);
	complex<double> z2(z*z);
	complex<double> sum=(z-demi)*log(z)-z+demi*ARM_NumericConstants::ARM_LOG_2PI;
	double a[]={0.083333333333333333333,
				-0.0027777777777777777778,
				0.00079365079365079365079,
				-0.00059523809523809523810,
				0.00084175084175084175084,
				-0.0019175269175269175269,
				0.0064102564102564102564,
				-0.029550653594771241830,
				0.17964437236883057316,
				-1.3924322169059011164,
				13.402864044168391994,
				 -156.84828462600201731, 
				2193.1033333333333333,
				 -36108.771253724989357,
				691472.26885131306711,
				-1.52382215394074161922e+7,
				3.829007513914141414141e+8, 
				-1.08822660357843910890e+10,
				3.473202837650022522522e+11, 
				-1.23696021422692744542e+13
				};
	for(int i=0;i<19;i++)
	{	zn/=z2;
		sum+=zn*a[i];	
	}
	return sum;

}


#undef ARM_CF_EPS 


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
