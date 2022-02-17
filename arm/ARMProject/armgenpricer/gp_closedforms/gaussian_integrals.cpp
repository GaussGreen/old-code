/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_integrals.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>
#include <complex>
#include <vector>

#include "gpclosedforms/long_double.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/erf.h"
#include "gpclosedforms/normal.h"

#include "gpbase/numericconstant.h"

#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 4000


long double atanh(long double y)
{
	return 0.5*log((1+y)/(1-y));
}


complex<long double> atanh(complex<long double> y)
{
	complex<long double> half(0.5,0.0);
	complex<long double> un(1.0,0.0);
	return half*log((un+y)/(un-y));
}

complex<long double> arctan(complex<long double> y)
{
	complex<long double> I(0.0,1.0);
	return -I*atanh(I*y);
}



long double acoth(long double y)
{
	return 0.5*log((1+y)/(y-1));
}

inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;} 

///////////////////////////////////////////////////////////////////////
///
/// les coefficients (abcisses et poids, commencent a i=0 
///
///////////////////////////////////////////////////////////////////////


Orthogonal_Polynomial_Coefficients::Orthogonal_Polynomial_Coefficients(int n)
:	weigths(n), points(n)
{

}



Orthogonal_Polynomial_Coefficients* Orthogonal_Polynomial_Coefficients::clone()
{
	int n=get_order();
	Orthogonal_Polynomial_Coefficients* newcoefs=new Orthogonal_Polynomial_Coefficients(n);
	for(int i=0;i<n;i++)
	{		
		newcoefs->set_points(i,points[i]) ;
		newcoefs->set_weigths(i,weigths[i]) ;
	}
	return newcoefs;
}

/////////////////////////////////////////////////////////////////////////
///
/// Gauss hermite coefficients: can be used to integrales like
/// Integral{-infinity, +infinity} of exp(-x^2)*f(x)
///
/////////////////////////////////////////////////////////////////////////



GaussHermite_Coefficients::GaussHermite_Coefficients(int n)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	GaussHermite_Coefficients_repository* repo=GaussHermite_Coefficients_repository::get_instance();
	Orthogonal_Polynomial_Coefficients* coefp=repo->get_coefficients(n);
	if(!coefp)
	{
		int i,its,j,m;
		long double p1,p2,p3,pp,z,z1;
		m=(n+1)/2;
		for (i=1;i<=m;i++) 
		{
			if (i == 1) {
// FIXMEFRED: mig.vc8 (23/05/2007 14:26:39):cast
				z=sqrt((long double)(2*n+1))-1.85575*pow((long double)(2*n+1),static_cast<long double>(-0.16667));
			} else if (i == 2) {
// FIXMEFRED: mig.vc8 (23/05/2007 14:26:42):cast
				z -= 1.14*pow((long double)n,static_cast<long double>(0.426))/z;
			} else if (i == 3) {
				z=1.86*z-0.86*points[0]; 
			} else if (i == 4) {
				z=1.91*z-0.91*points[1]; 
			} else {
				z=2.0*z-points[i-3];  
			}
			for (its=1;its<=ARM_CF_MAXIT;its++) {
				p1=ARM_NumericConstants::ARM_INVSQRTSQRTPI;
				p2=0.0;
				for (j=1;j<=n;j++) {
					p3=p2;
					p2=p1;
					p1=z*sqrt(2.0/j)*p2-sqrt(((long double)(j-1))/j)*p3;
				}
				pp=sqrt((long double)2*n)*p2;
				z1=z;
				z=z1-p1/pp;
				if (fabs(z-z1) <= ARM_CF_EPS) break;
			}
			if (its > ARM_CF_MAXIT) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"too many iterations in GaussHermite_Coefficients" );
			
			points[i-1]=z;		
			points[n-i] = -z; 
			weigths[i-1]=2.0/(pp*pp);  
			weigths[n-i]=weigths[i-1];  
		}
		repo->get_repository()->insert(make_pair(n,clone())); // we store it now in  the repository
	}
	else
	{
		for(int i=0;i<n;i++)		// the order n has already benn computed and is stored in the repository
		{
			points[i]=coefp->get_point(i);
			weigths[i]=coefp->get_weight(i);
		}
	}
}

/////////////////////////////////////////////////////////////////////////
///
/// Reduced Gauss hermite coefficients: can be used to integrales like
/// Integral{-infinity, +infinity} of exp(-x^2/2)/sqrt(2pi) *f(x)
///
/////////////////////////////////////////////////////////////////////////



ReducedGaussHermite_Coefficients::ReducedGaussHermite_Coefficients(int n)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	GaussHermite_Coefficients coefs(n);
	int i;
	for (i=0;i<n;i++)
	{	
		points[i]=coefs.get_point(i)*ARM_NumericConstants::ARM_SQRT_2;
		weigths[i]=coefs.get_weight(i)*ARM_NumericConstants::ARM_INVSQRTPI;
	}
}
/////////////////////////////////////////////////////////////////////////
///
///  ReducedGaussHermite_Coefficients : can be used to integrales like
///			Integral{0, +infinity} of x^alpha*exp(-x)*f(x)
///		simulate ReducedGaussHermite_Coefficients, but within boundaries (really based on legandres coefficients)
///
/////////////////////////////////////////////////////////////////////////
ReducedGaussHermite_Coefficients::ReducedGaussHermite_Coefficients(GaussLegendre_Coefficients* glcoeffs_ptr, double a,double b):
Orthogonal_Polynomial_Coefficients(glcoeffs_ptr->get_order())
{
	int n=glcoeffs_ptr->get_order(),i;
		for (i=0;i<n;i++) {
			points[i]=((glcoeffs_ptr->get_point(i))*(b-a)+a+b)/2.;
			weigths[i]=(glcoeffs_ptr->get_weight(i))*exp(-(points[i]*points[i])/2.)*(b-a)/(ARM_NumericConstants::ARM_SQRT_8_PI);
		}

}

/////////////////////////////////////////////////////////////////////////
///
///  Gauss Laguerre coefficients: can be used to integrales like
/// Integral{0, +infinity} of x^alpha*exp(-x)*f(x)
///
/////////////////////////////////////////////////////////////////////////
GaussLaguerre_Coefficients::GaussLaguerre_Coefficients(long double alpha, int n)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	int i,its,j;
	double ai,p1,p2,p3,pp,z,z1;
	long_double r(0,1);
	for (i=0;i<n;i++) {
		if (i == 0) {
			z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
		} else if (i == 1) {
			z += (15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*n);
		} else {
			ai=i-1;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alpha/
				(1.0+3.5*ai))*(z-points[i-2])/(1.0+0.3*alpha);
		}
		for (its=0;its<ARM_CF_MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j+1+alpha-z)*p2-(j+alpha)*p3)/(j+1);
			}
			pp=(n*p1-(n+alpha)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= ARM_CF_EPS) break;
		}
		if (its >= ARM_CF_MAXIT)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"too many iterations in GaussLaguerre_Coefficients" );
		
		}
		points[i]=z;
		r=gammaLD(alpha+n)/gammaLD((double)n);
		weigths[i] = -r.todouble()/(pp*n*p2);
	}
}



/////////////////////////////////////////////////////////////////////////
///
/// Gauss legendre coefficients: can be used to integrals like
/// Integral{-1, +1} of f(x)
///
/////////////////////////////////////////////////////////////////////////


GaussLegendre_Coefficients::GaussLegendre_Coefficients(int n)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	GaussLegendre_Coefficients_repository* repo=GaussLegendre_Coefficients_repository::get_instance();
	Orthogonal_Polynomial_Coefficients* coefp=repo->get_coefficients(n);
	if(!coefp)
	{
		int m,j,i;
		long double z1,z,pp,p3,p2,p1;
		m=(n+1)/2;
		for (i=0;i<m;i++) 
		{
			z=cos(ARM_NumericConstants::ARM_PI*(i+0.75)/(n+0.5));
			do {
				p1=1.0;
				p2=0.0;
				for (j=0;j<n;j++) 
				{
					p3=p2;
					p2=p1;
					p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
				}
				pp=n*(z*p1-p2)/(z*z-1.0);
				z1=z;
				z=z1-p1/pp;
			} while (fabs(z-z1) > ARM_CF_EPS);
			points[i]=-z;
			points[n-1-i]=z;
			weigths[i]=2.0/((1.0-z*z)*pp*pp);
			weigths[n-1-i]=weigths[i];
		}
		repo->get_repository()->insert(make_pair(n,clone())); // we store it now in  the repository
	}
	else
	{
		for(int i=0;i<n;i++)		// the order n has already benn computed and is stored in the repository
		{
			points[i]=coefp->get_point(i);
			weigths[i]=coefp->get_weight(i);
		}
	}

}

GaussLegendre_Coefficients::GaussLegendre_Coefficients(GaussLegendre_Coefficients* coefs,double x1,double x2) 
:	Orthogonal_Polynomial_Coefficients(coefs->get_order())
{
	int n=coefs->get_order();
	int m,i;
	double xm,xl;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	m=(n+1)/2;
	for (i=0;i<m;i++) {	
		points[i]=xm+xl*coefs->get_point(i);
		points[n-1-i]=xm+xl*coefs->get_point(n-1-i);
		weigths[i]=coefs->get_weight(i)*xl;
		weigths[n-1-i]=coefs->get_weight(n-1-i)*xl;
	}
}

void GaussLegendre_Coefficients::AdaptCoefficients(GaussLegendre_Coefficients* coefs2,double x1,double x2) 
{
	int n=get_order();
	int m,i;
	double xm,xl;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	m=(n+1)/2;
	for (i=0;i<m;i++) {	
		coefs2->set_points(i,xm+xl*get_point(i));
		coefs2->set_points(n-1-i,xm+xl*get_point(n-1-i));
		coefs2->set_weigths(i,get_weight(i)*xl);
		coefs2->set_weigths(n-1-i,get_weight(n-1-i)*xl);
	}
}


/*
GaussLegendre_Coefficients::GaussLegendre_Coefficients(int n,double x1,double x2)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	GaussLegendre_Coefficients_repository* repo=GaussLegendre_Coefficients_repository::get_instance();
	Orthogonal_Polynomial_Coefficients* coefp=repo->get_coefficients(n);
	if(!coefp)
	{
		int m,j,i;
		double z1,z,xm,xl,pp,p3,p2,p1;
		xm=0.5*(x2+x1);
		xl=0.5*(x2-x1);
		m=(n+1)/2;
		for (i=0;i<m;i++) {
			z=cos(ARM_NumericConstants::ARM_PI*(i+0.75)/(n+0.5));
			do {
				p1=1.0;
				p2=0.0;
				for (j=0;j<n;j++) {
					p3=p2;
					p2=p1;
					p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
				}
				pp=n*(z*p1-p2)/(z*z-1.0);
				z1=z;
				z=z1-p1/pp;
			} while (fabs(z-z1) > ARM_CF_EPS);
			points[i]=xm-xl*z;
			points[n-1-i]=xm+xl*z;
			weigths[i]=2.0*xl/((1.0-z*z)*pp*pp);
			weigths[n-1-i]=weigths[i];
		}
		repo->get_repository()->insert(make_pair(n,clone())); // we store it now in  the repository
	}
	for(int i=0;i<n-1;i++)		// the order n has already benn computed and is stored in the repository
	{
		points[i]=coefp->get_point(i);
		weigths[i]=coefp->get_weight(i);
	}
}

*/
GaussLegendre_Coefficients::GaussLegendre_Coefficients(int n,double x1,double x2)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	GaussLegendre_Coefficients coefs(n);
	int m,i;
	double xm,xl;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	m=(n+1)/2;
	for (i=0;i<m;i++) {	
		points[i]=xm+xl*coefs.get_point(i);
		points[n-1-i]=xm+xl*coefs.get_point(n-1-i);
		weigths[i]=coefs.get_weight(i)*xl;
		weigths[n-1-i]=coefs.get_weight(n-1-i)*xl;
	}
}



/////////////////////////////////////////////////////////////////////////
///
/// ConstantSteps: can be used to integrales like
/// Integral{0, Xmax} of f(x)
///
/////////////////////////////////////////////////////////////////////////

ConstantSteps_Coefficients::ConstantSteps_Coefficients(double xmax, int n)
:	Orthogonal_Polynomial_Coefficients(n)  
{
	int i;
	for (i=0;i<n;i++) {
		points[i]=(i+0.5)*xmax/n;
		weigths[i]=xmax/n;
	}
}

/////////////////////////////////////////////////////////////////////////
///
/// Stratified (Hermite-Legendre) coefficients: can be used to integrales like
/// Integral{a, b} of exp(-x^2/2)/sqrt(2pi) *f(x)
///
/////////////////////////////////////////////////////////////////////////


GaussStratifiedHermiteLegendre_Coefficients::GaussStratifiedHermiteLegendre_Coefficients(GaussLegendre_Coefficients* coefs,int p,double a,double b)
:	Orthogonal_Polynomial_Coefficients((coefs->get_order())*p)  
{
	int n=coefs->get_order();
	int i,j;
	double a1,b1,c1;
	vector<double> ai(p+1);
	ai[0]=a;
	double na=NormalCDF(a);
	double delta=(NormalCDF(b)-na)/p;
	for(j=1;j<p;j++)
	{
		ai[j]=NormalCDFInverse(na+j*delta);
	}
	ai[p]=b;
	for(j=0;j<p;j++)
	{
		a1=ai[j];
		b1=ai[j+1];
		for(i=0;i<n;i++)
		{
			c1=((b1-a1)*coefs->get_point(i)+b1+a1)/2.;
			points[j*n+i]=c1;
			
			weigths[j*n+i]=coefs->get_weight(i)*exp(-c1*c1/2.)*(b1-a1)/(2.*ARM_NumericConstants::ARM_SQRT_2_PI);
		}
	}
	
}




/////////////////////////////////////////////////////////////////////////
///
/// Initialization of  the statics of the repository
/// 
///
/////////////////////////////////////////////////////////////////////////

GaussHermite_Coefficients_repository* GaussHermite_Coefficients_repository::instance=NULL;

GaussLegendre_Coefficients_repository* GaussLegendre_Coefficients_repository::instance=NULL;

double LegendreInt(const DoubleToDoubleFunc& func, double a, double b, int nbpts)
{
	GaussLegendre_Coefficients c(nbpts);

	double x, scale = b - a, sum = 0.;
	for(int i = 0; i < nbpts; i++)
	{
		x = a + scale * (0.5 + 0.5 * c.get_point(i));
		sum += func(x) * 0.5 * c.get_weight(i);
	}
	
	sum *= scale;

	return sum;
}

CC_END_NAMESPACE()


 

#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/