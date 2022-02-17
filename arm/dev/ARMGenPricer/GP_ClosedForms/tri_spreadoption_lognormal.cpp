/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file change_numeraire.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpnumlib/gaussiananalytics.h"

#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/spreadoption_lognormal.h"
#include "gpclosedforms/bi_spreadoption_lognormal.h"
#include "gpclosedforms/tri_spreadoption_lognormal.h"

#include "gpbase/numericconstant.h"

/// kernel
#include "glob/expt.h"

CC_BEGIN_NAMESPACE(ARM)


inline double cmin(double x, double y) {return ((x<=y) ? x : y);}
inline double cmax(double x, double y) {return ((x<=y) ? y : x);}





///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///			Computation of a gaussian integral equal to  : 
///			Integral[phi(x) phi(y) phi(z) 1{a0+a2*S2+a3*S3 >0} dxdydz
///	
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///   LogNormalIntegralType2: Auxiliary functions
///
//////////////////////////////////////////////////////////////////////////////////////////////



double LogNormalIntegralType2_Integrand_H(double a0,double d,double a1,double b1,double c1,
								   double a2,double b2,double c2,
								   double a3,double b3,double c3,double x,double y)
{
	return (1/d)*log( ((-a1)*exp(b1*x+c1*y)+(-a2)*exp(b2*x+c2*y)+(-a3)*exp(b3*x+c3*y))/a0 );

}

double LogNormalIntegralType2_Integrand_H2(double a0,double d,double a1,double b1,double c1,
								   double a2,double b2,double c2,
								   double a3,double b3,double c3,double x)
{
	return (log(-a1*exp(b1*x)/a3-a2*exp(b2*x)/a3)-b3*x)/(c3-c1);

}


/////////////////////////////////////////////////////////////////////////////////
///
///  LogNormalIntegralType2: first integration 
///
/////////////////////////////////////////////////////////////////////////////////
double LogNormalIntegralType2_Integrand001(double a0,double d,double a1,double b1,double c1,
								   double a2,double b2,double c2,
								   double a3,double b3,double c3,double x,GaussLegendre_Coefficients* glcoeffs_ptr)
{
	
	double resultat,epsd,epsc3c1,gaussmax,linf,lsup,effectivelinf,effectivelsup;
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	if(d>=0) {epsd=1.;} else {epsd=-1.;}
	if(c3-c1>=0) {epsc3c1=1.;} else {epsc3c1=-1.;}
	gaussmax=8.;
	/// cas n1
	if(((a1 > 0) && (a2 < 0) && (a3 < 0))|| ((a1 < 0) && (a2 > 0) && (a3 < 0)))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-gaussmax,gaussmax);
		
		resultat=0;
		for (i=0;i<n;i++) {
			resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType2_Integrand_H(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
		}
	}

	if(((a1 < 0) && (a2 > 0) && (a3 > 0) ) || ((a1 > 0) && (a2 < 0) && (a3 > 0) ))
	{	
		if((c3-c1)>0)
		{
				linf=-gaussmax;
				lsup=LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x);
		}
		else
		{
				linf=LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x);
				lsup=gaussmax;
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat=0;
		if(effectivelinf<effectivelsup)
		{	
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType2_Integrand_H(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
		}
		resultat+=ARM_GaussianAnalytics::cdfNormal(-epsc3c1*LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x));
	}
		/// clean up of the vectors
	if(range !=0) delete range;

		return  resultat;

}

double LogNormalIntegralType2_Integrand002(double a0,double d,double a1,double b1,double c1,
								   double a2,double b2,double c2,
								   double a3,double b3,double c3,double x,GaussLegendre_Coefficients* glcoeffs_ptr)
{
	
	double resultat,epsd,epsc3c1,gaussmax,linf,lsup,effectivelinf,effectivelsup;
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	if(d>=0) {epsd=1.;} else {epsd=-1.;}
	if(c3-c1>=0) {epsc3c1=1.;} else {epsc3c1=-1.;}
	gaussmax=8.;
	/// cas n1
	if(((a1>0)&&(a2<0)&&(a3<0) ) || ((a1<0)&&(a2>0)&&(a3<0) ) || ((a1>0)&&(a2>0)&&(a3<0)))
	{	
		if((c3-c1)>0)
		{
			linf=LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x);
			lsup=gaussmax;
		}
		else
		{
			linf=-gaussmax;
			lsup=LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x);
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat=0;
		if(effectivelinf<effectivelsup)
		{	
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType2_Integrand_H(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
		}
		resultat+=ARM_GaussianAnalytics::cdfNormal(epsc3c1*LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x));
	}
	if(((a1 < 0) && (a2 < 0) && (a3 > 0)))
	{	
		if((c3-c1)>0)
		{
			linf=-gaussmax;
			lsup=LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x);
		}
		else
		{
			linf=LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x);
			lsup=gaussmax;
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat=0;
		if(effectivelinf<effectivelsup)
		{	
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType2_Integrand_H(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
		}
		resultat+=ARM_GaussianAnalytics::cdfNormal(-epsc3c1*LogNormalIntegralType2_Integrand_H2(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x));
	}
	if(((a1 < 0) && (a2 < 0) && (a3 < 0) ))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-gaussmax,gaussmax);
		
		resultat=0;
		for (i=0;i<n;i++) {
			resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType2_Integrand_H(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
		}
	}

		if(range !=0) delete range;

		return  resultat;
}

/////////////////////////////////////////////////////////////////////////////////
///
///  LogNormalIntegralType2: second integration 
///
///  retourne la valeur de :
///  Int{phi(x)phi(y)phi(z)1{a0*exp(d*z)+a1*exp(b1*x+c1*y)+a2*exp(b2*x+c2*y)+a3*exp(b3*x+c3*y)>=0}}
///
/////////////////////////////////////////////////////////////////////////////////

double LogNormalIntegralType2_Compute(double a0,double d,double a1,double b1,double c1,
								   double a2,double b2,double c2,
								   double a3,double b3,double c3,
								   GaussLegendre_Coefficients* glcoeffs_ptr)
{
	double resultat1,epsb1b2,resultat2,gaussmax,linf,lsup,effectivelinf,effectivelsup;
	if(b1-b2>=0) {epsb1b2=1.;} else {epsb1b2=-1.;}
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	gaussmax=8.;
	
	if(((a1 > 0) && (a2 < 0) && (a3 < 0) )) {
		if((b1-b2)>0)
		{
			linf=-gaussmax;
			lsup=log(-a2/a1)/(b1-b2);
		}
		else
		{
			linf=log(-a2/a1)/(b1-b2);
			lsup=gaussmax;
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat1=0;
		if(effectivelinf<effectivelsup)
			resultat1=0;
		for (i=0;i<n;i++) {
			resultat1+=LogNormalIntegralType2_Integrand001(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
		if((b1-b2)>0)
		{
			linf=log(-a2/a1)/(b1-b2);
			lsup=gaussmax;
		}
		else
		{
			linf=-gaussmax;
			lsup=log(-a2/a1)/(b1-b2);
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		if(range !=0) delete range;
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat2=0;
		if(effectivelinf<effectivelsup)
			resultat2=0;
		for (i=0;i<n;i++) {
			resultat2+=LogNormalIntegralType2_Integrand002(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
	}
	if(((a1 < 0) && (a2 > 0) && (a3 < 0) )) {
		if((b1-b2)>0)
		{
			linf=log(-a2/a1)/(b1-b2);
			lsup=gaussmax;
		}
		else
		{
			linf=-gaussmax;
			lsup=log(-a2/a1)/(b1-b2);
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat1=0;
		if(effectivelinf<effectivelsup)
			resultat1=0;
		for (i=0;i<n;i++) {
			resultat1+=LogNormalIntegralType2_Integrand001(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
		if((b1-b2)>0)
		{
			linf=-gaussmax;
			lsup=log(-a2/a1)/(b1-b2);
		}
		else
		{
			linf=log(-a2/a1)/(b1-b2);
			lsup=gaussmax;
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		if(range !=0) delete range;
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat2=0;
		if(effectivelinf<effectivelsup)
			resultat2=0;
		for (i=0;i<n;i++) {
			resultat2+=LogNormalIntegralType2_Integrand002(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
	}
	if(((a1<0)&&(a2<0)&&(a3>0)) || ((a1>0)&&(a2>0)&&(a3<0)) || ((a1<0)&&(a2<0)&&(a3<0))) {
		resultat1=0;
		linf=-gaussmax;
		lsup=gaussmax;
		if(range !=0) delete range;
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,linf,lsup);
		resultat2=0;
		if(effectivelinf<effectivelsup)
			resultat2=0;
		for (i=0;i<n;i++) {
			resultat2+=LogNormalIntegralType2_Integrand002(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
	}
    if(((a1 < 0) && (a2 > 0) && (a3 > 0))) {
	
		resultat1=ARM_GaussianAnalytics::cdfNormal(epsb1b2*log(-a2/a1)/(b1-b2));
	
		if((b1-b2)>0)
		{
			linf=log(-a2/a1)/(b1-b2);
			lsup=gaussmax;
		}
		else
		{
			linf=-gaussmax;
			lsup=log(-a2/a1)/(b1-b2);
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		if(range !=0) delete range;
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat2=0;
		if(effectivelinf<effectivelsup)
			resultat2=0;
		for (i=0;i<n;i++) {
			resultat2+=LogNormalIntegralType2_Integrand001(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
	}
	if(((a1 > 0) && (a2 < 0) && (a3 > 0) )) {
	
		resultat1=ARM_GaussianAnalytics::cdfNormal(-epsb1b2*log(-a2/a1)/(b1-b2));
	
		if((b1-b2)>0)
		{
			linf=-gaussmax;
			lsup=log(-a2/a1)/(b1-b2);
		}
		else
		{
			linf=log(-a2/a1)/(b1-b2);
			lsup=gaussmax;
		}
		effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
		effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
		if(range !=0) delete range;
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
		resultat2=0;
		if(effectivelinf<effectivelsup)
			resultat2=0;
		for (i=0;i<n;i++) {
			resultat2+=LogNormalIntegralType2_Integrand001(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);	
		}
	}
	if(((a1 > 0) && (a2 > 0) && (a3 > 0) )) {
		resultat1=0;
		resultat2=1;
	}
	if(range !=0) delete range;
	return 	cmin(cmax(resultat1+resultat2,0.),1.);

}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///  double digital option where the condition are : 
///			a0 + a1*S1+a2*S2+a3*S3 >0
///		and b0+b1*S1>0
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///
///  Packaging of the digitale in fonction of the natural parameters 
///
/////////////////////////////////////////////////////////////////////////////////  


double TriSpreadDigitalCall_aux(double l1,double l2,double l3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double alpha0,double alpha1,double alpha2,double alpha3,
								   double T,GaussLegendre_Coefficients* glcoeffs_ptr)
{
	if(( 1.+2*r12*r13*r23-r12*r12-r13*r13-r23*r23)<=0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TriSpreadDigitalOption2 : bad correlations :");
	}
	
	double a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3;
	double rom=(r23-r12*r13)/sqrt(1-r12*r12);
	double sqrtT=sqrt(T);
	a0=alpha3*l3*exp((m3-v3*v3/2.)*T);
	d=v3*sqrt(1.-r13*r13-rom*rom)*sqrtT;
	a1=alpha0;
	b1=-r13*v3*sqrtT;
	c1=-v3*rom*sqrtT;
	a2=alpha1*l1*exp((m1-v1*v1/2.)*T);
	b2=(v1-v3*r13)*sqrtT;
	c2=-v3*rom*sqrtT;
	a3=alpha2*l2*exp((m2-v2*v2/2.)*T);
	b3=(v2*r12-v3*r13)*sqrtT;
	c3=(v2*sqrt(1-r12*r12)-v3*rom)*sqrtT;

	if(a0>0)
	{
	return LogNormalIntegralType2_Compute(a0,d,a1,b1,c1,a2,b2,c2,a3,b3,c3,glcoeffs_ptr);
	}
	else
	{
		return 1.-LogNormalIntegralType2_Compute(-a0,d,-a1,b1,c1,-a2,b2,c2,-a3,b3,c3,glcoeffs_ptr);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   TriSpreadDigitalCall :  
///	  which pays 1 if a0+a1*S1+a2*S2+a3*S3>0
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


double TriSpreadDigitalCall(double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a1,double a2,double a3,
								   double T,int n)
{
	if(1-r12*r12-r13*r13-r23*r23+2.*r12*r13*r23<0) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TriSpreadDigitalCall : bad correlations: 1-r1^2-r2^2-r3^2-2r1*r2*r3<0");
	}
		
	GaussLegendre_Coefficients glcoeffs(n);
	return TriSpreadDigitalCall_aux(S1,S2,S3,v1,v2,v3,m1,m2,m3,0.99999*r12,0.99999*r13,0.99999*r23,a0,a1,a2,a3,T,&glcoeffs);

}




//////////////////////////////////////////////////////////////////////////////////////////
///
///   Call/Put  a0<0, Nuance Introductors
///
//////////////////////////////////////////////////////////////////////////////////////////

double TriSpreadDigitalOption(double S1,double S2,double S3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a1,double a2,double a3,
								   double T,int callput,int n)
{
	switch (callput)
	{
	case K_CALL :
		{
			return TriSpreadDigitalCall(S1,S2,S3,v1,v2,v3,m1,m2,m3,r12,r13,r23,a0,a1,a2,a3,T,n);	
			break;
		}
	case K_PUT :
		{
			return 1.-TriSpreadDigitalCall(S1,S2,S3,v1,v2,v3,m1,m2,m3,r12,r13,r23,a0,a1,a2,a3,T,n);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TriSpreadDigitalOption : callput , bad input :");
			break;
		}
	}
	
}




///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///			Computation of a gaussian integral equal to  : 
///			Integral[phi(x) phi(y) phi(z) 1{a0 + a1*S1+a2*S2+a3*S3 >0}1{b0+b1*S1>0} dxdydz
///	
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
///
///  LogNormalIntegralType3 : Auxiliary functions
///
//////////////////////////////////////////////////////////////////////////////////////////////



double LogNormalIntegralType3_Integrand_H(double d,double a1,double a2,double b2,double c2,
								   double a3,double b3,double c3,double x,double y)
{
	return (1/d)*log( ((-a2)*exp(b2*x+c2*y)+(-a3)*exp(b3*x+c3*y))/a1 );

}

double LogNormalIntegralType3_Integrand_H2(double d,double a1,double a2,double b2,double c2,
								   double a3,double b3,double c3,double x)
{
	return (-(b2-b3)*x+log(-a3/a2))/(c2-c3);

}

double LogNormalIntegralType3_Integrand_Q(double k,double l,double m)
{
	return log(-k/l)/m;
}

/////////////////////////////////////////////////////////////////////////////////
///
///  LogNormalIntegralType3: first integration 
///
/////////////////////////////////////////////////////////////////////////////////



double LogNormalIntegralType3_Integrand(double d,double a1,double a2,double b2,double c2,
								   double a3,double b3,double c3,
								   double k,double l,double m,
								   double x,GaussLegendre_Coefficients* glcoeffs_ptr)
{
	double resultat,h2,epsd,limsup;
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	if(d>=0) {epsd=1.;} else {epsd=-1.;}
	limsup=8.;
	if((a2!=0)&&(a3!=0)) {h2=LogNormalIntegralType3_Integrand_H2(d,a1,a2,b2,c2,a3,b3,c3,x);}
	/// cas n1
	if((a1 > 0) && (((a2 < 0) && (a3 > 0) && (c2 > c3)) || ((a2 > 0) && (a3 < 0) && (c2 < c3))))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,cmin(limsup,cmax(h2,-limsup)),limsup);
		if(h2>=limsup) resultat=0;
		else
		{
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType3_Integrand_H(d,a1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
			resultat+=ARM_GaussianAnalytics::cdfNormal(h2);
		}

	}
	/// cas n2
	if((a1 < 0) && (((a2 < 0) && (a3 > 0) && (c2 > c3)) || ((a2 > 0) && (a3 < 0) && (c2 < c3))))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-limsup, cmax(-limsup,cmin(h2,limsup)));
		if(h2<=-limsup) resultat=0;
		else
		{
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*LogNormalIntegralType3_Integrand_H(d,a1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
		}

	}
	/// cas n3
	if((a1 > 0) && (((a2 > 0) && (a3 < 0) && (c2 > c3)) || ((a2 < 0) && (a3 > 0) && (c2 < c3))))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-limsup,cmax(-limsup,cmin(h2,limsup)));
		if(h2<=-limsup) resultat=0;
		else
		{
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType3_Integrand_H(d,a1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
			resultat+=ARM_GaussianAnalytics::cdfNormal(-h2);
		}

	}
	/// cas n4
	if((a1 < 0) && (((a2 > 0) && (a3 < 0) && (c2 > c3)) || ((a2 < 0) && (a3 > 0) && (c2 < c3))))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,cmin(limsup,cmax(h2,-limsup)),limsup);
		if(h2>=limsup) resultat=0;
		else
		{
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*LogNormalIntegralType3_Integrand_H(d,a1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
			}
		}

	}
	/// cas n5
	if(((a1 > 0) && (a2 <= 0) && (a3 < 0)) || ((a1 > 0) && (a2 < 0) && (a3 <= 0)))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-limsup,limsup);
		resultat=0;
		for (i=0;i<n;i++) {
			resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*LogNormalIntegralType3_Integrand_H(d,a1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
		}

	}
	/// cas n6
	if(((a1 < 0) && (a2 >= 0) && (a3 > 0)) || ((a1 < 0) && (a2 > 0) && (a3 >= 0)))
	{	
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-limsup,limsup);
		resultat=0;
		for (i=0;i<n;i++) {
			resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*LogNormalIntegralType3_Integrand_H(d,a1,a2,b2,c2,a3,b3,c3,x,range->get_point(i)))*range->get_weight(i);
		}

	}
	/// cas n7
	if((a1 >= 0) && (a2 >= 0) && (a3 >= 0)) resultat=1;
	/// cas n8
	if((a1 < 0) && (a2 < 0) && (a3 < 0)) resultat=0;


	/// clean up of the vectors
	if(range !=0) delete range;

return resultat;
}


/////////////////////////////////////////////////////////////////////////////////
///
///  LogNormalIntegralType3: second integration 
///  retourne la valeur de :
///  Int{phi(x)phi(y)phi(z)1{a1*exp(d*z)+a2*exp(b2*x+c2*y)+a3*exp(b3*x+c3*y)>=0}1{k+l*exp(m*x)>=0}}
///
/////////////////////////////////////////////////////////////////////////////////


double LogNormalIntegralType3_Compute(double d,double a1,double a2,double b2,double c2,
								   double a3,double b3,double c3,
								   double k,double l,double m,
								   GaussLegendre_Coefficients* glcoeffs_ptr)
{
	double resultat,limsup,q;
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	limsup=8.;
	if (((l > 0) && (k < 0)) || ((l < 0) && (k > 0))) { q=LogNormalIntegralType3_Integrand_Q(k,l, m);}
	if(((m > 0) && (l > 0) && (k < 0)) || ((m < 0) && (l < 0) && (k > 0))) {
		range= new ReducedGaussHermite_Coefficients(glcoeffs_ptr,cmin(limsup,cmax(q,-limsup)),limsup);
	}
	if(((m < 0) && (l > 0) && (k < 0)) || ((m > 0) && (l < 0) && (k > 0))) {
		range= new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-limsup,cmax(-limsup,cmin(q,limsup)));
	}
	if((l >= 0) && (k >= 0)) {
		range= new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-limsup,limsup);
	}
	resultat=0;
	for (i=0;i<n;i++) {
		resultat+=LogNormalIntegralType3_Integrand(d,a1,a2,b2,c2,a3,b3,c3,k,l,m,range->get_point(i),glcoeffs_ptr)*range->get_weight(i);
	}

	if(range !=0) delete range;
	return cmin(cmax(resultat,0.),1.);

}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///  double digital option where the condition are : 
///			a0 + a1*S1+a2*S2+a3*S3 >0
///		and b0+b1*S1>0
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///
///  Packaging of the digitale in fonction of the natural parameters 
///
/////////////////////////////////////////////////////////////////////////////////

double TriSpreadDigitalCall2_aux(double l1,double l2,double l3,
							  double m1,double m2,double m3,
							  double v1,double v2,double v3,
							  double r12,double r13,double r23_0,
							  double alpha0,double alpha2,double alpha3,
							  double beta0,double beta1,double T,GaussLegendre_Coefficients* glcoeffs_ptr)
{
	if(r12>0.99999)	/// use of the 2 dim formula for degenerated cases rho12=1
	{
		return Generalized_BiSpreadDigitalCall_aux(
			l1, l3 ,m1, m3 ,v1, v3, r13, 
			alpha0, alpha2*exp(-T*(m1-v1*v1/2.)*v2/v1+T*(m2-v2*v2/2.))*pow(l1,-v2/v1)*l2,
			alpha3,beta0, beta1, T, 
			v2/v1, 1., 1.,glcoeffs_ptr);
	}
	if(r12<-0.99999)	/// use of the 2 dim formula for degenerated cases rho12=-1
	{
		return Generalized_BiSpreadDigitalCall_aux(
			l1, l3 ,m1, m3 ,v1, v3, r13, 
			alpha0, alpha2*exp(T*(m1-v1*v1/2.)*v2/v1+T*(m2-v2*v2/2.))*pow(l1,v2/v1)*l2,
			alpha3,beta0, beta1, T, 
			-v2/v1, 1., 1.,glcoeffs_ptr);
	}
	if(r13>0.99999)	/// use of the 2 dim formula for degenerated cases rho13=1
	{
		return Generalized_BiSpreadDigitalCall_aux(
			l1, l2 ,m1, m2 ,v1, v2, r12, 
			alpha0, alpha3*exp(-T*(m1-v1*v1/2.)*v3/v1+T*(m3-v3*v3/2.))*pow(l1,-v3/v1)*l3,
			alpha2,beta0, beta1, T, 
			v3/v1, 1., 1.,glcoeffs_ptr);
	}
	if(r13<-0.99999)	/// use of the 2 dim formula for degenerated cases rho13=-1
	{
		return Generalized_BiSpreadDigitalCall_aux(
			l1, l2 ,m1, m2 ,v1, v2, r12, 
			alpha0, alpha3*exp(T*(m1-v1*v1/2.)*v3/v1+T*(m3-v3*v3/2.))*pow(l1,v3/v1)*l3,
			alpha2,beta0, beta1, T, 
			-v3/v1, 1., 1.,glcoeffs_ptr);
	}
	double r23;		// on limite ici la valeur ce qui ne doit pas influer
	if(r23_0>0.99999) 
	{
		r23=0.99999; 
	}
	else 
	{
		if(r23_0<-0.99999) 
		{
			r23=-0.99999;
		}
		else 
		{
			r23=r23_0;
		}
	}

	double d,a1,a2,b2,c2,a3,b3,c3,k,l,m;
	a1=alpha3*l3*exp((m3-v3*v3/2.)*T);
	d=v3*sqrt((1.-r13*r13-(r23-r12*r13)*(r23-r12*r13)/(1-r12*r12))*T);
	a2=alpha0;
	b2=-r13*v3*sqrt(T);
	c2=-v3*(r23-r12*r13)*sqrt(T/(1-r12*r12));
	a3=alpha2*l2*exp((m2-v2*v2/2.)*T);
	b3=(r12*v2-r13*v3)*sqrt(T);
	c3=(v2*sqrt(1-r12*r12)-v3*(r23-r12*r13)/sqrt(1-r12*r12))*sqrt(T);
	k=beta0;
	l=beta1*l1*exp((m1-v1*v1/2.)*T);
	m=v1*sqrt(T);


	return LogNormalIntegralType3_Compute(d,a1,a2,b2,c2,a3,b3,c3,k,l,m,glcoeffs_ptr);

}


/////////////////////////////////////////////////////////////////////////////////
///
///   final Packaging of the digitale with computation of the Legendre points
///
/////////////////////////////////////////////////////////////////////////////////


double TriSpreadDigitalCall2(	   double l1,double l2,double l3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a2,double a3,
								   double b0,double b1,
								   double T,int n)
{
	if(1-r12*r12-r13*r13-r23*r23+2.*r12*r13*r23<0) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TriSpreadDigitalCall2 : bad correlations: 1-r1^2-r2^2-r3^2-2r1*r2*r3<0");
	}
/*	if ((r12==1.)||(r12==-1.)||(r13==1.)||(r13==-1.)||(r23==1.)||(r23==-1.))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TriSpreadDigitalCall2 : bad correlation =+or-1: you should use a bidim model");
	}
	*/

	GaussLegendre_Coefficients glcoeffs(n);
	return TriSpreadDigitalCall2_aux(l1,l2,l3,m1,m2,m3,v1,v2,v3,r12,r13,r23,a0,a2,a3,b0,b1,T,&glcoeffs);
}


//////////////////////////////////////////////////////////////////////////////////////////
///
///   Call/Put   Nuance Introductors
///
//////////////////////////////////////////////////////////////////////////////////////////

double TriSpreadDigitalOption2(double l1,double l2,double l3,
								   double v1,double v2,double v3,
								   double m1,double m2,double m3,
								   double r12,double r13,double r23,
								   double a0,double a2,double a3,
								   double b0,double b1,
								   double T,int callput,int n)
{
	switch (callput)
	{
	case K_CALL :
		{
			return TriSpreadDigitalCall2(l1,l2,l3,v1,v2,v3,m1,m2,m3,r12,r13,r23,a0,a2,a3,b0,b1,T,n);
			break;
		}
	case K_PUT :
		{
			return 1.-TriSpreadDigitalCall2(l1,l2,l3,v1,v2,v3,m1,m2,m3,r12,r13,r23,a0,a2,a3,b0,b1,T,n);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TriSpreadDigitalOption2 : callput , bad input :");
			break;
		}
	}
	
}


CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
