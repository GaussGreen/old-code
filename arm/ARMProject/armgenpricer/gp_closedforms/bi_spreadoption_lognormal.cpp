/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bi_spreadoption_lognormal.cpp
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

#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/spreadoption_lognormal.h"
#include "gpclosedforms/bi_spreadoption_lognormal.h"

#include "gpbase/numericconstant.h"

#include "expt.h"


CC_BEGIN_NAMESPACE(ARM)

inline double cmin(double x, double y) {return ((x<=y) ? x : y);}
inline double cmax(double x, double y) {return ((x<=y) ? y : x);}
inline double cnormalize(double x,double limsup) {return ((x<-limsup) ? -limsup : ((x> limsup) ? limsup : x));}


///////////////////////////////////////////////////////////////////////////////////////////////
///
///  double digital option where the condition are : 
///			a0 + a1*S1+a2*S2 >0
///		and b0+b1*S1>0
///
//////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
///
///  Auxiliary functions
///
//////////////////////////////////////////////////////////////////////////////////////////////



double biSpreadDigitalIntegrand_H(double d,double a1,double a2,double b2,
								   double a3,double b3,double x)
{
	return (1/d)*log( ((-a2)*exp(b2*x)+(-a3)*exp(b3*x))/a1 );

}

double biSpreadDigitalIntegrand_H2(double d,double a1,double a2,double b2,
								   double a3,double b3)
{
	return log(-a3/a2)/(b2-b3);

}

double biSpreadDigitalIntegrand_Q(double k,double l,double m)
{
	return log(-k/l)/m;
}

/////////////////////////////////////////////////////////////////////////////////
///
///  computation of the monodimensional integral
/// Int{phi(x)phi(y) 1{a1*exp(d*y)+a2*exp(b2*x)+a3*exp(b3*x)>=0}1{k+l*exp(m*x)>=0}}
///
/////////////////////////////////////////////////////////////////////////////////



double BilognormalIntegralType1_compute(double d,double a1,double a2,double b2,
										 double a3,double b3,
										 double k,double l,double m,
										GaussLegendre_Coefficients* glcoeffs_ptr)
{
	
	double resultat,epsd,h2,q,gaussmax,linf,lsup,effectivelinf,effectivelsup;
	if(d>=0) {epsd=1.;} else {epsd=-1.;}
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	gaussmax=8.;
	if(k!=0)
	{
		q=biSpreadDigitalIntegrand_Q(k,l,m);
		if((a1>0)&&(((a2<0)&&(a3>0)&&(b2>b3))||((a2>0)&&(a3<0)&&(b2<b3)))) {		///  cas 1
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			if(m>0)
			{
				linf=cmax(h2,q);lsup=gaussmax;
			}
			else
			{
				if (q>h2) 
				{
					linf=h2;lsup=q;
				}
				else
				{
					linf=gaussmax;lsup=gaussmax;
				}
			}
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
			if((m>0)&&(q<h2))
			{
				resultat+=ARM_GaussianAnalytics::cdfNormal(h2)-ARM_GaussianAnalytics::cdfNormal(q);
			}
			if(m<0)
			{
				resultat+=ARM_GaussianAnalytics::cdfNormal(cmin(q,h2));
			}
			
		}
		if((a1<0)&&(((a2<0)&&(a3>0)&&(b2>b3)) || ((a2>0)&&(a3<0)&&(b2<b3)))) {		///  cas 2
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			if(m>0)
			{
				if (q<h2) 
				{
					linf=q;lsup=h2;
				}
				else
				{
					linf=gaussmax;lsup=gaussmax;
				}
			}
			else
			{
				linf=-gaussmax;lsup=cmin(h2,q);
			}
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
		}
		if((a1>0)&&(((a2>0)&&(a3<0)&&(b2>b3)) || ((a2<0)&&(a3>0)&&(b2<b3)))) {		///  cas 3
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			if(m>0)
			{
				if (q<h2) 
				{
					linf=q;lsup=h2;
				}
				else
				{
					linf=gaussmax;lsup=gaussmax;
				}
			}
			else
			{
				linf=-gaussmax;lsup=cmin(h2,q);
			}
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
			if(m>0)
			{
				resultat+=ARM_GaussianAnalytics::cdfNormal(-cmax(q,h2));
			}
			if((m<0)&&(q>h2))
			{
				resultat+=ARM_GaussianAnalytics::cdfNormal(q)-ARM_GaussianAnalytics::cdfNormal(h2);
			}
			
			
		}
		if((a1<0)&&(((a2>0)&&(a3<0)&&(b2>b3)) || ((a2<0)&&(a3>0)&&(b2<b3)))) {		///  cas 4
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			if(m>0)
			{
				linf=cmax(h2,q);lsup=gaussmax;
			}
			else
			{
				if (q>h2) 
				{
					linf=h2;lsup=q;
				}
				else
				{
					linf=gaussmax;lsup=gaussmax;
				}
			}
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
		}
		if((((a1>0)&&(a2<=0)&&(a3<0)) || ((a1>0)&&(a2<0)&&(a3<=0)))) {				///  cas 5
			if(m>0)
			{
				linf=q;lsup=gaussmax;
			}
			else
			{
				linf=-gaussmax;lsup=q;
			}
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}
		}
		if((((a1<0)&&(a2>=0)&&(a3>0)) || ((a1<0)&&(a2>0)&&(a3>=0)))) {				///  cas 6
			if(m>0)
			{
				linf=q;lsup=gaussmax;
			}
			else
			{
				linf=-gaussmax;lsup=q;
			}
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}
		}
		if(((a1>=0) && (a2>=0) && (a3>=0))) {		///  cas 7
			resultat=ARM_GaussianAnalytics::cdfNormal(-q);
		}
		if(((a1<0) && (a2<0) && (a3<0))) {		///  cas 8
			resultat=0;
		}
	}
	if(k==0)
	{
		if(l<0) return 0;
		if((a1>0)&&(((a2<0)&&(a3>0)&&(b2>b3)) || ((a2>0)&&(a3<0)&&(b2<b3)))) {		///  cas 9
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			linf=h2;lsup=gaussmax;
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
			
		}
		if((a1<0)&&(((a2<0)&&(a3>0)&&(b2>b3)) || ((a2>0)&&(a3<0)&&(b2<b3)))) {		///  cas 10
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			linf=-gaussmax;lsup=h2;
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
			
		}
		if((a1>0)&&(((a2>0)&&(a3<0)&&(b2>b3)) || ((a2<0)&&(a3>0)&&(b2<b3)))) {		///  cas 11
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			linf=-gaussmax;lsup=h2;
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
			
		}
		if((a1<0)&&(((a2>0)&&(a3<0)&&(b2>b3)) || ((a2<0)&&(a3>0)&&(b2<b3)))) {		///  cas 12
			h2=biSpreadDigitalIntegrand_H2(d,a1,a2,b2,a3,b3);
			linf=h2;lsup=gaussmax;
			effectivelinf=cmin(gaussmax,cmax(linf,-gaussmax));
			effectivelsup=cmax(-gaussmax,cmin(lsup,gaussmax));
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,effectivelinf,effectivelsup);
			resultat=0;
			if(effectivelinf<effectivelsup){
				for (i=0;i<n;i++) {
					resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}}
			
		}
		if((((a1>0)&&(a2<=0)&&(a3<0)) || ((a1>0)&&(a2<0)&&(a3<=0)))) {				///  cas 13
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-gaussmax,gaussmax);
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(-epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}
		}
		if((((a1<0)&&(a2>=0)&&(a3>0)) || ((a1<0)&&(a2>0)&&(a3>=0)))) {				///  cas 14
			range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-gaussmax,gaussmax);
			resultat=0;
			for (i=0;i<n;i++) {
				resultat+=ARM_GaussianAnalytics::cdfNormal(epsd*biSpreadDigitalIntegrand_H(d,a1,a2,b2,a3,b3,range->get_point(i)))*range->get_weight(i);
			}
		}
	}
	return resultat;
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///  double digital option where the condition are : 
///			a0 + a1*S1+a2*S2 >0
///		and b0+b1*S1>0
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///
///              implementation of the generalized digital: 
///                  payoff =1	only if 	a0 + a1*S1^g1+a2*S2^g2 >0
///											and b0+b1*S1^gb>0
///
/////////////////////////////////////////////////////////////////////////////////

double Generalized_BiSpreadDigitalCall_aux(double l1,double l2,
							  double m1,double m2,
							  double sig1,double sig2,
							  double r12,
							  double alpha0,double alpha1,double alpha2,
							  double beta0,double beta1,double T,
							  double g1,double g2,double gb,
							  GaussLegendre_Coefficients* glcoeffs_ptr)
{
	double d,a1,a2,b2,a3,b3,k,l,m;
	a1=alpha2*pow(l2,g2)*exp(g2*(m2-sig2*sig2/2.)*T);
	d=g2*sig2*sqrt((1.-r12*r12)*T);
	a2=alpha0;
	b2=-g2*r12*sig2*sqrt(T);
	a3=alpha1*pow(l1,g1)*exp(g1*(m1-sig1*sig1/2.)*T);
	b3=(g1*sig1-g2*r12*sig2)*sqrt(T);
	k=beta0;
	l=beta1*pow(l1,gb)*exp(gb*(m1-sig1*sig1/2.)*T);
	m=sig1*sqrt(T);
	if(l>0)
	{
		return BilognormalIntegralType1_compute(d,a1,a2,b2,a3,b3,k,l,m,glcoeffs_ptr);
	}
	else
	{
		double inverseoption=BilognormalIntegralType1_compute(d,a1,a2,b2,a3,b3,-k,-l,m,glcoeffs_ptr);
		
	a1=alpha2*pow(l2,g2)*exp(g2*(m2-sig2*sig2/2.)*T);
	d=g2*sig2*sqrt((1.-r12*r12)*T);
	a2=alpha0;
	b2=-g2*r12*sig2*sqrt(T);
	a3=alpha1*pow(l1,g1)*exp(g1*(m1-sig1*sig1/2.)*T);
	b3=(g1*sig1-g2*r12*sig2)*sqrt(T);
	k=1;
	l=1;
	m=sig1*sqrt(T);
		
		return BilognormalIntegralType1_compute(d,a1,a2,b2,a3,b3,k,l,m,glcoeffs_ptr)-inverseoption;
	}
}

/////////////////////////////////////////////////////////////////////////////////
///
/// packaging : association of the digitale with computation of the Legendre points
///
/////////////////////////////////////////////////////////////////////////////////


double BiSpreadDigitalCall(	   double l1,double l2,
								   double v1,double v2,
								   double m1,double m2,
								   double r12,
								   double a0,double a1,double a2,
								   double b0,double b1,
								   double T,int n)
{

	GaussLegendre_Coefficients glcoeffs(n);
	return Generalized_BiSpreadDigitalCall_aux(l1,l2,m1,m2,v1,v2,r12,a0,a1,a2,b0,b1,T,1.,1.,1.,&glcoeffs);
}


//////////////////////////////////////////////////////////////////////////////////////////
///
///  packaging: Call/Put   Nuance Introductors
///
//////////////////////////////////////////////////////////////////////////////////////////

double BiSpreadDigitalOption(double l1,double l2,
								   double v1,double v2,
								   double m1,double m2,
								   double r12,
								   double a0,double a1,double a2,
								   double b0,double b1,
								   double T,int callput,int n)
{
	switch (callput)
	{
	case K_CALL :
		{
			return BiSpreadDigitalCall(l1,l2,v1,v2,m1,m2,r12,a0,a1,a2,b0,b1,T,n);
			break;
		}
	case K_PUT :
		{
			return 1.-BiSpreadDigitalCall(l1,l2,v1,v2,m1,m2,r12,a0,a1,a2,b0,b1,T,n);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSpreadDigitalOption : callput , bad input :");
			break;
		}
	}
	
}

CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
