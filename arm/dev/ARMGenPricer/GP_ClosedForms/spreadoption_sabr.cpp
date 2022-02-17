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


//////////////////////////////////////////////////////////////////////////////////////////
///
///   Regular Spreadoption Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////

/// we use generic pricing mechanism based on structure::PowerSpreadOption_Pricing 
/// no specific functions here

//////////////////////////////////////////////////////////////////////////////////////////
///
///    Spreadoption paying S3 Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////

#include "gpclosedforms/sabrvanilla.h"

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/numerics_Interface.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/spreadoption_sabr.h"
#include "gpclosedforms/distribution_interface.h"

#include <glob/expt.h>


CC_BEGIN_NAMESPACE(ARM)
 

double GaussianSABRDigitalCall(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	double nbstddev=3.;
	int nbsteps=120;int i;
	GaussLegendre_Coefficients c_root(legendreNb);
	double	Sum=0.;double x,y1,y2,y3,y4;
	GaussLegendre_Coefficients c1(&c_root,-nbstddev*f1,nbstddev*f1);
				for(i=0;i<legendreNb;i++){
					x=c1.get_point(i);
						y1=Export_SABR_Distribution(f1,x+K,T,alpha1,beta1,rho1,nu1,SABRFlag1,nbsteps, alpha_exp, alpha_tanh, kb_tanh);
						y2=Export_SABR_Distribution(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh);
						y3=Export_SABR_Density(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh);
						y4=DyBivariateCopula(y1,y2,rho);
					Sum+= (1.-DyBivariateCopula(
						Export_SABR_Distribution(f1,x+K,T,alpha1,beta1,rho1,nu1,SABRFlag1,nbsteps, alpha_exp, alpha_tanh, kb_tanh),
						Export_SABR_Distribution(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh),
						rho))*
						Export_SABR_Density(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh)*
						c1.get_weight(i);
				}
	return Sum;

}

/*
double GaussianSABRDigitalCallPayingS1(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	double nbstddev=3.;
	int nbsteps=120;int i;
	GaussLegendre_Coefficients c_root(legendreNb);
	double	Sum=0.;double x;
	GaussLegendre_Coefficients c1(&c_root,-nbstddev*f1,nbstddev*f1);
				for(i=0;i<legendreNb;i++){
					x=c1.get_point(i);
					Sum+= ((x+K)*DxBivariateCopula(
						Export_SABR_Distribution(f1,x+K,T,alpha1,beta1,rho1,nu1,SABRFlag1,nbsteps, alpha_exp, alpha_tanh, kb_tanh),
						Export_SABR_Distribution(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh),
						rho))*
						Export_SABR_Density(f1,x+K,T,alpha1,beta1,rho1,nu1,SABRFlag1,nbsteps, alpha_exp, alpha_tanh, kb_tanh)*
						c1.get_weight(i);
				}
	return Sum;

}
*/

double GaussianSABRDigitalCallPayingS2(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	double nbstddev=3.;
	int nbsteps=120;int i;
	GaussLegendre_Coefficients c_root(legendreNb);
	double	Sum=0.;double x;
	GaussLegendre_Coefficients c1(&c_root,-nbstddev*f1,nbstddev*f1);
				for(i=0;i<legendreNb;i++){
					x=c1.get_point(i);
					Sum+= x*(1.-DyBivariateCopula(
						Export_SABR_Distribution(f1,x+K,T,alpha1,beta1,rho1,nu1,SABRFlag1,nbsteps, alpha_exp, alpha_tanh, kb_tanh),
						Export_SABR_Distribution(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh),
						rho))*
						Export_SABR_Density(f2,x,T,alpha2,beta2,rho2,nu2,SABRFlag2,nbsteps, alpha_exp, alpha_tanh, kb_tanh)*
						c1.get_weight(i);
				}
	return Sum;

}

double GaussianSABRDigitalCallPayingS1(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	return f1-GaussianSABRDigitalCallPayingS2( f2,
							    alpha2,
							    beta2,
							    rho2,
							    nu2,
							    SABRFlag2,
							    f1,
							    alpha1,
							    beta1,
							    rho1,
							    nu1,
							    SABRFlag1,
							    rho,
							    -K,
							    T,
							    legendreNb, alpha_exp, alpha_tanh, kb_tanh
							   );
}

/*
double b= rho23*sigma3*f3*f1beta1*alpha1-rho13*sigma3*f3*rho12*f2beta2*alpha2;
	b=b/(f2beta2*f2beta2*alpha2*alpha2-rho12*f2beta2*alpha2/(alpha1*f1beta1));
	double a=(rho13*sigma3*f3-b*rho12*f2beta2*alpha2)/(alpha1*f1beta1);
	*/
double GaussianSABRDigitalCallPayingS3(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double f3,
							   double sigma3,
							   double rho12,
							   double rho13,
							   double rho23,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	double X1=pow(f1,beta1)*alpha1;
	double X2=pow(f2,beta2)*alpha2;
	double X3=f3*sigma3;
	double a=X3/X1*(rho12*rho23-rho13)/(rho12*rho12-1.);
	double b=X2/X1*(rho12*rho13-rho13)/(rho12*rho12-1.);
	double D0=GaussianSABRDigitalCall(f1,alpha1,beta1,rho1,nu1,SABRFlag1,
							    f2,alpha2,beta2,rho2,nu2,SABRFlag2,
							    rho12,K,T,legendreNb, alpha_exp, alpha_tanh, kb_tanh);
	double D1=GaussianSABRDigitalCallPayingS1(f1,alpha1,beta1,rho1,nu1,SABRFlag1,
							    f2,alpha2,beta2,rho2,nu2,SABRFlag2,
							    rho12,K,T,legendreNb, alpha_exp, alpha_tanh, kb_tanh);
	double D2=GaussianSABRDigitalCallPayingS2(f1,alpha1,beta1,rho1,nu1,SABRFlag1,
							    f2,alpha2,beta2,rho2,nu2,SABRFlag2,
							    rho12,K,T,legendreNb, alpha_exp, alpha_tanh, kb_tanh);

	return f3*D0+a*(D1-f1*D0)+b*(D2-f2*D0);
}

CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/