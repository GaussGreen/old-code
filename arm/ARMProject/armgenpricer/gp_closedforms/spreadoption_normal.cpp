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

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/spreadoption_normal_formula.h"

#include "expt.h"



CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////

double NormalExpIntegral(double a,double b,GaussLegendre_Coefficients* glcoeffs_ptr)
{
	double resultat,y,gaussmax,sum=0.;
	int i,n=glcoeffs_ptr->get_order();
	ReducedGaussHermite_Coefficients* range=0;
	gaussmax=8.;
		range=new ReducedGaussHermite_Coefficients(glcoeffs_ptr,-gaussmax,gaussmax);
		
		resultat=0;
		for (i=0;i<n;i++) {
			y=range->get_point(i);
			resultat+=exp(-y*y/2.)*ARM_GaussianAnalytics::cdfNormal(a+b*y)*range->get_weight(i);
		}
		return sum/ARM_NumericConstants::ARM_SQRT_2_PI;
}


double NormalSpreadOptionIntegral_eta(double a,double b,double c)
{
	return ARM_GaussianAnalytics::cdfNormal(a/sqrt(b*b+c*c));
}

double NormalSpreadOptionIntegral_phi(double a,double b,double c)
{
	double signb;
	if (b>=0) signb=1.; else signb= -1.;
	return signb*exp(-a*a/(2.*(b*b+c*c)))/sqrt(1.+c*c/(b*b))*ARM_NumericConstants::ARM_INVSQRT2PI;

}


double SpreadDigitalOption_N(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callOrPut,int optiontype)
{
	int nb_legendre_points=30;
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_N_Formula::DIGITALOPTION :
		{
			double standardDev=sqrt(sig1*sig1+sig2*sig2-2.*sig1*sig2*rho);
			return VanillaDigitalOption_N(S1-S2,standardDev,k,t,callOrPut);
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::SPREADOPTION :
		{
			double standardDev=sqrt(sig1*sig1+sig2*sig2-2.*sig1*sig2*rho);
			return VanillaOption_N(S1-S2,standardDev,k,t,callOrPut);
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::PAYFIRST :
		{	
			
			double a=S1-S2-k;
			double b=sig1-sig2*rho;
			double c=-sig2*sqrt(1.-rho*rho);
			return S1*NormalSpreadOptionIntegral_eta(a,b,c)+sig1*NormalSpreadOptionIntegral_phi(a,b,c);
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::PAYSECOND :
		{
			
			double a=S1-S2-k;
			double b=sig1-sig2*rho;
			double c=-sig2*sqrt(1.-rho*rho);
			return S2*NormalSpreadOptionIntegral_eta(a,b,c)+
				   sig2*rho*NormalSpreadOptionIntegral_phi(a,b,c)+
				   sig2*sqrt(1.-rho*rho)*NormalSpreadOptionIntegral_phi(a,c,b);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Normal_SpreadOption : incorrect optiontype toggle ");
		}
	}

	
}

double Vega1SpreadDigitalOption_N(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callOrPut,int optiontype)
{
	int nb_legendre_points=30;
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_N_Formula::DIGITALOPTION :
		{
			double standardDev=sqrt(sig1*sig1+sig2*sig2-2.*sig1*sig2*rho);
			double standardDev_Vega1=(sig1-rho*sig2)/standardDev;
			return VegaVanillaDigitalOption_N(S1-S2,standardDev,k,t,callOrPut)*standardDev_Vega1;
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::SPREADOPTION :
		{
			double standardDev=sqrt(sig1*sig1+sig2*sig2-2.*sig1*sig2*rho);
			double standardDev_Vega1=(sig1-rho*sig2)/standardDev;
			return VegaVanillaOption_N(S1-S2,standardDev,k,t,callOrPut)*standardDev_Vega1;
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::PAYFIRST :
		{	
			
			double shift=0.001;
			double a=S1-S2-k;
			double b=sig1-sig2*rho;
			double c=-sig2*sqrt(1.-rho*rho);
			double b_shifted_v1=b+shift;
			double v0=S1*NormalSpreadOptionIntegral_eta(a,b,c)+sig1*NormalSpreadOptionIntegral_phi(a,b,c);
			double v1=S1*NormalSpreadOptionIntegral_eta(a,b_shifted_v1,c)+sig1*NormalSpreadOptionIntegral_phi(a,b_shifted_v1,c);
			return (v1-v0)/shift;
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::PAYSECOND :
		{
		
			double shift=0.001;
			double a=S1-S2-k;
			double b=sig1-sig2*rho;
			double c=-sig2*sqrt(1.-rho*rho);
			double b_shifted_v1=b+shift;
			double v0=  S2*NormalSpreadOptionIntegral_eta(a,b,c)+
						sig2*rho*NormalSpreadOptionIntegral_phi(a,b,c)+
						sig2*sqrt(1.-rho*rho)*NormalSpreadOptionIntegral_phi(a,c,b);
			double v1=  S2*NormalSpreadOptionIntegral_eta(a,b_shifted_v1,c)+
						sig2*rho*NormalSpreadOptionIntegral_phi(a,b_shifted_v1,c)+
						sig2*sqrt(1.-rho*rho)*NormalSpreadOptionIntegral_phi(a,c,b_shifted_v1);
			return (v1-v0)/shift;
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Normal_SpreadOption : incorrect optiontype toggle ");
		}
	}

	
}


double Vega2SpreadDigitalOption_N(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callOrPut,int optiontype)
{
	int nb_legendre_points=30;
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_N_Formula::DIGITALOPTION :
		{
			double standardDev=sqrt(sig1*sig1+sig2*sig2-2.*sig1*sig2*rho);
			double standardDev_Vega2=(sig2-rho*sig1)/standardDev;
			return VegaVanillaDigitalOption_N(S1-S2,standardDev,k,t,callOrPut)*standardDev_Vega2;
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::SPREADOPTION :
		{
			double standardDev=sqrt(sig1*sig1+sig2*sig2-2.*sig1*sig2*rho);
			double standardDev_Vega2=(sig2-rho*sig1)/standardDev;
			return VegaVanillaOption_N(S1-S2,standardDev,k,t,callOrPut)*standardDev_Vega2;
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::PAYFIRST :
		{	
			
			double shift=0.001;
			double a=S1-S2-k;
			double b=sig1-sig2*rho;
			double c=-sig2*sqrt(1.-rho*rho);
			double b_shifted_v2=b-rho*shift;
			double c_shifted_v2=c-sqrt(1.-rho*rho)*shift;

			double v0=	S1*NormalSpreadOptionIntegral_eta(a,b,c)+
						sig1*NormalSpreadOptionIntegral_phi(a,b,c);
			double v1=	S1*NormalSpreadOptionIntegral_eta(a,b_shifted_v2,c_shifted_v2)+
						sig1*NormalSpreadOptionIntegral_phi(a,b_shifted_v2,c_shifted_v2);
			return (v1-v0)/shift;
			break;
		}
	case ARM_CF_SpreadDigitalOption_N_Formula::PAYSECOND :
		{
			
			double shift=0.001;
			double a=S1-S2-k;
			double b=sig1-sig2*rho;
			double c=-sig2*sqrt(1.-rho*rho);
			double b_shifted_v2=b-rho*shift;
			double c_shifted_v2=c-sqrt(1.-rho*rho)*shift;
			double v0=	S2*NormalSpreadOptionIntegral_eta(a,b,c)+
						sig2*rho*NormalSpreadOptionIntegral_phi(a,b,c)+
						sig2*sqrt(1.-rho*rho)*NormalSpreadOptionIntegral_phi(a,c,b);
			double v1=	S2*NormalSpreadOptionIntegral_eta(a,b_shifted_v2,c_shifted_v2)+
						sig2*rho*NormalSpreadOptionIntegral_phi(a,b_shifted_v2,c_shifted_v2)+
						sig2*sqrt(1.-rho*rho)*NormalSpreadOptionIntegral_phi(a,c_shifted_v2,b_shifted_v2);
			return (v1-v0)/shift;
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Normal_SpreadOption : incorrect optiontype toggle ");
		}
	}

	
}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

