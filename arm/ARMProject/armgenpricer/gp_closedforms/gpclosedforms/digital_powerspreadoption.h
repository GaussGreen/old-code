/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file digital_powerspreadoption.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_DIGITAL_POWERSPREADOPTION_H
#define _GP_CF_DIGITAL_POWERSPREADOPTION_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpclosedforms/templated_pricing.h"
#include "gpclosedforms/generic_templated_pricing.h"
#include "gpclosedforms/digital_templated_pricing.h"
//#include <iostream>
//#include <iomanip>
//#include <fstream>





CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  ARM_CF_DigitalPowerSpreadOption_Formula: implements the specific payoff :
///
///				cashflow=  1 if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename  smile,typename copula>
struct ARM_CF_DigitalPowerSpreadOption_Formula : public Templated_Pricing<smile,copula>
{
	////////////////////////////////////////////////////////////
	///
	///		external interface
	///
	////////////////////////////////////////////////////////////
	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a20,double b20,double k20);

	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const copula copula_instance,double t,int n,
				   double a20,double b20,double k20);
	////////////////////////////////////////////////////////////
	///
	///		Simple digital spreadoption : payoff 1 if S1-S2-K>0
	///
	////////////////////////////////////////////////////////////
	static double Simple_Digital_Pricing_aux(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double K);
	static double Simple_Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double K);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula: implements the specific payoff :
///
///				cashflow=  S2 if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// 
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename  smile,typename copula>
struct ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula : public Templated_Pricing<smile,copula>
{
	////////////////////////////////////////////////////////////
	///
	///		external interface
	///
	////////////////////////////////////////////////////////////
	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a20,double b20,double k20);

	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const copula copula_instance,double t,int n,
				   double a20,double b20,double k20);
	////////////////////////////////////////////////////////////
	///
	///		Simple digital spreadoption : payoff 1 if S1-S2-K>0
	///
	////////////////////////////////////////////////////////////
	static double Simple_Digital_Pricing_aux(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double K);
	static double Simple_Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double K);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula: implements the specific payoff :
///
///				cashflow=  S1 if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// used by spreadoption formula for SABR distribution because the integrated formula given by change of numeraire do not exist
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename  smile,typename copula>
struct ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula : public Templated_Pricing<smile,copula>
{
	////////////////////////////////////////////////////////////
	///
	///		external interface
	///
	////////////////////////////////////////////////////////////
	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a20,double b20,double k20);

	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const copula copula_instance,double t,int n,
				   double a20,double b20,double k20);
	////////////////////////////////////////////////////////////
	///
	///		Simple digital spreadoption : payoff 1 if S1-S2-K>0
	///
	////////////////////////////////////////////////////////////
	static double Simple_Digital_Pricing_aux(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double K);
	static double Simple_Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double K);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  ARM_CF_GSpreadOption_Formula: implements the specific payoff :
///
///				cashflow=  S1-S2-K if {S1-S2-k2>0},  and 0 otherwise
///
/// 
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename  smile,typename copula>
struct ARM_CF_GSpreadOption_Formula : public Templated_Pricing<smile,copula>
{
	////////////////////////////////////////////////////////////
	///
	///		external interface
	///
	////////////////////////////////////////////////////////////
	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double k20);
	static double Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double k20);


};




///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  1 if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// used by spreadoption formula for SABR distribution because the integrated formula given by change of numeraire do not exist
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double a20,double b20,double k20)
{		
			return Digital_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,a20,b20,k20);
}


///  introduce the homothetic transformations to transform the problem a*S1+b*S2 -K into S1-S1 -K

template<class smile, class copula>
double ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula,double t,int n,double a20,double b20,double k20)
{
	if ((fabs(a20)<1e-20) || (fabs(b20)<1e-20))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing   : a20 or b20 almost 0 !");
		
	}
	ArgumentList* Renormalized_Underlying1=smile::HomotheticTransformation(&Underlying1,fabs(a20));
	ArgumentList* Renormalized_Underlying2=smile::HomotheticTransformation(&Underlying2,fabs(b20));

	if((a20>0) && (b20<0))
	{
		return Simple_Digital_Pricing(*Renormalized_Underlying1,*Renormalized_Underlying2,copula, t, n,k20);
	}
	else if ((b20>0) && (a20<0))
	{
		return Simple_Digital_Pricing(*Renormalized_Underlying2,*Renormalized_Underlying1,copula, t, n,k20);
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing   : S1 + S2 -K> not yet implemented !");
		
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
///
///
/// This is where the action is for computation of digital options, for K>0 only
///	 K is changed of sign with respect to the mathamatica implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
/*
template<class smile, class copula>
double ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing_aux(
					const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula0,double t,int n, double K)
{
	double T=1.;
	double rho=copula0[0];
	double origlimK=copula::marginal_quantile(copula0,smile::probability_distribution(Underlying2,K,T),T);
	double limK=(origlimK>-10.)? origlimK: -10.;
	ReducedGaussHermite_Coefficients c(n);
	double Sum=0.;
	int i;
	double xpoint,xweight,I1,I2,I3,x1,x2,x3,x4,z;
	if (limK>=0)
	{
		for(i=0;i<n;i++)
		{
			xpoint=c.get_point(i);
			xweight=c.get_weight(i);
			Sum+=copula::Q_Integral(copula0,limK,xpoint)*exp(xpoint*xpoint/2.)*xweight;
		}
		I1=Sum/ARM_NumericConstants::ARM_INVSQRT2PI+0.5;
	}
	else
	{	
		for(i=0;i<n;i++)
		{
			xpoint=c.get_point(i);
			xweight=c.get_weight(i);
			Sum+=copula::Q_Integral(copula0,-limK,xpoint)*exp(xpoint*xpoint/2.)*xweight;
		}
		I1=0.5-Sum/ARM_NumericConstants::ARM_INVSQRT2PI;
	}
	if(limK>=6.)
	{
		I2=0.;
		I3=0.;
	}
	else
	{
		GaussLegendre_Coefficients d(n, limK, 6.);
		Sum=0;
		double oneoverrac=1./sqrt(1.-rho*rho),isum;
		for(i=0;i<n;i++)
		{
			xpoint=d.get_point(i);
			xweight=d.get_weight(i);
			Sum+=copula::Q_Total_Integral(copula0,xpoint)*xweight;
		}
		I2=Sum;
		Sum=0;

		for(i=0;i<n;i++)
		{
			xpoint=d.get_point(i);
			xweight=d.get_weight(i);
			x1=copula::marginal_distribution(copula0,xpoint,0);
			x2=smile::quantile(Underlying2,x1,0)-K;
			x3=smile::probability_distribution(Underlying1,x2,0);
			x4=copula::marginal_quantile(copula0,x3,0);
			z=oneoverrac*(x4-rho*xpoint);
			isum=copula::Q_Integral(copula0,z,xpoint);
			Sum+=isum*xweight;

		}

		I3=Sum;
		
	}
		
		return 1.-(I1+I2-I3);
}
*/
template<class smile, class copula>
double ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing_aux(
					const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula0,double t,int n, double K)
{
	double T=1.;
	double rho=copula0[0];
	double origlimK=copula::marginal_quantile(copula0,smile::probability_distribution(Underlying2,K,T),T);
	double limK=(origlimK>-10.)? origlimK: -10.;
	double Sum=0.;
	int i;
	double xpoint,xweight,I3,x1,x2,x3,x4,z;	
	GaussLegendre_Coefficients d(n, -10, 6.);
	Sum=0;
		double oneoverrac=1./sqrt(1.-rho*rho),isum;
	for(i=0;i<n;i++)
	{
		xpoint=d.get_point(i);
		xweight=d.get_weight(i);
		x1=copula::marginal_distribution(copula0,xpoint,0);
		x2=smile::quantile(Underlying2,x1,0)-K;
		if(x2<=0) 
		{
			x3=0;
		}
		else
		{
			x3=smile::probability_distribution(Underlying1,x2,0);
		}
		if(x3>=0.9999999999)
		{
			isum=copula::Q_Total_Integral(copula0,xpoint);
		}
		else
		{
			x4=copula::marginal_quantile(copula0,x3,0);
			z=oneoverrac*(x4-rho*xpoint);
			isum=copula::Q_Integral(copula0,z,xpoint);	
		}
		Sum+=isum*xweight;
	}
	I3=Sum;
	return 0.5-I3;
}

template<class smile, class copula>
double ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula0,double t,int n, double K)
{
	if (K>=0.)
	{
		return ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing_aux(
					Underlying1,Underlying2,copula0, t, n,  K);
	}
	else
	{
		return 1.-ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing_aux(
					Underlying2,Underlying1,copula0, t, n,  -K);
	}

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  S2 if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// used by spreadoption formula for SABR distribution because the integrated formula given by change of numeraire do not exist
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double a20,double b20,double k20)
{		
			return Digital_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,a20,b20,k20);
}



///  introduce the homothetic transformations to transform the problem a*S1+b*S2 -K into S1-S1 -K

template<class smile, class copula>
double ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula,double t,int n,double a20,double b20,double k20)
{
	if ((fabs(a20)<1e-20) || (fabs(b20)<1e-20))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing   : a20 or b20 almost 0 !");
		
	}
	ArgumentList* Renormalized_Underlying1=smile::HomotheticTransformation(&Underlying1,fabs(a20));
	ArgumentList* Renormalized_Underlying2=smile::HomotheticTransformation(&Underlying2,fabs(b20));

	if((a20>0) && (b20<0))
	{
		return Simple_Digital_Pricing(*Renormalized_Underlying1,*Renormalized_Underlying2,copula, t, n,k20);
	}
	else if ((b20>0) && (a20<0))
	{
		return Simple_Digital_Pricing(*Renormalized_Underlying2,*Renormalized_Underlying1,copula, t, n,k20);
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing   : S1 + S2 -K> not yet implemented !");
		
	}
}

template<class smile, class copula>
double ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula0,double t,int n, double K)
{
	double T=1.;
	double rho=copula0[0];
	GaussLegendre_Coefficients d(n, -1., 1.);
	double Sum=0.;
	int i;
	double xpoint,xweight,I1,I2,I3,x1,x2,x3,x4,z,limK;
	if ( K<0)
	{
		I1=0;
		GaussLegendre_Coefficients d2(&d, -10., 6.);
		Sum=0;
		double oneoverrac=1./sqrt(1.-rho*rho),isum;
		for(i=0;i<n;i++)
		{
			xpoint=d2.get_point(i);
			xweight=d2.get_weight(i);
			x1=copula::marginal_distribution(copula0,xpoint,0);
			x2=smile::quantile(Underlying2,x1,0);
			Sum+=x2*copula::Q_Total_Integral(copula0,xpoint)*xweight;
		}
		I2=Sum;
		Sum=0;
		for(i=0;i<n;i++)
		{
			xpoint=d2.get_point(i);
			xweight=d2.get_weight(i);
			x1=copula::marginal_distribution(copula0,xpoint,0);
			x2=smile::quantile(Underlying2,x1,0)-K;
			if(x2<=0)
			{
				x3=0;
			}
			else
			{
				x3=smile::probability_distribution(Underlying1,x2,0);
			}
			if(x3>=0.9999999999)
			{
				isum=copula::Q_Total_Integral(copula0,xpoint);
			}
			else
			{
				x4=copula::marginal_quantile(copula0,x3,0);
				z=oneoverrac*(x4-rho*xpoint);
				isum=copula::Q_Integral(copula0,z,xpoint);
			}
			x2+=K;;
			Sum+=x2*isum*xweight;
		}
		I3=Sum;
		
		
	}
	else			/// K>0
	{
		limK=copula::marginal_quantile(copula0,smile::probability_distribution(Underlying2,K,T),T);
		limK=(limK>-10.)? limK: -10.;
		Sum=0.;
		GaussLegendre_Coefficients d1(&d, -10., limK);
		for(i=0;i<n;i++)
		{
			xpoint=d1.get_point(i);
			xweight=d1.get_weight(i);
			x1=copula::marginal_distribution(copula0,xpoint,0);
			x2=smile::quantile(Underlying2,x1,0);
			Sum+=x2*copula::Q_Total_Integral(copula0,xpoint)*xweight;
		}
		I1=2.*Sum;
		GaussLegendre_Coefficients d2(&d, limK, 6.);
		Sum=0;
		double oneoverrac=1./sqrt(1.-rho*rho),isum;
		for(i=0;i<n;i++)
		{
			xpoint=d2.get_point(i);
			xweight=d2.get_weight(i);
			x1=copula::marginal_distribution(copula0,xpoint,0);
			x2=smile::quantile(Underlying2,x1,0);
			Sum+=x2*copula::Q_Total_Integral(copula0,xpoint)*xweight;
		}
		I2=Sum;
		Sum=0;
		for(i=0;i<n;i++)
		{
			xpoint=d2.get_point(i);
			xweight=d2.get_weight(i);
			x1=copula::marginal_distribution(copula0,xpoint,0);
			x2=smile::quantile(Underlying2,x1,0)-K;
			if(x2<0) 
			{
				isum=0;
			}
			else
			{
				x3=smile::probability_distribution(Underlying1,x2,0);
				if(x3>0.99999999)
				{
					isum=copula::Q_Total_Integral(copula0,xpoint);
				}
				else
				{
					x4=copula::marginal_quantile(copula0,x3,0);
					z=oneoverrac*(x4-rho*xpoint);
					isum=copula::Q_Integral(copula0,z,xpoint);
				}
			}
			Sum+=(x2+K)*isum*xweight;
		}
		I3=Sum;
		
	}
	
	return I1+I2-I3;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  S1 if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// used by spreadoption formula for SABR distribution because the integrated formula given by change of numeraire do not exist
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double a20,double b20,double k20)
{		
			return Digital_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,a20,b20,k20);
}


///  introduce the homothetic transformations to transform the problem a*S1+b*S2 -K into S1-S1 -K

template<class smile, class copula>
double ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula,double t,int n,double a20,double b20,double k20)
{
	if ((fabs(a20)<1e-20) || (fabs(b20)<1e-20))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing   : a20 or b20 almost 0 !");
		
	}
	ArgumentList* Renormalized_Underlying1=smile::HomotheticTransformation(&Underlying1,fabs(a20));
	ArgumentList* Renormalized_Underlying2=smile::HomotheticTransformation(&Underlying2,fabs(b20));

	if((a20>0) && (b20<0))
	{
		return Simple_Digital_Pricing(*Renormalized_Underlying1,*Renormalized_Underlying2,copula, t, n,k20);
	}
	else if ((b20>0) && (a20<0))
	{
		return Simple_Digital_Pricing(*Renormalized_Underlying2,*Renormalized_Underlying1,copula, t, n,k20);
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula<smile,copula>::Digital_Pricing   : S1 + S2 -K> not yet implemented !");
		
	}
}

template<class smile, class copula>
double ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula0,double t,int n, double K)
{
	///  We recompute the expectation of the First index, but this one can be available somewhere ...
	double T=1.;
	double rho=copula0[0];
	GaussLegendre_Coefficients d(n, 0., 1.);
	double Sum=0.;
	int i;
	double xpoint,xweight,S1;
	for(i=0;i<n;i++)
	{
		xpoint=d.get_point(i);
		xweight=d.get_weight(i);
		Sum+=smile::quantile(Underlying1,xpoint,0)*xweight;
	}
	S1=Sum;
	/// We use the correspondnance with the S2 analytics:
	double Complement=ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					  Underlying2,Underlying1,copula0, t, n, -K);
	return S1-Complement;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  SpreadOption: implements the specific payoff :
///
///				cashflow=  S1-S2-K if {S1-S2-k2>0},  and 0 otherwise
///
/// used by spreadoption formula for SABR distribution because the integrated formula given by change of numeraire do not exist
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_GSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double k20)
{				
			return Digital_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,k20);
}

template<class smile, class copula>
double ARM_CF_GSpreadOption_Formula<smile,copula>::Digital_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList& copula0,double t,int n,double K)
{		
			double payS1= ARM_CF_Index1PayingDigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					Underlying1, Underlying2,copula0, t, n,  K);
			double payS2= ARM_CF_Index2PayingDigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					Underlying1, Underlying2,copula0, t, n,  K);
			double payK=-K*ARM_CF_DigitalPowerSpreadOption_Formula<smile,copula>::Simple_Digital_Pricing(
					Underlying1,Underlying2,copula0, t, n, K);
			return payS1-payS2-payK;
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

