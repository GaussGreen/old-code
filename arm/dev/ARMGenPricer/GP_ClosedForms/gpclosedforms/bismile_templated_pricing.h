/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file powerspreadoption.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_BISMILE_TEMPLATEDPRICING_H
#define _GP_CF_BISMILE_TEMPLATEDPRICING_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/templated_pricing.h"
#include "gpclosedforms/generic_templated_pricing.h"
#include "gpclosedforms/digital_templated_pricing.h"



CC_BEGIN_NAMESPACE(ARM)

template<typename  smile1,typename smile2,typename copula>
struct Bismiled_Generic_Pricing
{
	///   Generic methods (relative to the smile distribution and the copula

	static double Generic_Pricing_aux(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);
	static double Generic_Pricing_aux2(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);



	/// used					in ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing
	/// and therefore  used		by	ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula::value
	static double Generic_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);
	static double Generic_Pricing2(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);


};







///////////////////////////////////////////////////////////////////////////
///
///
///
///
///			Two pricing algo, 
///
///				Generic_Pricing_aux  et Generic_Pricing_aux2 does the job
///
///////////////////////////////////////////////////////////////////////////

template<class smile1,class smile2,class copula>
double Bismiled_Generic_Pricing<smile1,smile2,copula>::Generic_Pricing_aux(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						BivariateToDoubleFunc* option_payoff)
{
	double left_limit=copula::marginal_left_limit(copula_arg);
	double right_limit=copula::marginal_right_limit(copula_arg);
	GaussLegendre_Coefficients c(n,left_limit,right_limit);
	double Sum=0;
	double Sum0=0;
	double x,y,x2,y2,x3,y3;
	int i,j;
	for(i=0;i<n;i++){
		x=c.get_point(i);
		x2=copula::marginal_distribution(copula_arg,x,t);
		x3=smile1::quantile(Underlying1,x2,t);
		for(j=0;j<n;j++){
			y=c.get_point(j);
			y2=copula::marginal_distribution(copula_arg,y,t);
			y3=smile2::quantile(Underlying2,y2, t);
			double dens=copula::multivariate_density(copula_arg,x,y);;
			Sum+= (*option_payoff)(x3,y3)*dens*c.get_weight(i)*c.get_weight(j);
			Sum0+=c.get_weight(i)*c.get_weight(j);	
		}	
	}
	return Sum;
}


template<class smile1,class smile2,class copula>
double Bismiled_Generic_Pricing<smile1,smile2,copula>::Generic_Pricing_aux2(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						BivariateToDoubleFunc* option_payoff)
{
	GaussHermite_Coefficients c(n);
	double Sum=0;
	double Sum0=0;
	double x,y,ycorrelated,f1,f2,paid;
	double Sum00=0.0;
	int i,j;
	for(i=0;i<n;i++) for(j=0;j<n;j++) Sum00 +=c.get_weight(i)*c.get_weight(j);
	
	for(i=0;i<n;i++){
		x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2;
		if((x<4.)&&(x>-4))
		{
			f1=smile1::gaussian_to_distribution(Underlying1, x, t);
			for(j=0;j<n;j++){
				
				y=c.get_point(j)*ARM_NumericConstants::ARM_SQRT_2;
				ycorrelated=copula::gaussian_mix(copula_arg,x,y);
				if((ycorrelated<4.)&&(ycorrelated>-4.))
				{
					f2=smile2::gaussian_to_distribution(Underlying2, ycorrelated, t);
					paid=(*option_payoff)(f1,f2);
					Sum+= paid*c.get_weight(i)*c.get_weight(j);
				}
				Sum0+=c.get_weight(i)*c.get_weight(j);	
			}
		}
		else
		{
			Sum0+=c.get_weight(i);	
		}
	}
	return Sum/ARM_NumericConstants::ARM_PI;
}

template<class smile1,class smile2, class copula>
double Bismiled_Generic_Pricing<smile1,smile2,copula>::Generic_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff)
{
	Swappedptr_funBivariateToDouble swp_payoff(option_payoff);		/// if necessary
	
	return Bismiled_Generic_Pricing<smile1,smile2,copula>::Generic_Pricing_aux(Underlying1,Underlying2,copula_arg,t,n,option_payoff);
	
}

template<class smile1,class smile2, class copula>
double Bismiled_Generic_Pricing<smile1,smile2,copula>::Generic_Pricing2(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff)
{
	Swappedptr_funBivariateToDouble swp_payoff(option_payoff);		/// if necessary
	
	return Bismiled_Generic_Pricing<smile1,smile2,copula>::Generic_Pricing_aux2(Underlying1,Underlying2,copula_arg,t,n,option_payoff);
	
}


///    




template<typename  smile1,typename  smile2,typename copula>
struct ARM_CF_Bismiled_PowerSpreadOption_Formula : public Bismiled_Generic_Pricing<smile1,smile2,copula>
{
	static double Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const copula& copula0,double t,int n,
				   double a10,double b10,double k10,double a20,double b20,double k20,int algorithm);

};




///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
/// used by spreadoption formula for SABR distribution because the integrated formula given by change of numeraire do not exist
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile1,class smile2, class copula>
double ARM_CF_Bismiled_PowerSpreadOption_Formula<smile1,smile2,copula>::Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				const copula& copula_instance,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20,int algorithm)
{

	const char* typecopula=typeid(copula).name();		// the selection is done here by rtti and string, because the partial template
	// feature is not implemented in VC6
	std::string gaussiancopula("struct ARM::GaussianCopula");
	std::string typecopulastr(typecopula);
	
	switch (algorithm)
	{
	case GENERIC_SPREADOPTION :
		{
			struct Ptr_PowerSpreadOption_C : BivariateToDoubleFunc 
			{
				Ptr_PowerSpreadOption_C(double a100,double b100,double k100,double a200,double b200,double k200) 
					: a1(a100),a2(a200),b1(b100),b2(b200),k1(k100),k2(k200) {}
				
				virtual double operator()(double s1,double s2) const 
				{
					double cashflow=a1*s1+b1*s2-k1;
					double decision=a2*s1+b2*s2-k2;
					double kk1=k1;
					double kk2=k2;
					double aa1=a1;
					double aa2=a2;
					double bb1=b1;
					double bb2=b2;
					if((a1*s1+b1*s2-k1>=0) && (a2*s1+b2*s2-k2 >=0))
					{
						return a1*s1+b1*s2-k1;
					}
					else
					{
						return 0.;
					}
				}
				
			private:
				double a1;
				double b1;
				double k1;
				double a2;
				double b2;
				double k2;
			};
			
			Ptr_PowerSpreadOption_C payoff(a10,b10,k10,a20,b20,k20);
			double value,intrinseq;
			double s1=Underlying1[0];
			double s2=Underlying2[0];

			if (typecopulastr == gaussiancopula)
			{
				value= Generic_Pricing2(Underlying1,Underlying2,*(copula_instance.structure), t, n,&payoff);
			}
			else
			{
				value= Generic_Pricing2(Underlying1,Underlying2,*(copula_instance.structure), t, n,&payoff);
			}
			if((a10*s1+b10*s2-k10>=0)&&(a20==a10)&&(b20==b10)&&(k20==k10))
			{
				intrinseq= a10*s1+b10*s2-k10;
			}
			else
			{
				intrinseq= 0.;
			}
			if(value>intrinseq) return value; else return intrinseq;
			
		}
	case DIGITAL_SPREADOPTION :
		{	
			if (typecopulastr == gaussiancopula)
			{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"ARM_CF_Bismiled_PowerSpreadOption_Formula<smile1,smile2,copula>::Pricing :  gaussian copula for digital ,not yet implemented ");
			break;
			}
			else
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"ARM_CF_Bismiled_PowerSpreadOption_Formula<smile1,smile2,copula>::Pricing : non gaussian copula for digital ,not yet implemented ");
			break;
			}
			
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"ARM_CF_Bismiled_PowerSpreadOption_Formula<smile1,smile2,copula>::Pricing : ALGORITHM : bad value ");
			break;
		}
	}
	
}



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


