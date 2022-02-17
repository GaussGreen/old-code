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
 
#ifndef _GP_CF_POWERSPREADOPTION_H
#define _GP_CF_POWERSPREADOPTION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/templated_pricing.h"
#include "gpclosedforms/generic_templated_pricing.h"
#include "gpclosedforms/digital_templated_pricing.h"



CC_BEGIN_NAMESPACE(ARM)

enum Spreadoption_AlgorithmType
	{
		GENERIC_SPREADOPTION,
		DIGITAL_SPREADOPTION
	};


template<typename  smile,typename copula>
struct ARM_CF_PowerSpreadOption_Formula : public Templated_Pricing<smile,copula>
{
	static double PowerSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a10,double b10,double k10,double a20,double b20,double k20,int algorithm);
	static double PowerSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   copula copula0,double t,int n,
				   double a10,double b10,double k10,double a20,double b20,double k20,int algorithm);

	static double PowerSpreadOption_Pricing(int i, const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a10,double b10,double k10,double a20,double b20,double k20,int algorithm);
	static double PowerSpreadOption_Pricing_With_Limits(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
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
template<class smile, class copula>
double ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20,int algorithm)
{
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

			const char* typecopula=typeid(copula).name();		// the selection is done here by rtti ans string, because the partial template
														// feature is not implemented in VC6
			std::string gaussiancopula("struct ARM::GaussianCopula");
			std::string typecopulastr(typecopula);

			
			if (typecopulastr == gaussiancopula)
			{
				return Generic_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,&payoff);
			}
			else
			{
				return Generic_Pricing_GeneralCopula(Underlying1,Underlying2,*(copula_instance.structure), t, n,&payoff);
			}

		}
	case DIGITAL_SPREADOPTION :
		{		
			return Digital_Generic_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,k10,a20,b20,k20);
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing : ALGORITHM : bad value ");
			break;
		}
	}
	
}


/// used by ARM_CF_PowerSpreadOption_Formula<ShiftedLogNormal_Smile,GaussianCopula>

template<class smile, class copula>
double ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20,int algorithm)
{
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
			
			return Generic_Pricing(Underlying1,Underlying2,copula, t, n,&payoff);	
			
		}
	case DIGITAL_SPREADOPTION :
		{		
			return Digital_Generic_Pricing(Underlying1,Underlying2,copula, t, n,k10,a20,b20,k20);
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing : ALGORITHM : bad value ");
			break;
		}
	}
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///      Derivatives of Basic Pricing of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing(int i_deriv,const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20,int algorithm)
{
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
			
			return Generic_Pricing(Underlying1,Underlying2,copula, t, n,&payoff);
		}
	case DIGITAL_SPREADOPTION :
		{
			struct Ptr_PowerSpreadOption_C : DoubleToDoubleFunc 
			{
				Ptr_PowerSpreadOption_C(double k100,double a200,double b200,double k200) 
					: a2(a200),b2(b200),k1(k100),k2(k200) {}
				
				virtual double operator()(double x) const 
				{
						return -(a2*x-k2)/b2 ;
				}
				virtual double cash() const
				{
					return k1;
				}
				
			private:
				double k1;
				double a2;
				double b2;
				double k2;
			};
			
			Ptr_PowerSpreadOption_C payoff(k10,a20,b20,k20);
			
			return Digital_Generic_Pricing(i,Underlying1,Underlying2,copula, t, n,k10,a20,b20,k20);


		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing : ALGORITHM : bad value ");
			break;
		}
	}
	
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing N2 of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
///		provides an information on when option =0 ( S2 < Ldown(S1)  )  through the function Ldown(S1)

///   used by ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing_With_Limits (const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20,int algorithm)
{
	struct Ptr_PowerSpreadOption_PayOff : BivariateToDoubleFunc 
	{
		Ptr_PowerSpreadOption_PayOff(double a100,double b100,double k100,double a200,double b200,double k200) 
			: a1(a100),a2(a200),b1(b100),b2(b200),k1(k100),k2(k200) {}

		virtual double operator()(double s1,double s2) const 
		{
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
	struct Ptr_PowerSpreadOption_Limit_S2 : UnivariateToDoubleFunc 
	{
		Ptr_PowerSpreadOption_Limit_S2(double a100,double b100,double k100,double a200,double b200,double k200) 
			: a1(a100),a2(a200),b1(b100),b2(b200),k1(k100),k2(k200) {}

		virtual double operator()(double S1) const 
		{	
			return (k2-a2*S1)/b2;
		}
	
		int IntegrationDirection() const
		{
			if (b2>0)
				return UnivariateToDoubleFunc::UPWARD;
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			///
			///		because in the lognormal case gaussian_to_distribution(x) is an increasing function
			//		and GaussianCopula::gaussian_mix is an increasing fucntion so in summary 
			///		smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,z), t) is an increasing function of z
			///		so when S2> (k2-a*S1)/b, we have z<za where 
			///		smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,z), t)=(k2-a*S1)/b
			///
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			else
				return UnivariateToDoubleFunc::DOWNWARD;
		}
		
	private:
		double a1;
		double b1;
		double k1;
		double a2;
		double b2;
		double k2;
	};
	struct Ptr_PowerSpreadOption_Limit_S1 : UnivariateToDoubleFunc 
	{
		Ptr_PowerSpreadOption_Limit_S1(double a100,double b100,double k100,double a200,double b200,double k200) 
			: a1(a100),a2(a200),b1(b100),b2(b200),k1(k100),k2(k200) {}

		virtual double operator()(double S2) const 
		{	
			return (k2-b2*S2)/a2;
		}
		int IntegrationDirection() const
		{
			if (a2>0)
				return UnivariateToDoubleFunc::DOWNWARD;
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			///
			///		because in the lognormal case gaussian_to_distribution(x) is a decreasing function
			//		and GaussianCopula::gaussian_mix is an increasing fucntion so in summary 
			///		smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,z), t) is a decreasing function of z
			///		so when S2> (k2-a*S1)/b, we have z<za where 
			///		smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,z), t)=(k2-a*S1)/b
			///
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			else
				return UnivariateToDoubleFunc::UPWARD;
		}
	private:
		double a1;
		double b1;
		double k1;
		double a2;
		double b2;
		double k2;
	};
	
	Ptr_PowerSpreadOption_PayOff payoff(a10,b10,k10,a20,b20,k20);
	if((b20==0)&&(a20==0))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing2   : bad arguments : a2 or b2");
	}
	else
	{
		if
			(b20!=0)	// l'integrale sur S2 est interne et limite par ls
		{
			Ptr_PowerSpreadOption_Limit_S2 ls(a10,b10,k10,a20,b20,k20);
			return Generic_Pricing_S2(Underlying1,Underlying2,copula, t, n,&payoff,&ls);
		}
		else		// So a2 should be here !=0	 , ici l'integrale sur S1 est interne et limite par ls
		{
			Ptr_PowerSpreadOption_Limit_S1 ls(a10,b10,k10,a20,b20,k20);
			return Generic_Pricing_S1(Underlying1,Underlying2,copula, t, n,&payoff,&ls);
		}
	}
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

