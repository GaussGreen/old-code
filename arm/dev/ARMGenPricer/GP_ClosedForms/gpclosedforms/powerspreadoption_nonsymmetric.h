/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file powerspreadoption_nonsymmetric.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_POWERSPREADOPTION_NONSYMMETRIC_H
#define _GP_CF_POWERSPREADOPTION_NONSYMMETRIC_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/templated_pricing.h"
#include "gpclosedforms/templated_pricing_nonsymmetric.h"
#include "gpclosedforms/generic_templated_pricing.h"
#include "gpclosedforms/digital_templated_pricing.h"



CC_BEGIN_NAMESPACE(ARM)



///********************************************************** Templated pricing of a gaussian power spreadoption *****************
///
///
///
///
///****************************************************************************************************************
template<typename  smile1,typename  smile2>
struct ARM_CF_PowerSpreadOption_GuaussianNonSymmetric_Formula : public Templated_Pricing_NonSymmetricGaussian<smile1,smile2>
{
	static double PowerSpreadOption_GaussianPricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
		double copula_correlation,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20);

};


template<class smile1,class smile2>
double ARM_CF_PowerSpreadOption_GuaussianNonSymmetric_Formula<smile1,smile2>::PowerSpreadOption_GaussianPricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   double copula_correlation,double t,int n,double a10,double b10,double k10,double a20,double b20,double k20)
{
	struct Ptr_PowerSpreadOption_PayOff : BivariateToDoubleFunc 
	{
		Ptr_PowerSpreadOption_PayOff(double a100,double b100,double k100,double a200,double b200,double k200) 
			: a1(a100),a2(a200),b1(b100),b2(b200),k1(k100),k2(k200) {}

		virtual double operator()(double s1,double s2) const 
		{
			double x=a1*s1+b1*s2-k1;
			if(x>=0)
			{
				return x;
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
	Ptr_PowerSpreadOption_PayOff payoff(a10,b10,k10,a20,b20,k20);
	double K=k20/a20;
	double epsilon=b20/a20;
	
	return Templated_Pricing_NonSymmetricGaussian<smile1,smile2>::Generic_SpreadoptionPricing(
		Underlying1,
		Underlying2,
		copula_correlation,
		K,
		t,
		n,
		epsilon,a20,
		&payoff);
	
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

