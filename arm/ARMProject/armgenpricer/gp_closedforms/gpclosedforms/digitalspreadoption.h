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
 
#ifndef _GP_CF_DIGITALSPREADOPTION_H
#define _GP_CF_DIGITALSPREADOPTION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/templated_pricing.h"
#include "gpclosedforms/generic_templated_pricing.h"
#include "gpclosedforms/digital_templated_pricing.h"



CC_BEGIN_NAMESPACE(ARM)

template<typename  smile,typename copula>
struct ARM_CF_DigitalSpreadOption_Formula : public Templated_Pricing<smile,copula>
{
	static double DigitalSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a10,double b10,double k10);
	static double DigitalSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   copula copula0,double t,int n,
				   double a10,double b10,double k10);

	static double DigitalSpreadOption_Pricing(int i, const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,
				   double a10,double b10,double k10);

};




///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Basic Pricing of  DigitalSpreadOption: implements the specific payoff :
///
///				cashflow=  1 if {a*S1+b*S2-k>0},  and 0 otherwise
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// avec copula copula_instance
template<class smile, class copula>
double ARM_CF_DigitalSpreadOption_Formula<smile,copula>::DigitalSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    copula copula_instance,double t,int n,double a10,double b10,double k10)
{
			
return Digital_Generic_Pricing(Underlying1,Underlying2,*(copula_instance.structure), t, n,1.0,a10,b10,k10);
	
}


/// avec ArgumentList&  copula
template<class smile, class copula>
double ARM_CF_DigitalSpreadOption_Formula<smile,copula>::DigitalSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				    const ArgumentList&  copula,double t,int n,double a10,double b10,double k10)
{
			return Digital_Generic_Pricing(Underlying1,Underlying2,copula, t, n,1.0,a10,b10,k10);	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///      Derivatives of Basic Pricing of  PowerSpreadOption: implements the specific payoff :
///
///				cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double ARM_CF_DigitalSpreadOption_Formula<smile,copula>::DigitalSpreadOption_Pricing(int i_deriv,const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,double a1,double b1,double k1)
{

			struct Ptr_DigitalSpreadOption_C : DoubleToDoubleFunc 
			{
				Ptr_DigitalSpreadOption_C(double a0,double b0,double k0) 
					: a(a0),b(b0),k(k0)) {}
				
				virtual double operator()(double x) const 
				{
						return -(a*x-k)/b ;
				}
				virtual double cash() const
				{
					return 1.0;
				}
				
			private:
				double a;
				double b;
				double k;
			};
			
			Ptr_DigitalSpreadOption_C payoff(a1,b1,k1);
			
			return Digital_Generic_Pricing(i,Underlying1,Underlying2,copula, t, n,k10,a20,b20,k20);
	
}



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

