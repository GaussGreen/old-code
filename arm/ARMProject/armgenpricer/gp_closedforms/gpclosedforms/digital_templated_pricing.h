/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file templated_pricing.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_DIGITAL_TEMPLATED_PRICING_H
#define _GP_CF_DIGITAL_TEMPLATED_PRICING_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>
#include <typeinfo>
#include <string>

#include "gpnumlib/gaussiananalytics.h"

#include "gpclosedforms/product.h"
#include "gpclosedforms/sum.h"
#include "gpclosedforms/difference.h"
#include "gpclosedforms/add_submodel.h"
#include "gpclosedforms/change_numeraire.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/smile_shiftedlognormal.h"
#include "gpclosedforms/gaussian_copula.h"
#include "gpclosedforms/ZaInverseLimit.h"
#include "gpclosedforms/templated_pricing.h"
#include "gpclosedforms/gaussian_copula.h"

#include "gpbase/numericconstant.h"


CC_BEGIN_NAMESPACE(ARM)

using std::string;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Templated_Pricing<smile,copula>::Digital_Generic_Pricing
///
///					computes the option
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
template<class smile,class copula>
double Templated_Pricing<smile,copula>::Digital_Generic_Pricing_1(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						double k1,double a2,double b2,double k2)
{
	 const char* typecopula=typeid(copula).name();		// the selection is done here by rtti ans string, because the partial template
														// feature is not implemented in VC6
	 std::string gaussiancopula("struct ARM::GaussianCopula");
	 std::string typecopulastr(typecopula);

	 if(typecopulastr == gaussiancopula)
	 {
		 ReducedGaussHermite_Coefficients c(n);
		 double Sum=0.,z_d;
		 int i;
		 double correl=copula_arg[0];
		 double h2,f2,n2,z,z_aux,z2,f2correl,sqrt1minuscorr;
		 sqrt1minuscorr=sqrt(1.-correl*correl);
		 for(i=0;i<n;i++)
		 {
			 z=c.get_point(i);
			 double xweight=c.get_weight(i);
			 z_d=ARM_GaussianAnalytics::cdfNormal(z);
			 if((z_d > 1e-15)&&(z_d<0.9999999999999999))
			 {
				 f2=smile::quantile(Underlying1,z_d,t);
				 f2correl=-(a2*f2-k2)/b2;
				 h2=smile::probability_distribution(Underlying2, f2correl, t);
				 
				 if((h2 > 1e-15)&&(h2<0.9999999999999999))
				 {
					 z2=ARM_GaussianAnalytics::cdfNormal_Inv(h2);
					 z_aux=(z2-correl*z)/sqrt1minuscorr;
					 n2=ARM_GaussianAnalytics::cdfNormal(z_aux);
					 Sum+=n2*xweight;
				 }
			 }
		 }
		 if(b2<0) return Sum;
		 else return 1.-Sum;
	 }
	 else
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Templated_Pricing<smile,copula>::Digital_Generic_Pricing_1   : bad copula");
	 }
}

*/

template<class smile,class copula>
double Templated_Pricing<smile,copula>::Digital_Generic_Pricing_1(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						double k1,double a2,double b2,double k2)
{
	 const char* typecopula=typeid(copula).name();		// the selection is done here by rtti ans string, because the partial template
														// feature is not implemented in VC6
	 std::string gaussiancopula("struct ARM::GaussianCopula");
	 std::string typecopulastr(typecopula);

	 if(typecopulastr == gaussiancopula)
	 {
		 ReducedGaussHermite_Coefficients c(n);
		 double Sum=0.;
		 int i;
		 double correl=copula_arg[0];
		 double f2,n2,z,z_aux,z2,f2correl,sqrt1minuscorr;
		 sqrt1minuscorr=sqrt(1.-correl*correl);
		 for(i=0;i<n;i++)
		 {
			 z=c.get_point(i);
			 double xweight=c.get_weight(i);
			 f2=smile::gaussian_to_distribution(Underlying1,z,t);
			 f2correl=-(a2*f2-k2)/b2;
			 if((smile::BoundedfromBelow()==1)&&(f2correl>= smile::LowerBoundary(Underlying2)))
			 {
				 z2=smile::distribution_to_gaussian(Underlying2,f2correl,t);
				 z_aux=(z2-correl*z)/sqrt1minuscorr;
				 n2=ARM_GaussianAnalytics::cdfNormal(z_aux);
				 Sum+=n2*xweight;
			 }
		 }
		 if(b2<0) return Sum;
		 else return 1.-Sum;
	 }
	 else
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Templated_Pricing<smile,copula>::Digital_Generic_Pricing_1   : bad copula");
	 }
}

template<class smile,class copula>
double Templated_Pricing<smile,copula>::Digital_Generic_Pricing(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						double k1,double a2,double b2,double k2)
{
	if(fabs(a2)<=fabs(b2))
	{
		return  Templated_Pricing<smile,copula>::Digital_Generic_Pricing_1(
						Underlying1,
						Underlying2,
						copula_arg,
						t,
						n,
						k1,a2,b2,k2);
	}
	else
	{
		return  Templated_Pricing<smile,copula>::Digital_Generic_Pricing_1(
						Underlying2,
						Underlying1,
						copula_arg,
						t,
						n,
						k1,b2,a2,k2);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Templated_Pricing<smile,copula>::Digital_Generic_Pricing_aux_First_Derivative
///
///					computes the option
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
template<class smile,class copula>
double Templated_Pricing<smile,copula>::Digital_Generic_Pricing_First_Derivative_2(
						int der_i,
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						double k1,double a2,double b2,double k2)
/// here i follows the conventions of the ARM_CF_SABR_VanillaOption_Formula :
///		FORWARD=0,
///		STRIKE=1,
///		MATURITY=2,   : not usable here
///		ALPHA=3,
///		BETA=4,
///		RHO=5,
///		NU=6
{
	 const char* typecopula=typeid(copula).name();		// the selection is done here by rtti ans string, because the partial template
														// feature is not implemented in VC6
	 std::string gaussiancopula("struct ARM::GaussianCopula");
	 std::string typecopulastr(typecopula);

	 if(typecopulastr == gaussiancopula)
	 {
		 ReducedGaussHermite_Coefficients c(n);
		 double Sum=0.,z_d,h2_deriv;
		 int i;
		 double correl=copula_arg[0];
		 double h2,f2,n2,z,z_aux,z_aux2,z2,f2correl,sqrt1minuscorr;
		 sqrt1minuscorr=sqrt(1.-correl*correl);
		 for(i=0;i<n;i++)
		 {
			 z=c.get_point(i);
			 z_d=ARM_GaussianAnalytics::cdfNormal(z);
			 if ((z_d > 1e-15)&&(z_d<0.9999999999999999))
			 {
				 f2=smile::quantile(Underlying1,z_d,t);
				 f2correl=-(a2*f2-k2)/b2;
				 h2=smile::probability_distribution(Underlying2, f2correl, t);
				 if ((h2 > 1e-15)&&(h2<0.9999999999999999))
				 {
					 h2_deriv=smile::probability_distribution_First_Derivative(der_i,Underlying2, f2correl, t);
					 z2=ARM_GaussianAnalytics::cdfNormal_Inv(h2);
					 z_aux=(z2-correl*z)/sqrt1minuscorr;
					 z_aux2=z_aux*z_aux;
					 n2=exp(-(z_aux2-z2*z2)/2.);
					 Sum+=n2*h2_deriv*c.get_weight(i);
				 }
			 }
		 }
		 if(b2<0) return Sum/sqrt1minuscorr;
		 else return -Sum/sqrt1minuscorr;
	 }
	 else
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Templated_Pricing<smile,copula>::Digital_Generic_Pricing_First_Derivative   : bad copula");
	 }
}
*/

						template<class smile,class copula>
double Templated_Pricing<smile,copula>::Digital_Generic_Pricing_First_Derivative_2(
						int der_i,
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						double k1,double a2,double b2,double k2)
/// here i follows the conventions of the ARM_CF_SABR_VanillaOption_Formula :
///		FORWARD=0,
///		STRIKE=1,
///		MATURITY=2,   : not usable here
///		ALPHA=3,
///		BETA=4,
///		RHO=5,
///		NU=6
{
	 const char* typecopula=typeid(copula).name();		// the selection is done here by rtti ans string, because the partial template
														// feature is not implemented in VC6
	 std::string gaussiancopula("struct ARM::GaussianCopula");
	 std::string typecopulastr(typecopula);

	 if(typecopulastr == gaussiancopula)
	 {
		 ReducedGaussHermite_Coefficients c(n);
		 double Sum=0.,h2_deriv;
		 int i;
		 double correl=copula_arg[0];
		 double f2,n2,z,z_aux,z_aux2,z2,f2correl,sqrt1minuscorr;
		 sqrt1minuscorr=sqrt(1.-correl*correl);
		 for(i=0;i<n;i++)
		 {
			 z=c.get_point(i);
			 f2=smile::gaussian_to_distribution(Underlying1,z,t);
			 f2correl=-(a2*f2-k2)/b2;
			 if((smile::BoundedfromBelow()==1)&&(f2correl>= smile::LowerBoundary(Underlying2)))
			 {
				 z2=smile::distribution_to_gaussian(Underlying2,f2correl,t);
				 h2_deriv=smile::probability_distribution_First_Derivative(der_i,Underlying2, f2correl, t);
				 z_aux=(z2-correl*z)/sqrt1minuscorr;
				 z_aux2=z_aux*z_aux;
				 n2=exp(-(z_aux2-z2*z2)/2.);
				 Sum+=n2*h2_deriv*c.get_weight(i);
			 }
		 }
		 if(b2<0) return Sum/sqrt1minuscorr;
		 else return -Sum/sqrt1minuscorr;
	 }
	 else
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Templated_Pricing<smile,copula>::Digital_Generic_Pricing_First_Derivative   : bad copula");
	 }
}


template<class smile,class copula>
double Templated_Pricing<smile,copula>::Digital_Generic_Pricing_First_Derivative_1(
						int der_i,
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						double k1,double a2,double b2,double k2)
{
	return  Templated_Pricing<smile,copula>::Digital_Generic_Pricing_First_Derivative_2(
						 der_i,
						 Underlying2,
						 Underlying1,
						 copula_arg,
						 t,
						 n,
						 k1, b2, a2, k2);
}



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

