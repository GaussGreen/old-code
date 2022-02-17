/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file templated_pricing_nonsymmetric.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_TEMPLATED_PRICING_NONSYMMETRIC_H
#define _GP_CF_TEMPLATED_PRICING_NONSYMMETRIC_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/product.h"
#include "gpclosedforms/sum.h"
#include "gpclosedforms/difference.h"
#include "gpclosedforms/add_submodel.h"
#include "gpclosedforms/change_numeraire.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/smile_shiftedlognormal.h"
#include "gpclosedforms/gaussian_copula.h"
#include "gpclosedforms/student_copula.h"
#include "gpclosedforms/ZaInverseLimit.h"
#include "gpclosedforms/templated_pricing.h"

#include "gpbase/numericconstant.h"


CC_BEGIN_NAMESPACE(ARM)


/// #define ARM_CF_INVSQRTPI 0.5641895835477563


///********************************************************** Templated pricing of a gaussian spreadoption *****************
///
///
///
///
///****************************************************************************************************************
template<typename  smile1,typename smile2>
struct Templated_Pricing_NonSymmetricGaussian
{
	///   Generic methods (relative to the smile distribution and the copula
	static double Generic_SpreadoptionPricing(
								const ArgumentList& Underlying1,
								const ArgumentList& Underlying2,
								double copula_correlation,
								double K,
								double t,
								int n,
								double epsilon,double a20,
								BivariateToDoubleFunc* option_payoff); /// payoff  if(S1 -epsilon*S2-K>0)
	
};	
	
template<class smile1,class smile2>
double Templated_Pricing_NonSymmetricGaussian<smile1,smile2>::Generic_SpreadoptionPricing(
								const ArgumentList& Underlying1,
								const ArgumentList& Underlying2,
								double copula_correlation,
								double K,
								double t,
								int n,
								double epsilon,
								double a20,  /// gives the direction of integration
								BivariateToDoubleFunc* option_payoff)
{
	ReducedGaussHermite_Coefficients c(n);
	GaussLegendre_Coefficients c2(n);
	double Sum=0;
	double x,y,Sy,Sx,xcorrelated;
	double sqrt_1_rho2=sqrt(1.0-copula_correlation*copula_correlation);
			for(int i=0;i<n;i++)
			{
				y=c.get_point(i);
				if((y<4.)&&(y>-4))
				{
					Sy=smile2::gaussian_to_distribution(Underlying2, y, t);
					// we assume the the support for x is {-infinity,infinity}
					// we want to integrate from -epsilon*Sy+K  to xmax or xmin depending on a20 

			
					for(int j=0;j<n;j++)  /// only for x greater than xlimit
					{
						x=c.get_point(j);
						xcorrelated=x*sqrt_1_rho2+copula_correlation*y;
						if((xcorrelated<4.)&&(xcorrelated>-4.))
						{
							Sx=smile1::gaussian_to_distribution(Underlying1, xcorrelated, t);
							double pay=(*option_payoff)(Sx,Sy);
							Sum+= pay*c.get_weight(i)*c.get_weight(j);
						}
					}
				}
			}

	return Sum	;
}




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

