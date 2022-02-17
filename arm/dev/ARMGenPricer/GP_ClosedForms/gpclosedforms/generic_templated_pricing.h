/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file generictemplated_pricing.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_GENERICTEMPLATED_PRICING_H
#define _GP_CF_GENERICTEMPLATED_PRICING_H

#include <glob/firsttoinc.h>
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
#include "gpclosedforms/ZaInverseLimit.h"
#include "gpclosedforms/templated_pricing.h"

#include "gpbase/numericconstant.h"


CC_BEGIN_NAMESPACE(ARM)


#define ARM_CF_INVSQRTPI 0.5641895835477563



///////////////////////////////////////////////////////////////////////////
///
///
///
///
/// for student and other copulas
///
///				Generic_Pricing_aux_GeneralCopula does the job
///
///////////////////////////////////////////////////////////////////////////

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_aux_GeneralCopula(
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
		x3=smile::quantile(Underlying1,x2,t);
		for(j=0;j<n;j++){
			y=c.get_point(j);
			y2=copula::marginal_distribution(copula_arg,y,t);
			y3=smile::quantile(Underlying2,y2, t);
			double dens=copula::multivariate_density(copula_arg,x,y);;
			Sum+= (*option_payoff)(x3,y3)*dens*c.get_weight(i)*c.get_weight(j);
			Sum0+=c.get_weight(i)*c.get_weight(j);	
		}	
	}
	return Sum;
}

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_GeneralCopula(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff)
{
	Swappedptr_funBivariateToDouble swp_payoff(option_payoff);		/// if necessary
	
	return Templated_Pricing<smile,copula>::Generic_Pricing_aux_GeneralCopula(Underlying1,Underlying2,copula_arg,t,n,option_payoff);
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Generic implementation S2 of the pricing of an instrument with a bivariate payoff
///   Do the computation in a gaussian space (space of variable with a gaussian marginal density, but
///			the codependence can be non gaussian )
///
///		Payoff : option_payoff should be a fonctor that return the price conditional to S1=S1 and S2=S2
///
///		smile : the smile of the the underlying, caracterized by the member function gaussian_to_distribution
///
///		copula : the codependence of the underlying is given by the mapping in gaussian space : gaussian_mix
///
///		uses an information on when option =0   
///		through the function domain_limits(S1) and the direction domain_limits.IntegrationDirection()
///
///
///
///				FOR   GAUSSIAN   COPULA   ONLY    !!!!!!!!!!!!!
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_S2_GeneralCopula(
								const ArgumentList& Underlying1,
								const ArgumentList& Underlying2,
								const ArgumentList& copula_arg,
								double t,
								int n,
								BivariateToDoubleFunc* option_payoff,
								UnivariateToDoubleFunc* domain_limits)
{
	GaussHermite_Coefficients c(n);
	GaussLaguerre_Coefficients c2(0.,n);
	double Sum=0;
	double x,y,ycorrelated,f1,f2,a;
	switch (domain_limits->IntegrationDirection())
	{
	case UnivariateToDoubleFunc::UPWARD :
		{
			for(int i=0;i<n;i++)
			{
				x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2	;
				if((x<4.)&&(x>-4))
				{
					f1=smile::gaussian_to_distribution(Underlying1, x, t);
					a=(*domain_limits)(f1);
					double za=0;
					za=DistributionZaInverseLimit<smile,copula>::result(Underlying2,copula_arg,x,t,a);
					// za soulution de smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,za), t)==a
					// ici sera donc solution de student to distribution (za...)=a
					// on va integrer a partir de za
					for(int j=0;j<n;j++)
					{
						y=za+c2.get_point(j);
						ycorrelated=copula::gaussian_mix(copula_arg,x,y);
						if((ycorrelated<4.)&&(ycorrelated>-4.))
						{
							f2=smile::gaussian_to_distribution(Underlying2, ycorrelated, t);
							Sum+= (*option_payoff)(f1,f2)*exp(-y*y/2+c2.get_point(j))*c.get_weight(i)*c2.get_weight(j);
						}
					}
				}
			}
			break;
		}
		
	case UnivariateToDoubleFunc::DOWNWARD :
		{
			for(int i=0;i<n;i++)
			{
				x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2	;
				if((x<4.)&&(x>-4))
				{
					f1=smile::gaussian_to_distribution(Underlying1, x, t);
					a=(*domain_limits)(f1);
					double za=0;
					za=DistributionZaInverseLimit<smile,copula>::result(Underlying2,copula_arg,x,t,a);
					// za soulution de smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,za), t)==a
					for(int j=0;j<n;j++)
					{		
						y=za-c2.get_point(j); 
						ycorrelated=copula::gaussian_mix(copula_arg,x,y);
						if((ycorrelated<4.)&&(ycorrelated>-4.))
						{
							f2=smile::gaussian_to_distribution(Underlying2, ycorrelated, t);
							Sum+= (*option_payoff)(f1,f2)*exp(-y*y/2+c2.get_point(j))*c.get_weight(i)*c2.get_weight(j);
						}
					}
				}
			}
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Generic_Pricing_S2 : domain_limits->IntegrationDirection() : bad value ");
			break;
		}
	}
	return Sum/(ARM_NumericConstants::ARM_SQRT_2*ARM_NumericConstants::ARM_PI)	;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Generic implementation S1 of the pricing of an instrument with a bivariate payoff
///   Do the computation in a gaussian space (space of variable with a gaussian marginal density, but
///			the codependence can be non gaussian )
///
///		Payoff : option_payoff should be a fonctor that return the price conditional to S1=S1 and S2=S2
///
///		smile : the smile of the the underlying, caracterized by the member function gaussian_to_distribution
///
///		copula : the codependence of the underlying is given by the mapping in gaussian space : gaussian_mix
///
///		uses an information on when option =0 ( S2 < Ldown(S1)  )  
///		through the function domain_limits(S1) and the direction domain_limits.IntegrationDirection()
///
///
///				FOR   GAUSSIAN   COPUL   ONLY    !!!!!!!!!!!!!
///
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_S1_GeneralCopula(
								const ArgumentList& Underlying1,
								const ArgumentList& Underlying2,
								const ArgumentList& copula_arg,
								double t,
								int n,
								BivariateToDoubleFunc* option_payoff,
								UnivariateToDoubleFunc* domain_limits)
{
	GaussHermite_Coefficients c(n);
	GaussLaguerre_Coefficients c2(0.,n);
	double Sum=0;
	double x,y,xcorrelated,f1,f2,a,za;
	switch (domain_limits->IntegrationDirection())
	{
	case UnivariateToDoubleFunc::UPWARD :
		{
			for(int i=0;i<n;i++)
			{
				y=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2	;
				if((y<4.)&&(y>-4))
				{
					f2=smile::gaussian_to_distribution(Underlying2, y, t);
					a=(*domain_limits)(f2);  // retourne la limite a partir de la quelle on integre
					// za soulution de smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,za), t)==a
					za=DistributionZaInverseLimit<smile,copula>::result(Underlying2,copula_arg,y,t,a);
					for(int j=0;j<n;j++)
					{
						x=za+c2.get_point(j);
						xcorrelated=copula::gaussian_mix(copula_arg,y,x);
						if((xcorrelated<4.)&&(xcorrelated>-4.))
						{
							f1=smile::gaussian_to_distribution(Underlying1, xcorrelated, t);
							Sum+= (*option_payoff)(f1,f2)*exp(-x*x/2+c2.get_point(j))*c.get_weight(i)*c2.get_weight(j);
						}
					}
				}
			}
			break;
		}
		
	case UnivariateToDoubleFunc::DOWNWARD :
		{
			for(int i=0;i<n;i++)
			{
				y=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2	;
				if((y<4.)&&(y>-4))
				{
					f2=smile::gaussian_to_distribution(Underlying2, y, t);
					a=(*domain_limits)(f2);
					double za=0;
					
					za=DistributionZaInverseLimit<smile,copula>::result(Underlying2,copula_arg,y,t,a);	
					// za soulution de smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,za), t)==a
					for(int j=0;j<n;j++)
					{		
						x=za-c2.get_point(j); 
						xcorrelated=copula::gaussian_mix(copula_arg,x,y);
						if((xcorrelated<4.)&&(xcorrelated>-4.))
						{
							f1=smile::gaussian_to_distribution(Underlying1, xcorrelated, t);
							Sum+= (*option_payoff)(f1,f2)*exp(-x*x/2+c2.get_point(j))*c.get_weight(i)*c2.get_weight(j);
						}
					}
				}
			}
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Generic_Pricing_S2 : domain_limits->IntegrationDirection() : bad value ");
			break;
		}
	}
	return Sum/(ARM_NumericConstants::ARM_SQRT_2*ARM_NumericConstants::ARM_PI)	;
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

