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
 
#ifndef _GP_CF_TEMPLATED_PRICING_H
#define _GP_CF_TEMPLATED_PRICING_H

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

#include "gpbase/numericconstant.h"


CC_BEGIN_NAMESPACE(ARM)


/// #define ARM_CF_INVSQRTPI 0.5641895835477563

////////////////////////////////////////////////////////////////////////
///
///     Generic Functorial representation of Univariate payoffs 
///
////////////////////////////////////////////////////////////////////////

struct UnivariateToDoubleFunc : public CC_NS(std,unary_function)<double,double> 
{
	virtual double operator()(double) const = 0;
	virtual int IntegrationDirection() const {return 0;};
	enum {UPWARD, DOWNWARD} Direction;
};


struct ptr_funUnivariateToDouble : UnivariateToDoubleFunc 
{
	typedef double (*funcUnivariateToDouble)(double );
	ptr_funUnivariateToDouble( funcUnivariateToDouble ptrF ) : itsPFunc(ptrF ) {}
	virtual double operator()(double x) const {return (*itsPFunc)(x); }
	
private:
	funcUnivariateToDouble itsPFunc;
};

////////////////////////////////////////////////////////////////////////
///
///     Generic Functorial representation of bivariate payoffs 
///
////////////////////////////////////////////////////////////////////////

struct BivariateToDoubleFunc : public CC_NS(std,binary_function)<double,double,double> 
{
	virtual double operator()(double,double) const = 0;
};


struct ptr_funBivariateToDouble : BivariateToDoubleFunc 
{
	typedef double (*funcBivariateToDouble)( double,double );
	ptr_funBivariateToDouble( funcBivariateToDouble ptrF ) : itsPFunc(ptrF ) {}
	virtual double operator()(double x,double y) const {return (*itsPFunc)(x,y); }

	
private:
	funcBivariateToDouble itsPFunc;
};

struct Swappedptr_funBivariateToDouble :BivariateToDoubleFunc 
{
	Swappedptr_funBivariateToDouble( BivariateToDoubleFunc* bdf) : itsPFunc(bdf) {}
	virtual double operator()(double x,double y) const {return (*itsPFunc)(y,x); }

	
private:
	BivariateToDoubleFunc * itsPFunc;
};

////////////////////////////////////////////////////////////////////////
///
///     Generic  representation of bivariate payoffs formula
///
////////////////////////////////////////////////////////////////////////

template<typename  smile,typename copula>
struct Templated_Pricing
{
	///   Generic methods (relative to the smile distribution and the copula
	static double Discriminant(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);

	static double Generic_Pricing_aux(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);

	/// used					in ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing
	/// and therefore  used		by	ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula::value
	static double Generic_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);

	static double Digital_Generic_Pricing_1(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double k1,double a2,double b2,double k2);
	static double Digital_Generic_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double k1,double a2,double b2,double k2);

	static double Digital_Generic_Pricing_First_Derivative_1(int i,const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double k1,double a2,double b2,double k2);

	static double Digital_Generic_Pricing_First_Derivative_2(int i,const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,double k1,double a2,double b2,double k2);


	/// used					in ARM_CF_PowerSpreadOption_Formula<smile,copula>::PowerSpreadOption_Pricing_With_Limits,
	/// and therefore	 used	by ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value

	static double Generic_Pricing_S2(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff, UnivariateToDoubleFunc* domain_limits);
	static double Generic_Pricing_S1(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff, UnivariateToDoubleFunc* domain_limits);

	static double Generic_Pricing_with_2Limits(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff,UnivariateToDoubleFunc*l1,UnivariateToDoubleFunc*l2);
	static double Certitude(const ArgumentList& copula0,double t,int n);

	///  New methods valid for student copula and others

	static double Generic_Pricing_aux_GeneralCopula(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);
	static double Generic_Pricing_GeneralCopula(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula0,double t,int n,BivariateToDoubleFunc* option_payoff);

	static double Generic_Pricing_S2_GeneralCopula(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff, UnivariateToDoubleFunc* domain_limits);
	static double Generic_Pricing_S1_GeneralCopula(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula_arg,double t,int n,BivariateToDoubleFunc* option_payoff, UnivariateToDoubleFunc* domain_limits);

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Generic implementation of the pricing of an instrument with a bivariate payoff
///   Do the computation in a gaussian space (space of variable with a gaussian marginal density, but
///			the codependence can be non gaussian )
///
///		Payoff : option_payoff should be a fonctor that return the price conditional to S1=S1 and S2=S2
///
///		smile : the smile of the the underlying, caracterized by the member function gaussian_to_distribution
///
///		copula : the codependence of the underlying is given by the mapping in gaussian space : gaussian_mix
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// used for SABR distribution and other distribution where an approach with an integration 
/// that start at a computable boundary is not possible.

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_aux(
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
	double x,y,ycorrelated,f1,f2;
	double Sum00=0.0;
	int i,j;
	for(i=0;i<n;i++) for(j=0;j<n;j++) Sum00 +=c.get_weight(i)*c.get_weight(j);

	for(i=0;i<n;i++){
		x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2;
		if((x<4.)&&(x>-4))
		{
			f1=smile::gaussian_to_distribution(Underlying1, x, t);
			for(j=0;j<n;j++){
				
				y=c.get_point(j)*ARM_NumericConstants::ARM_SQRT_2;
				ycorrelated=copula::gaussian_mix(copula_arg,x,y);
				if((ycorrelated<4.)&&(ycorrelated>-4.))
				{
					f2=smile::gaussian_to_distribution(Underlying2, ycorrelated, t);
					double paid=(*option_payoff)(f1,f2);
					Sum+= paid*c.get_weight(i)*c.get_weight(j);
					Sum0+=c.get_weight(i)*c.get_weight(j);
				}
			}
		}
	}
	return Sum/ARM_NumericConstants::ARM_PI;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Templated_Pricing<smile,copula>::Generic_Pricing
///
///					prepares the payoff fonctor and launch Generic_Pricing_aux that does the job
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/// for the spreadoption SABR+Gaussian
template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						BivariateToDoubleFunc* option_payoff)
{
	Swappedptr_funBivariateToDouble swp_payoff(option_payoff);
	double d00=Templated_Pricing<smile,copula>::Discriminant(Underlying1,Underlying2,copula_arg,t,n,option_payoff);
	double d12=Templated_Pricing<smile,copula>::Discriminant(Underlying2,Underlying1,copula_arg,t,n,&swp_payoff);

	if(d00>=d12) 
	{
		return Templated_Pricing<smile,copula>::Generic_Pricing_aux(Underlying1,Underlying2,copula_arg,t,n,option_payoff);
	}
	else
	{
		return Templated_Pricing<smile,copula>::Generic_Pricing_aux(Underlying2,Underlying1,copula_arg,t,n,&swp_payoff);
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Generic implementation of the pricing of an instrument with a bivariate payoff
/// in the case where the payoff is zero after or before za in the S1 dimension

///		In the S1 dimension , integrate after or before a coupure determined by za

///		za is determined by DistributionZaInverseLimit<smile,copula>::result(...)
///		the direction is upward if domain_limits->IntegrationDirection() == UnivariateToDoubleFunc::UPWARD
///		the direction is downward if domain_limits->IntegrationDirection() == UnivariateToDoubleFunc::DOWNWARD



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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_S1(
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
					a=(*domain_limits)(f2);
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



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///		Generic implementation of the pricing of an instrument with a bivariate payoff
///		in the case where the payoff is zero after or before za in the S1 dimension

///		In the S2 dimension , integrate after or before a coupure determined by za

///		za is determined by DistributionZaInverseLimit<smile,copula>::result(...)
///		the direction is upward if domain_limits->IntegrationDirection() == UnivariateToDoubleFunc::UPWARD
///		the direction is downward if domain_limits->IntegrationDirection() == UnivariateToDoubleFunc::DOWNWARD

///		Do the computation in a gaussian space (space of variable with a gaussian marginal density, but
///		the codependence can be non gaussian )

///		Payoff : option_payoff should be a fonctor that return the price conditional to S1=S1 and S2=S2
///
///		smile : the smile of the the underlying, caracterized by the member function gaussian_to_distribution
///
///		copula : the codependence of the underlying is given by the mapping in gaussian space : gaussian_mix
///
///		uses an information on when option =0 ( S2 < Ldown(S1)  )  through the function Ldown(S1)
///		through the function domain_limits(S1) and the direction domain_limits.IntegrationDirection()

///
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Generic_Pricing_S2(
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

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Discriminant(
						const ArgumentList& Underlying1,
						const ArgumentList& Underlying2,
						const ArgumentList& copula_arg,
						double t,
						int n,
						BivariateToDoubleFunc* option_payoff)
{
	double ycorrelated,f1,f2,tt;
			f1=smile::gaussian_to_distribution(Underlying1, 0, t);
			ycorrelated=copula::gaussian_mix(copula_arg,0,0);
			f2=smile::gaussian_to_distribution(Underlying2, ycorrelated, t);
			tt=(*option_payoff)(f1,f2);
			return (*option_payoff)(f1,f2);
}


////////////////////////////////////////////////////////////////////////
///
///      Certitude (probability of the domain taken in consideration for the pricing)
///
////////////////////////////////////////////////////////////////////////

template<class smile, class copula>
double Templated_Pricing<smile,copula>::Certitude(const ArgumentList& copula_arg,double t,int n)
{
	GaussHermite_Coefficients c(n);
	double Sum=0;
	double Sum0=0;
	double x,y,ycorrelated;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2	;
			y=c.get_point(j)*ARM_NumericConstants::ARM_SQRT_2	;
			ycorrelated=copula::gaussian_mix(copula_arg,x,y);
			if((x<4.)&&(x>-4.)&&(ycorrelated<4.)&&(ycorrelated>-4.))
			{
				Sum0+=c.get_weight(i)*c.get_weight(j);
			}
			
		}
	}
	return Sum0/ARM_NumericConstants::ARM_PI;
}





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

