/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file ZaInverseLimit.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_ZAINVERSELIMIT_H
#define _GP_CF_ZAINVERSELIMIT_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gaussian_copula.h"

#include <cmath>

#include "gpclosedforms/inverse.h" /// for the DoubleToDoubleFunc

#include "gpclosedforms/smile_2logsmiled.h"
#include "gpclosedforms/smile_bisabr.h"
#include "gpclosedforms/smile_glambda.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/smile_shiftedlognormal.h"



CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///			Computes  za soulution de smile::gaussian_to_distribution(Underlying2, copula::gaussian_mix(copula_arg,x,za), t)==a 
///				it is necessary for spliting the integration at singularities
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class smile, class copula>
struct DistributionZaInverseLimit
{
	static double result (const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double za;
		struct PricingFunctionToInverse : public DoubleToDoubleFunc 
		{
			const ArgumentList& copula_arg0;
			const ArgumentList& Underlying20;
			double x0;
			double t0;
			PricingFunctionToInverse(const ArgumentList& Underlying200,const ArgumentList& copula_arg00,double x00,double t00):
			Underlying20(Underlying200),copula_arg0(copula_arg00),x0(x00),t0(t00)
			{}
			virtual double operator() (double za0) const 
			{
				return smile::gaussian_to_distribution(Underlying20, copula::gaussian_mix(copula_arg0,x0,za0), t0);
			}
		};
		PricingFunctionToInverse func(Underlying2,copula_arg,x,t);
		za=Inverse(func,Inverse::REAL)(a,0.0,0.1,1e-12);
		return za;
	}
};

template<class smile1,class smile2, class copula>
struct DistributionZaInverseLimit_NonSymmetric
{
	static double result (const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double za;
		struct PricingFunctionToInverse : public DoubleToDoubleFunc 
		{
			const ArgumentList& copula_arg0;
			const ArgumentList& Underlying20;
			double x0;
			double t0;
			PricingFunctionToInverse(const ArgumentList& Underlying200,const ArgumentList& copula_arg00,double x00,double t00):
			Underlying20(Underlying200),copula_arg0(copula_arg00),x0(x00),t0(t00)
			{}
			virtual double operator() (double za0) const 
			{
				return smile2::gaussian_to_distribution(Underlying20, copula::gaussian_mix(copula_arg0,x0,za0), t0);
			}
		};
		PricingFunctionToInverse func(Underlying2,copula_arg,x,t);
		za=Inverse(func,Inverse::REAL)(a,0.0,0.1,1e-12);
		return za;
	}
};

// FIXMEFRED: mig.vc8 (22/05/2007 18:46:40):template<>
template<>
struct DistributionZaInverseLimit<ShiftedLogNormal_Smile,GaussianCopula>
{
	static double result(const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double f=Underlying2[0];
		double sigma=Underlying2[1];
		double corr=copula_arg[0];
		return (sigma*sigma*t-2*corr*sigma*sqrt(t)*x+2*log(a/f))/(2*sqrt(1-corr*corr)*sigma*sqrt(t));
	}
};
/// la non existence des templates partiels nous oblige a creer autant de structure que de type de smile
/// ceci devrait cesser dans une version plu recente du compilateur
// FIXMEFRED: mig.vc8 (22/05/2007 18:46:40):template<>
template<>
struct DistributionZaInverseLimit_NonSymmetric<ShiftedLogNormal_Smile,ShiftedLogNormal_Smile,GaussianCopula>
{
	static double result(const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double f=Underlying2[0];
		double sigma=Underlying2[1];
		double corr=copula_arg[0];
		return (sigma*sigma*t-2*corr*sigma*sqrt(t)*x+2*log(a/f))/(2*sqrt(1-corr*corr)*sigma*sqrt(t));
	}
};
// FIXMEFRED: mig.vc8 (22/05/2007 18:46:40):template<>
template<>
struct DistributionZaInverseLimit_NonSymmetric<SABR_smile,ShiftedLogNormal_Smile,GaussianCopula>
{
	static double result(const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double f=Underlying2[0];
		double sigma=Underlying2[1];
		double corr=copula_arg[0];
		return (sigma*sigma*t-2*corr*sigma*sqrt(t)*x+2*log(a/f))/(2*sqrt(1-corr*corr)*sigma*sqrt(t));
	}
};
// FIXMEFRED: mig.vc8 (22/05/2007 18:46:40):template<>
template<>
struct DistributionZaInverseLimit_NonSymmetric<GLambda_Smile,ShiftedLogNormal_Smile,GaussianCopula>
{
	static double result(const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double f=Underlying2[0];
		double sigma=Underlying2[1];
		double corr=copula_arg[0];
		return (sigma*sigma*t-2*corr*sigma*sqrt(t)*x+2*log(a/f))/(2*sqrt(1-corr*corr)*sigma*sqrt(t));
	}
};

// FIXMEFRED: mig.vc8 (22/05/2007 18:46:40):template<>
template<>
struct DistributionZaInverseLimit_NonSymmetric<BiSABR_smile,ShiftedLogNormal_Smile,GaussianCopula>
{
	static double result(const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double f=Underlying2[0];
		double sigma=Underlying2[1];
		double corr=copula_arg[0];
		return (sigma*sigma*t-2*corr*sigma*sqrt(t)*x+2*log(a/f))/(2*sqrt(1-corr*corr)*sigma*sqrt(t));
	}
};

// FIXMEFRED: mig.vc8 (22/05/2007 18:46:40):template<>
template<>
struct DistributionZaInverseLimit_NonSymmetric<Spread_Shifted2LogNormal_Smile,ShiftedLogNormal_Smile,GaussianCopula>
{
	static double result(const ArgumentList& Underlying2,const ArgumentList& copula_arg,double x,double t,double a)
	{
		double f=Underlying2[0];
		double sigma=Underlying2[1];
		double corr=copula_arg[0];
		return (sigma*sigma*t-2*corr*sigma*sqrt(t)*x+2*log(a/f))/(2*sqrt(1-corr*corr)*sigma*sqrt(t));
	}
};

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

