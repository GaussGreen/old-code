/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_student.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include <glob/expt.h>
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/smile_student.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/incompletebeta.h"

#include "gpclosedforms/smile_student.h"



CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Class : Student_Smile
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline double cmin(double x, double y) {return ((x<=y) ? x : y);}
inline double cmax(double x, double y) {return ((x<=y) ? y : x);}


double Student_Smile::probability_density(double nu,double shift, double sigma,double x)
{
	double z=(x-shift)/sigma;
	return gamma((nu+1.0)/2.0)/(gamma(nu/2.0)*ARM_NumericConstants::ARM_SQRT_PI*sqrt(nu))*pow((1.0+z*z/nu),-(nu+1.0)/2.0)/sigma;
}

double Student_Smile::probability_distribution(double nu,double shift, double sigma,double x)
{
	double z=(x-shift)/sigma;
	if (x>=0)
	 {
		 return 1.-IncompleteBeta(nu/2.,0.5, nu/(nu+z*z))/2. ;
	 }
	 else
		 {
		 return IncompleteBeta(nu/2.,0.5,nu/(nu+z*z))/2. ;
	 }
}



double Student_Smile::inverse_distribution(double nu,double shift, double sigma,double x)
{
	double sign=(x>0.5)?1.0 : -1.0;

	double inversebeta= sign*sqrt(nu*(-1.+1./IncompleteBetaInverse(nu/2.,0.5,1.,1.-2.*x)));

	return shift+sigma*inversebeta;
}

double Student_Smile::inverse_distribution2(double nu,double shift, double sigma,double K)
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double nu0;
		double shift0;
		double sigma0;
		DistributionToInverse(double nu1,double shift1, double sigma1):
		nu0(nu1),shift0(shift1),sigma0(sigma1) {}
		
		virtual double operator() (double K0)  const
		{
			if(K0<=1e-9)
			{
				return 0.;
			}
				
			return Student_Smile::probability_distribution(nu0,shift0,sigma0,K0);
		}
	};
	double s=cmax(fabs(shift),0.1);
	DistributionToInverse x(nu,shift,sigma);
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(K,s,s/10.,1e-12);  // suppose that  Y is always postive
}

double Student_Smile::gaussian_to_distribution(double nu,double shift, double sigma,double x)
{
	return inverse_distribution(nu,shift,sigma,ARM_GaussianAnalytics::cdfNormal(x));
}


////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a generic implementation that uses the following conventions :
/// 
///				For the SABR Distributions the  
///      Underlying[0],	//NU
///		 Underlying[1],	//SHIFT
///		 Underlying[2],	//SIGMA
///
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////

double Student_Smile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::gaussian_to_distribution  : bad argsize");
	}

		return inverse_distribution(
			Underlying[0],	//nu
			Underlying[1],	//shift
			Underlying[2],	//sigma
			ARM_GaussianAnalytics::cdfNormal(x)
			);

	
}
double Student_Smile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::gaussian_to_distribution  : bad argsize");
	}	
	return	ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying[0],
			Underlying[1],
			Underlying[2],
			x));
}

double Student_Smile::distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	double s=0.0001;
	ArgumentList* Underlying1;
	ArgumentList* Underlying2;
	switch(i)
	{	
	case 0 :										/// NU
		{
			Underlying1=new ArgumentList(Underlying[0]-s,Underlying[1],Underlying[2]);
			Underlying2=new ArgumentList(Underlying[0]+s,Underlying[1],Underlying[2]);
			break;
		}
	case 1 :										/// SHIFT 
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1]-s,Underlying[2]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1]+s,Underlying[2]);
			break;
		}
	case 2 :										/// SIGMA
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]-s);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]+s);
			break;
		}
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::distribution_to_gaussian_first_derivative : incorrect input");
	};
	double v1=distribution_to_gaussian( *Underlying1,  x,  t);
	double v2=distribution_to_gaussian( *Underlying2,  x,  t);
	delete Underlying1;
	delete Underlying2;
	return (v2-v1)/(2.*s);
				
}



double Student_Smile::quantile(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::quantile  : bad argsize");
	}
	return inverse_distribution(
		Underlying[0],	//nu
		Underlying[1],	//shift
		Underlying[2],	//sigma
		x
		);		
}

double Student_Smile::probability_density(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::quantile  : bad argsize");
	}
		return probability_density(
		Underlying[0],	//nu
		Underlying[1],	//shift
		Underlying[2],	//sigma
		x
		);		
}

double Student_Smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::quantile  : bad argsize");
	}
	ArgumentList arg(
		Underlying[0],	//nu
		Underlying[1],	//shift
		Underlying[2],	//sigma
		x
		);	

	return probability_distribution(
		Underlying[0],	//nu
		Underlying[1],	//shift
		Underlying[2],	//sigma
		x
		);		
}

double Student_Smile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	{
	double s=0.0001;
	ArgumentList* Underlying1;
	ArgumentList* Underlying2;
	switch(i)
	{	
	case 0 :										/// NU
		{
			Underlying1=new ArgumentList(Underlying[0]-s,Underlying[1],Underlying[2]);
			Underlying2=new ArgumentList(Underlying[0]+s,Underlying[1],Underlying[2]);
			break;
		}
	case 1 :										/// SHIFT 
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1]-s,Underlying[2]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1]+s,Underlying[2]);
			break;
		}
	case 2 :										/// SIGMA
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]-s);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]+s);
			break;
		}
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Student_Smile::distribution_to_gaussian_first_derivative : incorrect input");
	};
	double v1=probability_distribution( *Underlying1,  x,  t);
	double v2=probability_distribution( *Underlying2,  x,  t);
	delete Underlying1;
	delete Underlying2;
	return (v2-v1)/(2.*s);
				
}

}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/