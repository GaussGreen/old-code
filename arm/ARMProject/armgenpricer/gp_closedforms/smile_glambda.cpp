/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_glambda.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2005
 */

#include "firsttoinc.h"

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include "expt.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/smile_glambda.h"



CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001

inline double cmax(double x, double y) {return ((x<=y) ? y : x);}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Class : GLambda_Smile
///
///    Approach following the book 
///			Fitting statistical distributions
///			THe Generalized Lambda Distribution and Generalized Bootstrap Methods
///			by Z. A. Karian and E. J. Dudewicz

///		and
///
///		Statistical Modeling with Quantile Functions
///		by W. G. Gilchrist
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GLambda_Smile::quantile(double p,double l1,double l2,double l3, double l4, double l5,double l6)
{
	return l1*(1.-pow(1.-p,l2))+l3*pow(p,l4)+l5*p/pow(1.-p,l6);

}

double GLambda_Smile::quantile(const ArgumentList& Underlying, double x,double t)
{
	int argsize=Underlying.size();
	if ((argsize<6)||(argsize>6))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"GLambda_Smile::quantile  : bad argsize");
	}
	return quantile(
		x,
		Underlying[0],	//l1
		Underlying[1],	//l2
		Underlying[2],	//l3
		Underlying[3],	//l4
		Underlying[4],	//l5
		Underlying[5]	//l6
		);		
}

void GLambda_Smile::QuantileAndAllDerivatives(
					double p,
					double l1,
					double l2,
					double l3,
					double l4,
					double l5,
					double l6,
					double* digitalprice,
					double* der_l1,
					double* der_l2,
					double* der_l3,
					double* der_l4,
					double* der_l5,
					double* der_l6)
{
	double pow1moinsp_p6=pow(1.0-p,-l6);
	double pow1moinsp_p2=pow(1.0-p,l2);
	double log1moinsp=log(1.0-p);
	double powp_p4=pow(p,l4);
	*digitalprice=quantile(p,l1,l2,l3,l4,l5,l6);

	*der_l1=1.-pow1moinsp_p2;

	*der_l2=-l1*pow1moinsp_p2*log1moinsp;

	*der_l3=powp_p4;

	*der_l4=l3*powp_p4*log(p);

	*der_l5=p*pow1moinsp_p6;

	*der_l6=-l5*log1moinsp*p*pow1moinsp_p6;
	return;

}


double GLambda_Smile::call_option(double K,double l1,double l2,double l3, double l4, double l5,double l6,int n)
{
	GaussLegendre_Coefficients c(n,0.0,1.0);
	double Sum=0;
	int i;
	for(i=0;i<n;i++)
	{
			Sum+=cmax(quantile(c.get_point(i),l1,l2,l3,l4,l5,l6)-K,0);
	}
	return Sum;
}


double GLambda_Smile::digital_call_option(double K,double l1,double l2,double l3, double l4, double l5,double l6,int n)
{
	double shift=K/10000.;
	double result= (-GLambda_Smile::call_option(K+shift,l1,l2,l3,l4,l5,l6,n)+
		GLambda_Smile::call_option(K-shift,l1,l2,l3,l4,l5,l6,n))/(2.*shift);
	if (result<=0) return 1e-12;
	else return result;
}



double GLambda_Smile::inverse_distribution(double K,double l1,double l2,double l3, double l4, double l5,double l6)
{
	
	return quantile(K,l1,l2,l3,l4,l5,l6);  // suppose that  Y is always postive
}


double GLambda_Smile::gaussian_to_distribution(double x,double l1,double l2,double l3, double l4, double l5,double l6)
{
	return inverse_distribution(ARM_GaussianAnalytics::cdfNormal(x),l1,l2,l3,l4,l5,l6);
}


////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a generic implementation that uses the following conventions :
/// 
///				For the SABR Distributions the  
///      Underlying[0],	//l1
///		 Underlying[1],	//l2
///		 Underlying[2],	//l3
///		 Underlying[3],	//l4
///		 Underlying[4],	//l5
///		 Underlying[5]	//l6
///		 
///
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////

double GLambda_Smile::gaussian_to_distribution(const ArgumentList& Underlying, double x,double t)
{
	int argsize=Underlying.size();
	if ((argsize<6)||(argsize>6))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"GLambda_Smile::gaussian_to_distribution  : bad argsize");
	}

		return inverse_distribution(
			ARM_GaussianAnalytics::cdfNormal(x),
			Underlying[0],	//l1
			Underlying[1],	//l2
			Underlying[2],	//l3
			Underlying[3],	//l4
			Underlying[4],	//l5
			Underlying[5]	//l6
			);
}
double GLambda_Smile::distribution_to_gaussian(const ArgumentList& Underlying, double x,double t)
{
	int argsize=Underlying.size();
	if ((argsize<6)||(argsize>6))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"GLambda_Smile::gaussian_to_distribution  : bad argsize");
	}	
		return ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			x,t
			));
}

double GLambda_Smile::distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x,double t)
{
	double s=0.00001;
	ArgumentList* Underlying1;
	ArgumentList* Underlying2;
	switch(i)
	{	
	case 0 :										/// l1
		{
			Underlying1=new ArgumentList(Underlying[0]-s,Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]);
			Underlying2=new ArgumentList(Underlying[0]+s,Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]);
			break;
		}
	case 1 :											/// l2 
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1]-s,Underlying[2],Underlying[3],Underlying[4],Underlying[5]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1]+s,Underlying[2],Underlying[3],Underlying[4],Underlying[5]);
			break;
		}
	case 2 :											/// l3
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]-s,Underlying[3],Underlying[4],Underlying[5]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]+s,Underlying[3],Underlying[4],Underlying[5]);
			break;
		}
	case 3 :											/// l4
		{		
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3]-s,Underlying[4],Underlying[5]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3]+s,Underlying[4],Underlying[5]);
			break;
		}
	case 4 :											/// l5
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4]-s,Underlying[5]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4]+s,Underlying[5]);
			break;
		}
	case 5 :											/// l6
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]-s);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]+s);
			break;
		}
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::distribution_to_gaussian_first_derivative : incorrect input");
	};
	double v1=distribution_to_gaussian( *Underlying1,  x,t);
	double v2=distribution_to_gaussian( *Underlying2,  x,t);
	delete Underlying1;
	delete Underlying2;
	return (v2-v1)/(2.*s);
				
}

ArgumentList*  GLambda_Smile:: HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	return new ArgumentList(
		positivenumber*(*Underlying)[0],
		(*Underlying)[1],
		positivenumber*(*Underlying)[2],
		(*Underlying)[3],
		positivenumber*(*Underlying)[4],
		(*Underlying)[5]
		);
}




double GLambda_Smile::probability_distribution(const ArgumentList& Underlying, double K,double t)
{
	if (K==0) return 0;
	struct QuantileToInverse : public DoubleToDoubleFunc 
	{
		QuantileToInverse(double l0a,double l1a, double l2a, double l3a, double l4a, double l5a):
		l0(l0a),l1(l1a),l2(l2a),l3(l3a),l4(l4a),l5(l5a) {}
		virtual double operator() (double K0)  const
		{
			return GLambda_Smile::quantile(K0,l0,l1,l2,l3,l4,l5);
		}
		double l0;
		double l1;
		double l2;
		double l3;
		double l4;
		double l5;
	};
	
	QuantileToInverse x(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]);
	return Inverse(x,Inverse::BOUNDEDBY0AND1)(K,0.5,0.05,1e-12);
}

double GLambda_Smile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x,double t)
{
	return 0.0;  /// not implemented
}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/