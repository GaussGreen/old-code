/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_glambda.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SMILE_GLAMBDA_H
#define _GP_CF_SMILE_GLAMBDA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h"


CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the GLAMBDA Distributions the  
///      Underlying[0],	//l1
///		 Underlying[1],	//l2
///		 Underlying[2],	//l3
///		 Underlying[3],	//l4
///		 Underlying[4],	//l5
///		
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GLambda_Smile
{
	static double sigma_implicit(double K,double l1,double l2,double l3, double l4, double l5,double l6);
	static double call_option(double K,double l1,double l2,double l3, double l4, double l5,double l6,int n);
	static double digital_call_option(double K,double l1,double l2,double l3, double l4, double l5,double l6,int n);
	static double inverse_distribution(double K,double l1,double l2,double l3, double l4, double l5,double l6);
	static double gaussian_to_distribution(double K,double l1,double l2,double l3, double l4, double l5,double l6);
	static double quantile(double K,double l1,double l2,double l3, double l4, double l5,double l6);
	static void QuantileAndAllDerivatives(
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
					double* der_l6);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x,double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x,double t);
	static double quantile(const ArgumentList& Underlying1, double x,double t);
	static double probability_density(const ArgumentList& Underlying, double x,double t);
	static double probability_distribution(const ArgumentList& Underlying, double x,double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x,double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x,double t);
	static ArgumentList* HomotheticTransformation(const ArgumentList* arg, double positivenumber);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return 0.;}
};



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

