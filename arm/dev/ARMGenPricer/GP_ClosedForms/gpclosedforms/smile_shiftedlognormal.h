/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_shiftedlognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SMILE_SHIFTEDLOGNORMAL_H
#define _GP_CF_SMILE_SHIFTEDLOGNORMAL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"



CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the Shifted Loognormal Distributions the  
///		 Underlying[0],	//INDEX
///      Underlying[1],	//SIGMA
///		 Underlying[2],	//ALPHA
///		 
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ShiftedLogNormal_Smile
{
	static double call_option(double f,double K,double tex,double sigma, double alpha);
	static double digital_call_option(double f,double K,double tex,double sigma,  double alpha);
	static double inverse_distribution(double f,double K,double tex, double sigma, double alpha);
	static double gaussian_to_distribution(double f,double K,double tex, double sigma, double alpha);
	static void   volshift_calibration( double f,  double tex, 
										double K1, double price1,
										double K2, double price2,
										/// result
										double& sigma,
										double& alpha);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x, double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x, double t);
	static double quantile(const ArgumentList& Underlying1, double x, double t);
	static double probability_density(const ArgumentList& Underlying, double x, double t);
	static double probability_distribution(const ArgumentList& Underlying, double x, double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t);
	static ArgumentList* HomotheticTransformation(const ArgumentList* arg, double positivenumber);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return Underlying[1];}	
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

