/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_2logsmiled.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */
 
#ifndef _GP_CF_SMILE_2LOGSMILED_H
#define _GP_CF_SMILE_2LOGSMILED_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"



CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the Shifted Lognormal Distributions the  
///		 Underlying[0],	//INDEX1
///      Underlying[1],	//SIGMA1
///		 Underlying[2],	//INDEX2
///      Underlying[3],	//SIGMA2
///		 Underlying[4], //ALPHA
///		 Underlying[5],	//RHO
///		 Underlying[6],	//N
///		 
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Spread_Shifted2LogNormal_Smile
{
	static double call_option(double f1,double f2,double K,double tex,double sigma1,
		double sigma2,double alpha, double rho, int n);
	static double digital_call_option(double f1,double f2,double K,double tex,double sigma1, 
		double sigma2,double alpha,double rho, int n);
	static double inverse_distribution(double f1,double f2,double K,double tex,double sigma1,
		double sigma2, double alpha,double rho, int n);
	static double gaussian_to_distribution(double f1,double f2,double K,double tex,double sigma1, 
		double sigma2,double alpha,double rho, int n);


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

