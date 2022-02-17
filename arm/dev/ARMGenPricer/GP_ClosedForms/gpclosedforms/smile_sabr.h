/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_sabr.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SMILE_SABR_H
#define _GP_CF_SMILE_SABR_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h"


CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the SABR Distributions the  
///      Underlying[0],	//INDEX1
///		 Underlying[1],	//ALPHA1
///		 Underlying[2],	//BETA1
///		 Underlying[3],	//RHO1
///		 Underlying[4],	//NU1
///		 Underlying[5]	//FLAG
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SABR_smile
{
	static double sigma_implicit(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,double alpha_exp,double alpha_tanh,double kb_tanh);
	static double call_option(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,double alpha_exp,double alpha_tanh,double kb_tanh);
	static double digital_call_option(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,double alpha_exp,double alpha_tanh,double kb_tanh);
	static double inverse_distribution(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,double alpha_exp,double alpha_tanh,double kb_tanh);
	static double gaussian_to_distribution(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,double alpha_exp,double alpha_tanh,double kb_tanh);

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
	static double LowerBoundary(const ArgumentList& Underlying) {return 0.;}
};



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

