/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_bisabr.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date October 2006
 */
 
#ifndef _GP_CF_BISMILE_SABR_H
#define _GP_CF_BISMILE_SABR_H

#include "firsttoinc.h"
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
///      Underlying[5],	//INDEX2
///		 Underlying[6],	//ALPHA2
///		 Underlying[7],	//BETA2
///		 Underlying[8],	//RHO2
///		 Underlying[9],	//NU2
///		 Underlying[10]	//RHOS
///		 Underlying[11]	//RHOV
///		 Underlying[12]	//RHOC12
///		 Underlying[13]	//RHOC21
///		 Underlying[14]	//T
///		 Underlying[15]	//FLAG
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BiSABR_smile
{
	static double sigma_implicit(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,double alpha_exp,double alpha_tanh,double kb_tanh);
	static 	double BiSABR_smile::call_option(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double K,double T,int flag,double alpha_exp,double alpha_tanh,double kb_tanh);
	static double BiSABR_smile::digital_call_option(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double K,double T,int flag,double alpha_exp,double alpha_tanh,double kb_tanh);
	static double BiSABR_smile::inverse_distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double proba,double T,int flag,double alpha_exp,double alpha_tanh,double kb_tanh);
	static  double BiSABR_smile::gaussian_to_distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag,double alpha_exp,double alpha_tanh,double kb_tanh);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x, double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x, double t);
	static double quantile(const ArgumentList& Underlying1, double x, double t);
	static double probability_density(const ArgumentList& Underlying, double x, double t);
	static double probability_distribution(const ArgumentList& Underlying, double x, double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t);
	static ArgumentList* HomotheticTransformation(const ArgumentList* arg, double positivenumber);
	static ArgumentList* HomotheticTransformation2(const ArgumentList* arg, double positivenumber);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return 0.;}
};



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

