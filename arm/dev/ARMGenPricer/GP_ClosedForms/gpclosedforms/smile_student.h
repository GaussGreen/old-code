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
 
#ifndef _GP_CF_SMILE_STUDENT_H
#define _GP_CF_SMILE_STUDENT_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h"


CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the Student Distributions the  
///      Underlying[0],	//NU
///		 Underlying[1],	//SHIFT
///		 Underlying[2],	//SIGMA
///
///		it means (x-SHIFT)/SIGMA is a student variable of rank NU
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Student_Smile
{
	static double probability_density(double nu,double shift, double sigma, double x);
	static double probability_distribution(double nu,double shift, double sigma, double x);
	static double inverse_distribution(double nu,double shift, double sigma,double y);
	static double inverse_distribution2(double nu,double shift, double sigma,double y);
	static double gaussian_to_distribution(double nu,double shift,double sigma,double x);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x, double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x, double t);
	static double quantile(const ArgumentList& Underlying1, double x, double t);
	static double probability_density(const ArgumentList& Underlying, double x, double t);
	static double probability_distribution(const ArgumentList& Underlying, double x, double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return 0.;}
};



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

