/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file copula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/gaussian_copula.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/incompletebeta.h"

#include <glob/expt.h>

using namespace std;
CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////////////////////
///
/// Compute the gaussian mix , that mean the resulting point in the gaussian space
/// this means that if X ~ N[0,1] and Y ~ N[0,1] the resulting point Z will be such that
/// (X,Z) will have the searched codependance , stigmatized by the copula.
/// here the copula is the the gaussian copula


///			Contens of the argumentlist for the gaussian copula:
///
///		Underlying[0],	//Correlation
///
///////////////////////////////////////////////////////////////////////////////////////
double GaussianCopula::gaussian_mix(const ArgumentList& copula,double x,double y)
{
	 int argsize=copula.size();
	 if ((argsize<1)||(argsize>1))
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"GaussianCopula::gaussian_mix : bad argsize");
	 }
	 double copula_corr=copula[0];  // Gaussian Correlation
	 return copula_corr*x+sqrt(1-copula_corr*copula_corr)*y;

}


///////////////////////////////////////////////////////////////////////////////////////
///
/// 
/// here the copula is the the student copula

/// take gaussian points and add the copula dependency to return the density
///			Contens of the argumentlist for the gaussian copula:
///
///		Underlying[0],	//Correlation
///		Underlying[1],	//Degre
///
///////////////////////////////////////////////////////////////////////////////////////
double StudentCopula::multivariate_density(const ArgumentList& copula,double x,double y)
{
	 int argsize=copula.size();
	 if ((argsize<2)||(argsize>2))
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"StudentCopula::marginal_distribution : bad argsize");
	 }
	 double copula_corr=copula[0];  // Student Correlation
	 double copula_degre=copula[0];  // Student degre
	 return pow(1.+(x*x+y*y-2.*x*y*copula_corr)/(copula_degre*(1.-copula_corr*copula_corr)),-(copula_degre+2.)/2.)/
		 (sqrt(1-copula_corr*copula_corr)* ARM_NumericConstants::ARM_PI*2.);

}

double StudentCopula::marginal_distribution(const ArgumentList& copula,double x,double t)
{
	 int argsize=copula.size();
	 if ((argsize<2)||(argsize>2))
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"StudentCopula::marginal_distribution : bad argsize");
	 }
	 double copula_corr=copula[0];  // Student Correlation
	 double copula_degre=copula[0];  // Student degre
	 if (x>=0)
	 {
		 return 1.-IncompleteBeta(copula_degre/2.,0.5,copula_degre*t/(copula_degre*t+x*x))/2. ;
	 }
	 else
		 {
		 return 1.-IncompleteBeta(copula_degre/2.,0.5,copula_degre*t/(copula_degre*t+x*x))/2. ;
	 }

}

double StudentCopula::marginal_left_limit(const ArgumentList& copula)
{
return -10;
}

double StudentCopula::marginal_right_limit(const ArgumentList& copula)
{
return -10;
}

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
