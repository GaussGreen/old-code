/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_copula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/gaussian_copula.h"
#include "gpbase/numericconstant.h"


#include "expt.h"

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


double GaussianCopula::multivariate_density(const ArgumentList& copula,double x,double y)
{
	 int argsize=copula.size();
	 if ((argsize<1)||(argsize>1))
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"GaussianCopula::multivariate_density : bad argsize");
	 }
	 double copula_corr=copula[0];  //  Correlation
	
	 return 0;

}

double GaussianCopula::marginal_distribution(const ArgumentList& copula,double x,double t)
{
	 int argsize=copula.size();
	 if ((argsize<1)||(argsize>1))
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"GaussianCopula::marginal_distribution : bad argsize");
	 }
	 double copula_corr=copula[0];  //  Correlation

	 
	return 0;
	 

}

double GaussianCopula::marginal_left_limit(const ArgumentList& copula)
{
return -10;
}

double GaussianCopula::marginal_right_limit(const ArgumentList& copula)
{
return 10;
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
