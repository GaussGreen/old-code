/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file student_copula.cpp
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

#include "gpclosedforms/student_copula.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/incompletebeta.h"
#include "gpclosedforms/student_copula.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/gamma.h"

#include "expt.h"

using namespace std;
CC_BEGIN_NAMESPACE(ARM)



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
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"StudentCopula::multivariate_density : bad argsize");
	 }
	 double copula_corr=copula[0];  // Student Correlation
	 double copula_degre=copula[1];  // Student degre
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
	 double copula_degre=copula[1];  // Student degre
	 if (x>=0)
	 {
		 return 1.-IncompleteBeta(copula_degre/2.,0.5,copula_degre/(copula_degre+x*x))/2. ;
	 }
	 else
		 {
		 return IncompleteBeta(copula_degre/2.,0.5,copula_degre/(copula_degre+x*x))/2. ;
	 }

}

double StudentCopula::marginal_left_limit(const ArgumentList& copula)
{
return -15;
}

double StudentCopula::marginal_right_limit(const ArgumentList& copula)
{
return 15;
}

double StudentCopula::gaussian_mix(const ArgumentList& copula,double x,double y)
{
	 int argsize=copula.size();
	 if ((argsize<1)||(argsize>1))
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"StudentCopula::gaussian_mix : bad argsize");
	 }
	 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"StudentCopula::gaussian_mix : you should not be there !!!");
	 return 0;

}


double StudentCopula::Q_Integral(const ArgumentList& copula ,double H,double x)
{
	double nu=copula[1];
	return H*pow(nu/(nu+x*x),1.+nu/2.)*Hypergeometric2F1(0.5,1.+nu/2.,1.5,-H*H/(x*x+nu))/(2.*ARM_NumericConstants::ARM_PI);
}


double StudentCopula::Q_Total_Integral(const ArgumentList& copula ,double x)
{
	double nu=copula[1];
	return pow(nu,nu/2.)*pow(x*x+nu,-(1.+nu)/2.)*gamma((1.+nu)/2.)/(2.*gamma(nu/2.))/ARM_NumericConstants::ARM_SQRT_PI;
}

double StudentCopula::marginal_quantile(const ArgumentList& copula,double x,double t)
{
	double nu=copula[1];
	double x1=(x<1e-11)? 1e-11: x;
	double rr;
	if (x1>=0.5)
	{
		rr= IncompleteBetaInverse(nu/2.,0.5,1.,1.-2.*x1);
		return sqrt((nu-nu*rr)/rr);
	}
	else
	{
		rr= IncompleteBetaInverse(nu/2.,0.5,1.,2.*x1-1.);
		return -sqrt((nu-nu*rr)/rr);
	}

}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
