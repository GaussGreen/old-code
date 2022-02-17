/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file copula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_COPULA_H
#define _GP_CF_COPULA_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h"

CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////////////////////
///
///	Class  : Gaussian Copula : one arg  : correlation
/// 
///////////////////////////////////////////////////////////////////////////////////////

struct GaussianCopula
{
	/// this the generic interface to copula 
	static double gaussian_mix(const ArgumentList& copula,double x,double y);
};

///////////////////////////////////////////////////////////////////////////////////////
///
///	Class  : Student Copula : two args  : correlation and degré
/// 
///////////////////////////////////////////////////////////////////////////////////////

struct StudentCopula
{
	/// this the generic interface to copula 
	static double gaussian_mix(const ArgumentList& copula,double x,double y);
	static double multivariate_density(const ArgumentList& copula,double x,double y);
	static double marginal_distribution(const ArgumentList& copula,double x,double t);
	static double marginal_left_limit(const ArgumentList& copula);
	static double marginal_right_limit(const ArgumentList& copula);
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

