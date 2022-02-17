/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file student_copula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_STUDENTCOPULADESQCRIPTION_H
#define _GP_CF_STUDENTCOPULADESQCRIPTION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/generic_copula.h"

CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////////////////////
///
///	Class  : Student Copula : two args  : correlation and degré
/// 
///////////////////////////////////////////////////////////////////////////////////////

struct StudentCopula : public GenericCopula
{
	/// this the generic interface to copula 
	static  double gaussian_mix(const ArgumentList& copula,double x,double y);
	static  double multivariate_density(const ArgumentList& copula,double x,double y);
	static  double marginal_distribution(const ArgumentList& copula,double x,double t);
	static  double marginal_quantile(const ArgumentList& copula,double x,double t);
	static  double marginal_left_limit(const ArgumentList& copula);
	static  double marginal_right_limit(const ArgumentList& copula);

	///   for digital computations in the case of elliptic multivariate distribution
	/// integral up to H ofthe elliptic multivariate_density with 0 correlations with respect to S1
	static	double  Q_Integral(const ArgumentList& copula , double H, double z);
	/// limit of Q_Integral(H, z)  when H tend toward infinity
	static  double  Q_Total_Integral(const ArgumentList& copula ,  double z);

	StudentCopula(ArgumentList* structure0):GenericCopula(structure0) {}
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

