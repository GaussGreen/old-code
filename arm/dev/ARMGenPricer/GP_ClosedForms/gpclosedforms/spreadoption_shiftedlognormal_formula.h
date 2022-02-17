/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_shiftedlognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_SHIFTEDLOGNORMAL_FORMULA_H
#define _GP_CF_SPREADOPTION_SHIFTEDLOGNORMAL_FORMULA_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result

#include "powerspreadoption.h"		/// for instanciation of the templates
#include "smile_shiftedlognormal.h"
#include "gaussian_copula.h"



CC_BEGIN_NAMESPACE(ARM)


struct ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula
{
	typedef ARM_CF_PowerSpreadOption_Formula<ShiftedLogNormal_Smile,GaussianCopula> structure;
	enum ArgumentType
	{
		INDEX1,			//i=0
		INDEX2,			//i=1
		SIGMA1,			//i=2	
		ALPHA1,			//i=3
		SIGMA2,			//i=4
		ALPHA2,			//i=5
		CORRELATION,	//i=6
		TIMETOEXPIRATION,	//i=7
		A1,				//i=8
		B1,				//i=9	
		K1,				//i=10
		A2,				//i=11
		B2,				//i=12
		K2,				//i=13
		NBSTEPS			//i=14
	};
	enum 
	{ 
		Nb_Parameters =15
	};
	enum
	{ 
		Nb_Derivable_Parameters =14
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
	///   Generic methods (relative to the smile distribution and the copula
	static double Generic_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,BivariateToDoubleFunc* option_payoff);
	static double PowerSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,
				   double a10,double b10,double k10,double a20,double b20,double k20);
	static double Certitude(const ArgumentList& copula,double t,int n);

};




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

