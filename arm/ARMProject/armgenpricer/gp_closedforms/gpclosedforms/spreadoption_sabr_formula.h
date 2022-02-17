/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_sabr_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_SABR_FORMULA_H
#define _GP_CF_SPREADOPTION_SABR_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gpclosedforms/powerspreadoption.h"		/// for instanciation of the templates
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/bismile_templated_pricing.h"


CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula
///			Purpose : Evaluation of spread option cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///			Assumptions: Extended SABR hypothesis fo the factors and a gaussian copula
///
///////////////////////////////////////////////////////////////////////


struct ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula
{
///	typedef ARM_CF_PowerSpreadOption_Formula<SABR_smile,GaussianCopula> structure;
	typedef ARM_CF_Bismiled_PowerSpreadOption_Formula<SABR_smile,SABR_smile,GaussianCopula> structure;
	enum ArgumentType
	{
		INDEX1,
		INDEX2,
		ALPHA1,
		BETA1,
		RHO1,
		NU1,
		ALPHA2,
		BETA2,
		RHO2,
		NU2,
		CORRELATION,
		TIMETOEXPIRATION,
		A1,
		B1,
		K1,
		A2,
		B2,
		K2,
		ALPHA_EXP,
		ALPHA_TANH,
		KB_TANH,
		FLAG,
		NBSTEPS,

		ALGORITHM			/// internal
	};
	

	enum 
	{ 
		Nb_Parameters =24
	};
	enum
	{ 
		Nb_Derivable_Parameters =21
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value2(const ArgumentList& a);
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

