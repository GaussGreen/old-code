/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_glambda_student_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_GLAMBDA_STUDENT_FORMULA_H
#define _GP_CF_SPREADOPTION_GLAMBDA_STUDENT_FORMULA_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gpclosedforms/powerspreadoption.h"		/// for instanciation of the templates
#include "gpclosedforms/smile_glambda.h"


CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Glambda_Student_PowerSpreadOption_Formula
///			Purpose : Evaluation of spread option cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///			Assumptions: Extended SABR hypothesis fo the factors and a student copula
///
///////////////////////////////////////////////////////////////////////


struct ARM_CF_GLambda_Student_PowerSpreadOption_Formula
{
	typedef ARM_CF_PowerSpreadOption_Formula<GLambda_Smile,StudentCopula> structure;
	enum ArgumentType
	{
		FIRST_UNDERLYING_L1,
		FIRST_UNDERLYING_L2,
		FIRST_UNDERLYING_L3,
		FIRST_UNDERLYING_L4,
		FIRST_UNDERLYING_L5,
		SECOND_UNDERLYING_L1,
		SECOND_UNDERLYING_L2,
		SECOND_UNDERLYING_L3,
		SECOND_UNDERLYING_L4,
		SECOND_UNDERLYING_L5,
		CORRELATION,
		DEGRE,
		A1,
		B1,
		K1,
		A2,
		B2,
		K2,
		NBSTEPS
	};
	
	enum 
	{ 
		Nb_Parameters =19
	};
	enum
	{ 
		Nb_Derivable_Parameters =18
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

