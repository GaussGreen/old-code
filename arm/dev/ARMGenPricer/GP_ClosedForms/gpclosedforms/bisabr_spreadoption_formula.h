/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_spreadoption_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
 
#ifndef _GP_CF_BISABR_SPREADOPTION_FORMULA_H
#define _GP_CF_BISABR_SPREADOPTION_FORMULA_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_BiSABR_SpreadOption_Formula
///			Purpose : Evaluation of spread option under the BiSABR model
///			
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_BiSABR_SpreadOption_Formula
{
	enum ArgumentType
	{
		INDEX1,				// i=0
		ALPHA1,				// i=1
		BETA1,				// i=2
		RHO1,				// i=3
		NU1,				// i=4
		INDEX2,				// i=5
		ALPHA2,				// i=6
		BETA2,				// i=7
		RHO2,				// i=8
		NU2,				// i=9
		RHOS,				// i=10
		RHOV,				// i=11
		RHOC12,				// i=12
		RHOC21,				// i=13
		K,					// i=14
		T,					// i=15
		CALLORPUT,			// i=16
		FLAG				// i=17
	};
	enum 
	{ 
		Nb_Parameters =18
	};

	enum
	{ 
		Nb_Derivable_Parameters =16
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
