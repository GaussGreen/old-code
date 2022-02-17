/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file asian_lognormal_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_ASIAN_LOGNORMAL_FORMULA_H
#define _GP_CF_ASIAN_LOGNORMAL_FORMULA_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result

CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Asian_Vanilla_Lognormal_Formula
///			Purpose : Geman Yor evaluation of Vanilla Asian Option
///			Assumptions: LogNormal hypothesis
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Asian_Vanilla_Lognormal_Formula
{
	enum ArgumentType
	{
		INDEX,				// i=0
		STRIKE,				// i=1
		TIMETOMATURITY,		// i=2
		INTERESTRATE,		// i=3
		VOLATILITY,			// i=4
		CALLORPUT,			// i=5
		ALPHA,				// i=6
		NBTERMS				// i=7
	};
	enum 
	{ 
		Nb_Parameters =8
	};

	enum
	{ 
		Nb_Derivable_Parameters =5
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
