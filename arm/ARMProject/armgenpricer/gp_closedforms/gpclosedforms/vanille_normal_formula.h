/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file vanille_normal_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_VANILLE_NORMAL_FORMULA_H
#define _GP_CF_VANILLE_NORMAL_FORMULA_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Vanille_Normal_Formula
///			Purpose : Evaluation of vanilla option under the normal model
///			Assumptions: normal process for the underlying
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Vanille_Normal_Formula
{
	enum ArgumentType
	{
		INDEX,					// i=0
		VOLATILITY ,			// i=1
		STRIKE ,				// i=2
		TIMETOMATURITY,			// i=3
		CALLORPUT				// i=4
	};
	enum 
	{ 
		Nb_Parameters =5
	};

	enum
	{ 
		Nb_Derivable_Parameters =4
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
