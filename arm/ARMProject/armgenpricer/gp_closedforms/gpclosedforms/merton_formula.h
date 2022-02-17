/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file merton.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_MERTON_FORMULA_H
#define _GP_CF_MERTON_FORMULA_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Merton_JumpDiffusion_Formula
///			Purpose : Evaluation of vanilla option under the Merton model
///			Assumptions: Normal Jumps of the log of the underlying
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Merton_JumpDiffusion_Formula
{
	enum ArgumentType
	{
		INDEX,				// i=0
		STRIKE,				// i=1
		TIMETOMATURITY,		// i=2
		VOLATILITY,			// i=3
		JUMPPROBABILITY,	// i=4
		JUMPMEAN,			// i=5
		JUMPVOLATILITY,		// i=6
		CALLORPUT,			// i=7
		NBTERMS				// i=8
	};
	enum 
	{ 
		Nb_Parameters =9
	};

	enum
	{ 
		Nb_Derivable_Parameters =7
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
