/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_integrals.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_CEV_FORMULA_H
#define _GP_CF_CEV_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result

CC_BEGIN_NAMESPACE(ARM)

struct ARM_CF_CEV_VanillaOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,		// i=0
		INDEX,				// i=1
		STRIKE,				// i=2
		VOLATILITY,			// i=3
		BETA,				// i=4
		MU,					// i=5
		CALLORPUT,			// i=7
		NBTERMS				// i=8
	};
	enum 
	{ 
		Nb_Parameters =8
	};
	enum 
	{ 
		Nb_Derivable_Parameters =6		
	};

	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};


struct ARM_CF_CEV_DoubleBarrierOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,		// i=0
		INDEX,				// i=1
		STRIKE,				// i=2		
		VOLATILITY,			// i=3
		BETA,				// i=4
		MU,					// i=5
		BARRIERDOWN,		// i=6
		BARRIERUP,			// i=7
		CALLORPUT,			// i=8
		NBTERMS				// i=9
	};
	enum 
	{ 
		Nb_Parameters =10
	};
	enum 
	{ 
		Nb_Derivable_Parameters =8
	};

	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};

struct ARM_CF_CEV_BarrierOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,		// i=0
		INDEX,				// i=1
		STRIKE,				// i=2		
		VOLATILITY,			// i=3
		BETA,				// i=4
		MU,					// i=5
		BARRIER,			// i=6
		OPTIONTYPE,			// i=7
		CALLORPUT,			// i=8
		NBTERMS				// i=9
	};
	enum 
	{ 
		Nb_Parameters =10
	};
	enum
	{ 
		Nb_Derivable_Parameters =7
	};
	enum OptionType
	{
		DOWN_AND_IN,
		UP_AND_IN,
		DOWN_AND_OUT,
		UP_AND_OUT
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


