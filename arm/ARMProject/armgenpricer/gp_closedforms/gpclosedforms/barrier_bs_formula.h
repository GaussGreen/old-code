/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file barriere_bs.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_BARRIERE_BS_H
#define _GP_CF_BARRIERE_BS_H

#include "gpbase/port.h"
#include <vector>


CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_BS_SingleBarriere_Formula
///			Purpose : Single European Barriere Option
///			Assumptions: Black and Sholes 
///
///////////////////////////////////////////////////////////////////////


struct ARM_CF_BS_SingleBarriere_Formula
{
	enum ArgumentType
	{
			FORWARD,
			STRIKE,
			BARRIERE,
			REBATE,
			VOLATILITY,
			MATURITY,
			CALLPUT,
			OPTIONTYPE

	};
	enum 
	{ 
		Nb_Parameters =8
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
	static double value2(const ArgumentList& a);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);	

};


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_BS_EuropeanDoubleBarriere_Formula
///			Purpose : Single European Barriere Option
///			Assumptions: Black and Sholes 
///
///////////////////////////////////////////////////////////////////////


struct ARM_CF_BS_DoubleBarriere_Formula
{
	enum ArgumentType
	{
			FORWARD,
			STRIKE,
			UPBARRIERE,
			DOWNBARRIERE,
			VOLATILITY,
			MATURITY,
			CALLPUT

	};
	enum 
	{ 
		Nb_Parameters =7
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value2(const ArgumentList& a);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);	

};



///////////////////////////////////////////////////////////////////////
///  
///			class : BS_PartialTime_Start_SingleBarrier_Formula
///			Purpose : Single European Barriere Option with a premature end of the barrier
///			Assumptions: Black and Sholes 
///
///////////////////////////////////////////////////////////////////////

struct BS_PartialTime_Start_SingleBarrier_Formula
{
	enum ArgumentType
	{
			FORWARD,
			STRIKE,
			BARRIERE,
			REBATE,
			VOLATILITY,
			BARRIERENDTIME,
			MATURITY,
			CALLPUT,
			OPTIONTYPE

	};
	enum 
	{ 
		Nb_Parameters =8
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
	static double value2(const ArgumentList& a);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);	

};

///////////////////////////////////////////////////////////////////////
///  
///			class : BS_SingleBarrier_2Asset_Formula
///			Purpose :  European Option with a barrier on aonther asset
///			Assumptions: Black and Sholes 
///
///////////////////////////////////////////////////////////////////////

struct BS_SingleBarrier_2Asset_Formula
{
	enum ArgumentType
	{
			FORWARD1,
			STRIKE1,
			FORWARD2,
			STRIKE2,
			VOLATILITY1,
			VOLATILITY2,
			CORRELATION,
			MATURITY,
			CALLPUT,
			OPTIONTYPE

	};
	enum 
	{ 
		Nb_Parameters =10
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
	static double value2(const ArgumentList& a);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);	

};

///////////////////////////////////////////////////////////////////////
///  
///			class : BS_PartialTime_Start_SingleBarrier_Formula
///			Purpose : Single European Barriere Option with a premature end of the barrier
///			Assumptions: Black and Sholes 
///
///////////////////////////////////////////////////////////////////////

struct BS_PartialTime_End_SingleBarrier_Formula
{
	enum ArgumentType
	{
			FORWARD,
			STRIKE,
			BARRIERE,
			REBATE,
			VOLATILITY,
			BARRIERSTARTTIME,
			MATURITY,
			CALLPUT,
			OPTIONTYPE

	};
	enum 
	{ 
		Nb_Parameters =9
	};
	enum OptionType
	{
		CROSS_AND_IN,
		CROSS_AND_OUT,
		INSIDE_UP_AND_IN,
		INSIDE_UP_AND_OUT,
		INSIDE_DOWN_AND_IN,
		INSIDE_DOWN_AND_OUT
	};


	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value2(const ArgumentList& a);
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
