/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_LOGNORMAL_FORMULA_H
#define _GP_CF_SPREADOPTION_LOGNORMAL_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result

#include "product.h"			/// for typedef manipulation !
#include "sum.h"
#include "difference.h"
#include "add_submodel.h"
#include "change_numeraire.h"

CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_SpreadDigitalOption_Formula
///			Purpose : Evaluation of spread option cashflow=  1{S1-S2-K>0}
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_SpreadDigitalOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,
		INDEX1,
		INDEX2,
		VOLATILITY1,
		VOLATILITY2,
		CORRELATION,
        STRIKE,
		CALLORPUT,
		NBSTEPS
	};
	enum 
	{ 
		Nb_Parameters =9
	};
		enum 
	{ 
		Nb_Derivable_Parameters =7
	};
	
	enum DerivedStruct			// the derived type that are built using typedefs
	{
		DIGITALOPTION,
		SPREADOPTION,
		PAYFIRST,
		PAYSECOND
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Vega1SpreadOption_Formula
///			Purpose : Evaluation of spread option cashflow=  1{S1-S2-K>0}
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Vega1SpreadOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,
		INDEX1,
		INDEX2,
		VOLATILITY1,
		VOLATILITY2,
		CORRELATION,
        STRIKE,
		OPTIONTYPE,
		NBSTEPS
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

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Vega2SpreadOption_Formula
///			Purpose : Evaluation of spread option cashflow=  1{S1-S2-K>0}
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Vega2SpreadOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,
		INDEX1,
		INDEX2,
		VOLATILITY1,
		VOLATILITY2,
		CORRELATION,
        STRIKE,
		CALLORPUT,
		NBSTEPS
	};
	enum 
	{ 
		Nb_Parameters =9
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const  ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Smiled_SpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  1{S1-S2-K>0}
///									 cashflow put  =  1{S1-S2-K<0}
///			Assumptions: Lognormal hypothesis, 
///         take into acount the smile of the spread option
///
///////////////////////////////////////////////////////////////////////


typedef 	Sum3<Add_Input<ARM_CF_SpreadDigitalOption_Formula,2>,			// add the 2 smile slope parameter
			Product<
				Add_Input<ARM_CF_Vega1SpreadOption_Formula,2>,
				ARM_CF_SpreadDigitalOption_Formula::Nb_Parameters     // parameter 1 added = smile slope 1
				>,
			Product<
				Add_Input<ARM_CF_Vega2SpreadOption_Formula,2>,
				ARM_CF_SpreadDigitalOption_Formula::Nb_Parameters+1     // parameter 2 added = smile slope 2
				>
			>

								
			ARM_CF_Smiled_SpreadDigitalOption_Formula;
	   

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index1Paying_SpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S1*1{S1-S2-K>0}
///									 cashflow put  =  S1*1{S1-S2-K<0}
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////


typedef 	Product<Add_SubModel<Add_SubModel<ARM_CF_SpreadDigitalOption_Formula,
				Change_Numeraire1<		
					ARM_CF_SpreadDigitalOption_Formula::INDEX1,			// on modifie l'input S1		
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,	// multiplicativement par l'exp(sig1^2*T)
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_SpreadDigitalOption_Formula::INDEX2,			// on modifie l'input S2
					ARM_CF_SpreadDigitalOption_Formula::CORRELATION,	// multiplicativement par exp(sig1*sig2*corr*T)
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
								>,
				ARM_CF_SpreadDigitalOption_Formula::INDEX1>
								
	   ARM_CF_Index1Paying_SpreadDigitalOption_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index1Paying_Smiled_SpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S1*1{S1-S2-K>0}
///									 cashflow put  =  S1*1{S1-S2-K<0}
///			Assumptions: Lognormal hypothesis
///			smile taken into account : smile slopes should be provided
///
///////////////////////////////////////////////////////////////////////


typedef 	Product<Add_SubModel<Add_SubModel<ARM_CF_Smiled_SpreadDigitalOption_Formula,
				Change_Numeraire1<
					ARM_CF_SpreadDigitalOption_Formula::INDEX1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_SpreadDigitalOption_Formula::INDEX2,
					ARM_CF_SpreadDigitalOption_Formula::CORRELATION,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
								>,
				ARM_CF_SpreadDigitalOption_Formula::INDEX1>
								
	   ARM_CF_Index1Paying_Smiled_SpreadDigitalOption_Formula;


///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index2Paying_SpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S2*1{S1-S2-K>0}
///									 cashflow put  =  S2*1{S1-S2-K<0}
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

 
typedef 	Product<Add_SubModel<Add_SubModel<ARM_CF_SpreadDigitalOption_Formula,
				Change_Numeraire1<
					ARM_CF_SpreadDigitalOption_Formula::INDEX2,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_SpreadDigitalOption_Formula::INDEX1,
					ARM_CF_SpreadDigitalOption_Formula::CORRELATION,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
								>,
				ARM_CF_SpreadDigitalOption_Formula::INDEX2>
								
	   ARM_CF_Index2Paying_SpreadDigitalOption_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index2Paying_Smiled_SpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S2*1{S1-S2-K>0}
///									 cashflow put  =  S2*1{S1-S2-K<0}
///			Assumptions: Lognormal hypothesis
///			smile taken into account : smile slopes should be provided
///
///////////////////////////////////////////////////////////////////////

typedef 	Product<Add_SubModel<Add_SubModel<ARM_CF_Smiled_SpreadDigitalOption_Formula,
				Change_Numeraire1<
					ARM_CF_SpreadDigitalOption_Formula::INDEX2,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_SpreadDigitalOption_Formula::INDEX1,
					ARM_CF_SpreadDigitalOption_Formula::CORRELATION,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1> 
								>,
				ARM_CF_SpreadDigitalOption_Formula::INDEX2>
								
	   ARM_CF_Index2Paying_Smiled_SpreadDigitalOption_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_SpreadOption_Formula
///			Purpose : Evaluation of  cashflow call =  (S1-S2-K)*1{S1-S2-K>0}
///									 cashflow put  =  (K-S1+S2)*1{S1-S2-K<0}
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

 
typedef		Product<
				Difference<Difference<
							ARM_CF_Index1Paying_SpreadDigitalOption_Formula,
							ARM_CF_Index2Paying_SpreadDigitalOption_Formula>,
					   Product<
							ARM_CF_SpreadDigitalOption_Formula,
							ARM_CF_SpreadDigitalOption_Formula::STRIKE>
					  >,
				ARM_CF_SpreadDigitalOption_Formula::CALLORPUT>
			ARM_CF_SpreadOption_Formula;
		

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

