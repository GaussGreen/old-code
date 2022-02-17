/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_normal_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_normal.h"
#include "gpclosedforms/spreadoption_normal_formula.h"


#include <glob/expt.h>

CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_SpreadDigitalOption_N_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////

double ARM_CF_SpreadDigitalOption_N_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_SpreadDigitalOption_N_Formula::value,i,a,s);
}


double ARM_CF_SpreadDigitalOption_N_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_SpreadDigitalOption_N_Formula::value,i,j,a,s1,s2);
}


 double ARM_CF_SpreadDigitalOption_N_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_SpreadDigitalOption_N_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SpreadDigitalOption_N_Formula::value  : bad argsize");
	 }
	 double val= SpreadDigitalOption_N(
			 a[ARM_CF_SpreadDigitalOption_N_Formula::INDEX1],	//INDEX1
			 a[ARM_CF_SpreadDigitalOption_N_Formula::INDEX2],	//INDEX2
			 a[ARM_CF_SpreadDigitalOption_N_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_SpreadDigitalOption_N_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_SpreadDigitalOption_N_Formula::CORRELATION],	//CORRELATION
			 a[ARM_CF_SpreadDigitalOption_N_Formula::STRIKE],	//STRIKE
			 a[ARM_CF_SpreadDigitalOption_N_Formula::TIMETORESET],	//t
			 a[ARM_CF_SpreadDigitalOption_N_Formula::CALLORPUT],	//CALLORPUT
			 a[ARM_CF_SpreadDigitalOption_N_Formula::TYPE]	// TYPE
			);
		return val;
	}



 double ARM_CF_SpreadDigitalOption_N_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_SpreadDigitalOption_N_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_SpreadDigitalOption_N_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_SpreadDigitalOption_N_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_SpreadDigitalOption_N_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_SpreadDigitalOption_N_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_SpreadDigitalOption_N_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_SpreadDigitalOption_N_Formula::STRIKE :
			return 0.00001;		// Strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SpreadDigitalOption_N_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_SpreadDigitalOption_N_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_SpreadDigitalOption_N_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments : <7 ");
	 }
	 
	 if (a[ARM_CF_SpreadDigitalOption_N_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_SpreadDigitalOption_N_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_SpreadDigitalOption_N_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_SpreadDigitalOption_N_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1
	 	
	 if (a[ARM_CF_SpreadDigitalOption_N_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
		
	 if (abs(a[ARM_CF_SpreadDigitalOption_N_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION) not <=1"); //		Abs(CORRELATION)<=1
	
	 if ((a[ARM_CF_SpreadDigitalOption_N_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_SpreadDigitalOption_N_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_SpreadDigitalOption_N_Formula::check_dimension(int rank)
{
	 if ((rank<0)||(rank>>=ARM_CF_SpreadDigitalOption_N_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : Vega1SpreadDigitalOption_N_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////


 double ARM_CF_Vega1SpreadDigitalOption_N_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_SpreadDigitalOption_N_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Vega1SpreadDigitalOption_N_Formula::value  : bad argsize");
	 }
	 double val= Vega1SpreadDigitalOption_N(
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::INDEX1],	//INDEX1
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::INDEX2],	//INDEX2
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::CORRELATION],	//CORRELATION
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::STRIKE],	//STRIKE
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::TIMETORESET],	//t
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::CALLORPUT],	//CALLORPUT
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::TYPE]	//	TYPE
			);
		return val;
	}
	double ARM_CF_Vega1SpreadDigitalOption_N_Formula::value(int i,const ArgumentList& a, double s)
	{
		return standard_first_derivative(ARM_CF_Vega1SpreadDigitalOption_N_Formula::value,i,a,s);
	}
	double ARM_CF_Vega1SpreadDigitalOption_N_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
	{
		
		return standard_second_derivative(ARM_CF_Vega1SpreadDigitalOption_N_Formula::value,i,j,a,s1,s2);
		
	}
	

	double ARM_CF_Vega1SpreadDigitalOption_N_Formula::specific_shift(int i) 
		{
		switch(i)
		{
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_Vega1SpreadDigitalOption_N_Formula::STRIKE :
			return 0.00001;		// Strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Vega1SpreadDigitalOption_N_Formula::specific_shift : incorrect input");
		}
	}
	
	ArgumentList_Checking_Result ARM_CF_Vega1SpreadDigitalOption_N_Formula::check_argument(const ArgumentList& a)
		{
		if(a.size()!=ARM_CF_Vega1SpreadDigitalOption_N_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments : ");
	 }
	 
	 if (a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1
	 	
	 if (a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
		
	 if (abs(a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION) not <=1"); //		Abs(CORRELATION)<=1
	
	 if ((a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); //	
 	return ArgumentList_Checking_Result(true,string(""));
 }


	ArgumentList_Checking_Result ARM_CF_Vega1SpreadDigitalOption_N_Formula::check_dimension(int rank)
		{
		if ((rank<0)||(rank>=ARM_CF_Vega1SpreadDigitalOption_N_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}


//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_Vega2SpreadDigitalOption_N_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////
	double ARM_CF_Vega2SpreadDigitalOption_N_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_SpreadDigitalOption_N_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Vega2SpreadDigitalOption_N_Formula::value  : bad argsize");
	 }
	 double val= Vega2SpreadDigitalOption_N(
			  a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::INDEX1],	//INDEX1
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::INDEX2],	//INDEX2
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::CORRELATION],	//CORRELATION
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::STRIKE],	//STRIKE
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::TIMETORESET],	//t
			 a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::CALLORPUT],	//CALLORPUT
			 a[ARM_CF_Vega1SpreadDigitalOption_N_Formula::TYPE]	//	TYPE
			);
		return val;
	}


	double ARM_CF_Vega2SpreadDigitalOption_N_Formula::value(int i,const ArgumentList& a, double s)
	{
		return standard_first_derivative(ARM_CF_Vega2SpreadDigitalOption_N_Formula::value,i,a,s);
		
	}
	double ARM_CF_Vega2SpreadDigitalOption_N_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
	{
		
		return standard_second_derivative(ARM_CF_Vega2SpreadDigitalOption_N_Formula::value,i,j,a,s1,s2);
		
	}


	double ARM_CF_Vega2SpreadDigitalOption_N_Formula::specific_shift(int i) 
		{
		switch(i)
		{
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_Vega2SpreadDigitalOption_N_Formula::STRIKE :
			return 0.00001;		// Strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Vega2SpreadDigitalOption_N_Formula::specific_shift : incorrect input");
		}
	}
	
	ArgumentList_Checking_Result ARM_CF_Vega2SpreadDigitalOption_N_Formula::check_argument(const ArgumentList& a)
		{
		if(a.size()!=ARM_CF_Vega2SpreadDigitalOption_N_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments : !=8 ");
	 }
	 
	 if (a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1
	 	
	 if (a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
		
	 if (abs(a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION) not <=1"); //		Abs(CORRELATION)<=1
	
	 if ((a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_Vega2SpreadDigitalOption_N_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


	ArgumentList_Checking_Result ARM_CF_Vega2SpreadDigitalOption_N_Formula::check_dimension(int rank)
		{
		if ((rank<0)||(rank>=ARM_CF_Vega2SpreadDigitalOption_N_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}



CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

