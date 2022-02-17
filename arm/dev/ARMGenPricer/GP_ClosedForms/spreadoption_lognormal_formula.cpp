/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_lognormal.h"

#include <glob/expt.h>

 
CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_SpreadDigitalOption_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////

double ARM_CF_SpreadDigitalOption_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_SpreadDigitalOption_Formula::value,i,a,s);

}


double ARM_CF_SpreadDigitalOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_SpreadDigitalOption_Formula::value,i,j,a,s1,s2);

}


 double ARM_CF_SpreadDigitalOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_SpreadDigitalOption_Formula::Nb_Parameters) 
	 {
		 throw ("ARM_CF_SpreadDigitalOption_Formula::value  : bad argsize");
	 }
	 
	 double val= SpreadDigitalOption(
			 a[ARM_CF_SpreadDigitalOption_Formula::INDEX1],	//INDEX1
			 a[ARM_CF_SpreadDigitalOption_Formula::INDEX2],	//INDEX2
			 a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_SpreadDigitalOption_Formula::CORRELATION],	//CORRELATION
			 a[ARM_CF_SpreadDigitalOption_Formula::STRIKE],	//STRIKE
			 a[ARM_CF_SpreadDigitalOption_Formula::TIMETORESET],	//t
			 a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT],	//CALLPUT
			 a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS]	//NBSTEPS
			);
		return val;
	}



 double ARM_CF_SpreadDigitalOption_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_SpreadDigitalOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_SpreadDigitalOption_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_SpreadDigitalOption_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_SpreadDigitalOption_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_SpreadDigitalOption_Formula::STRIKE :
			return 0.00001;		// Strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SpreadDigitalOption_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_SpreadDigitalOption_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_SpreadDigitalOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_SpreadDigitalOption_Formula::TIMETORESET]<0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_SpreadDigitalOption_Formula::INDEX1]<0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_SpreadDigitalOption_Formula::INDEX2]<0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1
	 	
	 if (a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
		
	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," fabs(CORRELATION) not <=1"); //		fabs(CORRELATION)<=1

	 if ((a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_SpreadDigitalOption_Formula:CALLORPUT should be CALL or PUT"); //			Should be Call or Put

	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS])<=1) return ArgumentList_Checking_Result(false," NBSTEPS should be >= 1"); //		fabs(NBSTEPS)>=1
	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS])>=200) return ArgumentList_Checking_Result(false," NBSTEPS should be <= 200"); //		fabs(NBSTEPS)<=200


 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_SpreadDigitalOption_Formula::check_dimension(int rank)
{
	 if ((rank<0)||(rank>=ARM_CF_SpreadDigitalOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

#undef ARM_CF_SQRT_2_PI
#undef ARM_CF_SQRT2
#undef ARM_CF_INVSQRTPI


//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_Vega1SpreadOption_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////
	double ARM_CF_Vega1SpreadOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	if (argsize!=ARM_CF_Vega1SpreadOption_Formula::Nb_Parameters) 
	 {
		 throw ("ARM_CF_SpreadDigitalOption_Formula::value  : bad argsize");
	 }
	 double val= Vega1SpreadOption(
			 a[ARM_CF_SpreadDigitalOption_Formula::INDEX1],	//INDEX1
			 a[ARM_CF_SpreadDigitalOption_Formula::INDEX2],	//INDEX2
			 a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_SpreadDigitalOption_Formula::CORRELATION],	//CORRELATION
			 a[ARM_CF_SpreadDigitalOption_Formula::STRIKE],	//STRIKE
			 a[ARM_CF_SpreadDigitalOption_Formula::TIMETORESET],	//t
			 a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT],	//CALLPUT
			 a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS]	//NBSTEPS
			);
		return val;
	}
	double ARM_CF_Vega1SpreadOption_Formula::value(int i,const ArgumentList& a, double s)
	{
		return standard_first_derivative(ARM_CF_Vega1SpreadOption_Formula::value,i,a,s);

}
	double ARM_CF_Vega1SpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
		{
	
	return standard_second_derivative(ARM_CF_Vega1SpreadOption_Formula::value,i,j,a,s1,s2);

}


	double ARM_CF_Vega1SpreadOption_Formula::specific_shift(int i) 
		{
		switch(i)
		{
		case ARM_CF_Vega1SpreadOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_Vega1SpreadOption_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_Vega1SpreadOption_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_Vega1SpreadOption_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_Vega1SpreadOption_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_Vega1SpreadOption_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_Vega1SpreadOption_Formula::STRIKE :
			return 0.00001;		// Strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Vega1SpreadOption_Formula::specific_shift : incorrect input");
		}
	}
	
	ArgumentList_Checking_Result ARM_CF_Vega1SpreadOption_Formula::check_argument(const ArgumentList& a)
		{	 
	 if(a.size()!=ARM_CF_SpreadDigitalOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_SpreadDigitalOption_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_SpreadDigitalOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_SpreadDigitalOption_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1
	 	
	 if (a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
		
	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," fabs(CORRELATION) not <=1"); //		fabs(CORRELATION)<=1
	
	 if ((a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_Vega1SpreadOption_Formula:CALLORPUT should be CALL or PUT"); //			Call or Put Flag 

	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS])<=1) return ArgumentList_Checking_Result(false," NBSTEPS should be >= 1"); //		fabs(NBSTEPES)>=1
	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS])>=200) return ArgumentList_Checking_Result(false," NBSTEPS should be <= 200"); //		fabs(NBSTEPES)>=1

 	return ArgumentList_Checking_Result(true,string(""));
 }


	ArgumentList_Checking_Result ARM_CF_Vega1SpreadOption_Formula::check_dimension(int rank)
		{
	if ((rank<0)||(rank>=ARM_CF_SpreadDigitalOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}


	//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_Vega2SpreadOption_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////
	double ARM_CF_Vega2SpreadOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_Vega2SpreadOption_Formula::Nb_Parameters) 
	 {
		 throw ("ARM_CF_SpreadDigitalOption_Formula::value  : bad argsize");
	 }
	 double val= Vega2SpreadOption(
			a[ARM_CF_SpreadDigitalOption_Formula::INDEX1],	//INDEX1
			 a[ARM_CF_SpreadDigitalOption_Formula::INDEX2],	//INDEX2
			 a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_SpreadDigitalOption_Formula::CORRELATION],	//CORRELATION
			 a[ARM_CF_SpreadDigitalOption_Formula::STRIKE],	//STRIKE
			 a[ARM_CF_SpreadDigitalOption_Formula::TIMETORESET],	//t
			 a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT],	//CALLPUT
			 a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS]	//NBSTEPS
			);
		return val;
	}


	double ARM_CF_Vega2SpreadOption_Formula::value(int i,const ArgumentList& a, double s)
	{
		return standard_first_derivative(ARM_CF_Vega2SpreadOption_Formula::value,i,a,s);

}
	double ARM_CF_Vega2SpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
		{
	
	return standard_second_derivative(ARM_CF_Vega2SpreadOption_Formula::value,i,j,a,s1,s2);

}


	double ARM_CF_Vega2SpreadOption_Formula::specific_shift(int i) 
		{
		switch(i)
		{
		case ARM_CF_Vega2SpreadOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_Vega2SpreadOption_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_Vega2SpreadOption_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_Vega2SpreadOption_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_Vega2SpreadOption_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_Vega2SpreadOption_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_Vega2SpreadOption_Formula::STRIKE :
			return 0.00001;		// Strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Vega1SpreadOption_Formula::specific_shift : incorrect input");
		}
	}
	
	ArgumentList_Checking_Result ARM_CF_Vega2SpreadOption_Formula::check_argument(const ArgumentList& a)
		{
	 if(a.size()!=ARM_CF_SpreadDigitalOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_SpreadDigitalOption_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_SpreadDigitalOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_SpreadDigitalOption_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1
	 	
	 if (a[ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
		
	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," fabs(CORRELATION) not <=1"); //		fabs(CORRELATION)<=1
	
	 if ((a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_SpreadDigitalOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_Vega2SpreadOption_Formula:CALLORPUT should be CALL or PUT"); //			positivity of STRIKE

	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS])<=1) return ArgumentList_Checking_Result(false," NBSTEPS should be >= 1"); //		fabs(NBSTEPS)>=1
	 if (fabs(a[ARM_CF_SpreadDigitalOption_Formula::NBSTEPS])>=200) return ArgumentList_Checking_Result(false," NBSTEPS should be <= 200"); //		fabs(NBSTEPS)<=200

 	return ArgumentList_Checking_Result(true,string(""));
 }


	ArgumentList_Checking_Result ARM_CF_Vega2SpreadOption_Formula::check_dimension(int rank)
		{
	if ((rank<0)||(rank>=ARM_CF_SpreadDigitalOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




CC_END_NAMESPACE()
 
#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_SQRT2 
#undef ARM_CF_INVSQRTPI

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/