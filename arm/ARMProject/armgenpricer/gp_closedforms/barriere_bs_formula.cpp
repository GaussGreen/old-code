/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\Barriere_bs_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/Barriere_bs.h"
#include "gpclosedforms/Barriere_bs_formula.h"

CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_BS_SingleBarriere_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////



double ARM_CF_BS_SingleBarriere_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BS_SingleBarriere_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_SingleBarriere_Formula::value  : bad argsize");
	 }
	 double val= BS_SingleBarrierOption(
			 a[ARM_CF_BS_SingleBarriere_Formula::FORWARD],	//INDEX
			 a[ARM_CF_BS_SingleBarriere_Formula::STRIKE],		//STRIKE
			 a[ARM_CF_BS_SingleBarriere_Formula::BARRIERE],	//BARRIERE
			 a[ARM_CF_BS_SingleBarriere_Formula::REBATE],		//REBATE
			  a[ARM_CF_BS_SingleBarriere_Formula::MATURITY],	//MATURITY
			 a[ARM_CF_BS_SingleBarriere_Formula::VOLATILITY],	//VOLATILITY
			 a[ARM_CF_BS_SingleBarriere_Formula::DISCOUNT],		// DISCOUNT									//discountRate 
			 0.0,												//dividend
			 a[ARM_CF_BS_SingleBarriere_Formula::CALLPUT],	//CALLPUT
			 a[ARM_CF_BS_SingleBarriere_Formula::OPTIONTYPE]	//OPTIONTYPE
			);
		return val;
	}

double ARM_CF_BS_SingleBarriere_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_BS_SingleBarriere_Formula::value,i,a,s);

}


double ARM_CF_BS_SingleBarriere_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_BS_SingleBarriere_Formula::value,i,j,a,s1,s2);

}





 double ARM_CF_BS_SingleBarriere_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_BS_SingleBarriere_Formula::FORWARD :
			return 0.00001;		// time to maturity
		case ARM_CF_BS_SingleBarriere_Formula::STRIKE :
			return 0.00001;		// Forward index 1
		case ARM_CF_BS_SingleBarriere_Formula::BARRIERE :
			return 0.00001;		// Forward index 2
		case ARM_CF_BS_SingleBarriere_Formula::REBATE :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_BS_SingleBarriere_Formula::VOLATILITY :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_SingleBarriere_Formula::DISCOUNT :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_SingleBarriere_Formula::MATURITY :
			return 0.00001;		// Correlation between the two indexes
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_SingleBarriere_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_BS_SingleBarriere_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_BS_SingleBarriere_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 
	 if (a[ARM_CF_BS_SingleBarriere_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," FORWARD not positive"); //		positivity of FORWARD

	 if (a[ARM_CF_BS_SingleBarriere_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	 
	 if (a[ARM_CF_BS_SingleBarriere_Formula::BARRIERE]<=0) return ArgumentList_Checking_Result(false," BARRIERE not positive"); //			positivity of BARRIERE

	 if (a[ARM_CF_BS_SingleBarriere_Formula::REBATE ]<0) return ArgumentList_Checking_Result(false," REBATE not positive"); //		positivity of REBATE
	 	
	 if (a[ARM_CF_BS_SingleBarriere_Formula::VOLATILITY]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY
		
	 if (a[ARM_CF_BS_SingleBarriere_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," MATURITY not positive"); //		positivity of MATURITY

	 if ((a[ARM_CF_BS_SingleBarriere_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN) &&
		 (a[ARM_CF_BS_SingleBarriere_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN) &&
		 (a[ARM_CF_BS_SingleBarriere_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT) &&
		 (a[ARM_CF_BS_SingleBarriere_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT))
		  return ArgumentList_Checking_Result(false," OPTIONTYPE should be DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT or UP_AND_OUT"); //	checking of option type

	 if ((a[ARM_CF_BS_SingleBarriere_Formula::CALLPUT]!=K_CALL) &&
		  (a[ARM_CF_BS_SingleBarriere_Formula::CALLPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLPUT should be K_CALL or K_PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_BS_SingleBarriere_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>ARM_CF_BS_SingleBarriere_Formula::Nb_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_BS_DoubleBarriere_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////



double ARM_CF_BS_DoubleBarriere_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BS_DoubleBarriere_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_DoubleBarriere_Formula::value  : bad argsize");
	 }
	 double val= BS_DoubleBarrierOption(
			 a[ARM_CF_BS_DoubleBarriere_Formula::FORWARD],	//INDEX
			 a[ARM_CF_BS_DoubleBarriere_Formula::STRIKE],		//STRIKE
			 
			 a[ARM_CF_BS_DoubleBarriere_Formula::DOWNBARRIERE],	//DOWNBARRIERE
			 a[ARM_CF_BS_DoubleBarriere_Formula::UPBARRIERE],	//UPBARRIERE
			 a[ARM_CF_BS_DoubleBarriere_Formula::MATURITY],		//MATURITY
			 a[ARM_CF_BS_DoubleBarriere_Formula::VOLATILITY],	//VOLATILITY
			 a[ARM_CF_BS_DoubleBarriere_Formula::INTERESTRATE],												//discountRate
			 a[ARM_CF_BS_DoubleBarriere_Formula::DIVIDENDYIELD],												//dividend
			 0.0,												// convexity of the lower barrier
			 0.0,												// convexity of the upper barrier
			 a[ARM_CF_BS_DoubleBarriere_Formula::CALLPUT],	//CALLPUT
			 10													// nb of terms code en dur
			);
		return val;
	}

double ARM_CF_BS_DoubleBarriere_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_BS_DoubleBarriere_Formula::value,i,a,s);

}


double ARM_CF_BS_DoubleBarriere_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_BS_DoubleBarriere_Formula::value,i,j,a,s1,s2);

}





 double ARM_CF_BS_DoubleBarriere_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_BS_DoubleBarriere_Formula::FORWARD :
			return 0.00001;		// time to maturity
		case ARM_CF_BS_DoubleBarriere_Formula::STRIKE :
			return 0.00001;		// Forward index 1
		case ARM_CF_BS_DoubleBarriere_Formula::UPBARRIERE :
			return 0.00001;		// Forward index 2
		case ARM_CF_BS_DoubleBarriere_Formula::DOWNBARRIERE :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_BS_DoubleBarriere_Formula::VOLATILITY :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_DoubleBarriere_Formula::MATURITY :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_BS_DoubleBarriere_Formula::INTERESTRATE :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_BS_DoubleBarriere_Formula::DIVIDENDYIELD :
			return 0.00001;		// Correlation between the two indexes
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_DoubleBarriere_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_BS_DoubleBarriere_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_BS_DoubleBarriere_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 
	 if (a[ARM_CF_BS_DoubleBarriere_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," FORWARD not positive"); //		positivity of FORWARD

	 if (a[ARM_CF_BS_DoubleBarriere_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	 
	 if (a[ARM_CF_BS_DoubleBarriere_Formula::UPBARRIERE]<=0) return ArgumentList_Checking_Result(false," UPBARRIERE not positive"); //			positivity of UPBARRIERE
	
	 if (a[ARM_CF_BS_DoubleBarriere_Formula::DOWNBARRIERE ]<=0) return ArgumentList_Checking_Result(false," DOWNBARRIERE not positive"); //		positivity of DOWNBARRIERE
	 	
	 if (a[ARM_CF_BS_DoubleBarriere_Formula::VOLATILITY]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY
		
	 if (a[ARM_CF_BS_DoubleBarriere_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," MATURITY not positive"); //		positivity of MATURITY

	 if ((a[ARM_CF_BS_DoubleBarriere_Formula::CALLPUT]!=K_CALL) &&
		  (a[ARM_CF_BS_DoubleBarriere_Formula::CALLPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLPUT should be K_CALL or K_PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_BS_DoubleBarriere_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>ARM_CF_BS_DoubleBarriere_Formula::Nb_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}



//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_Start_SingleBarrier_Formula::value  : bad argsize");
	 }
	 double val=BS_PartialTime_Start_SingleBarrierOption(
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::FORWARD],	//INDEX
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::STRIKE],		//STRIKE
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::BARRIERE],	//BARRIERE
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::REBATE],		//REBATE
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::BARRIERENDTIME],	//MATURITY
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::MATURITY],	//MATURITY
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::VOLATILITY],	//VOLATILITY
			 0.0,												//discountRate 
			 0.0,												//dividend
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::CALLPUT],	//CALLPUT
			 a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::OPTIONTYPE]	//OPTIONTYPE
			);
		return val;
	}

double ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::value,i,a,s);

}


double ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::value,i,j,a,s1,s2);

}





 double ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::FORWARD :
			return 0.00001;		// time to maturity
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::STRIKE :
			return 0.00001;		// Forward index 1
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::BARRIERE :
			return 0.00001;		// Forward index 2
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::REBATE :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::VOLATILITY :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::BARRIERENDTIME :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::MATURITY :
			return 0.00001;		// Correlation between the two indexes
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 
	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," FORWARD not positive"); //		positivity of FORWARD

	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	 
	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::BARRIERE]<=0) return ArgumentList_Checking_Result(false," BARRIERE not positive"); //			positivity of BARRIERE
	
	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::REBATE ]<0) return ArgumentList_Checking_Result(false," REBATE not positive"); //		positivity of REBATE
	 	
	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::VOLATILITY]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY
		
	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," MATURITY not positive"); //		positivity of MATURITY

	 	 if (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::BARRIERENDTIME]<=0) return ArgumentList_Checking_Result(false," BARRIERENDTIME not positive"); //		positivity of BARRIERENDTIME

	 if ((a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::DOWN_AND_IN) &&
		 (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::UP_AND_IN) &&
		 (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::DOWN_AND_OUT) &&
		 (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::UP_AND_OUT))
		  return ArgumentList_Checking_Result(false," OPTIONTYPE should be DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT or UP_AND_OUT"); //	checking of option type

	 if ((a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::CALLPUT]!=K_CALL) &&
		  (a[ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::CALLPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLPUT should be K_CALL or K_PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::Nb_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

 //////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_BS_PartialTime_End_SingleBarrier_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_Start_SingleBarrier_Formula::value  : bad argsize");
	 }
	 double val=BS_PartialTime_End_SingleBarrierOption(
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::FORWARD],	//INDEX
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::STRIKE],		//STRIKE
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::BARRIERE],	//BARRIERE
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::REBATE],		//REBATE
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::BARRIERSTARTTIME],	//BARRIERSTARTTIME
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::MATURITY],	//MATURITY
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::VOLATILITY],	//VOLATILITY
			 0.0,												//discountRate 
			 0.0,												//dividend
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CALLPUT],	//CALLPUT
			 a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]	//OPTIONTYPE
			);
		return val;
	}

double ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::value,i,a,s);

}


double ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::value,i,j,a,s1,s2);

}





 double ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::FORWARD :
			return 0.00001;		// time to maturity
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::STRIKE :
			return 0.00001;		// Forward index 1
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::BARRIERE :
			return 0.00001;		// Forward index 2
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::REBATE :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::VOLATILITY :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::BARRIERSTARTTIME :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::MATURITY :
			return 0.00001;		// Correlation between the two indexes
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 
	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," FORWARD not positive"); //		positivity of FORWARD

	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	 
	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::BARRIERE]<=0) return ArgumentList_Checking_Result(false," BARRIERE not positive"); //			positivity of BARRIERE
	
	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::REBATE ]<0) return ArgumentList_Checking_Result(false," REBATE not positive"); //		positivity of REBATE
	 	
	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::VOLATILITY]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY
		
	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," MATURITY not positive"); //		positivity of MATURITY
	  
	 if (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::BARRIERSTARTTIME]<=0) return ArgumentList_Checking_Result(false," BARRIERSTARTTIME not positive"); //		positivity of BARRIERSTARTTIME

	 if ((a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_IN) &&
		 (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_OUT) &&
		 (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_IN) &&
		 (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_OUT) &&
		 (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_IN) &&
		 (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::OPTIONTYPE]!=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_OUT))
		  return ArgumentList_Checking_Result(false," OPTIONTYPE should be CROSS_AND_IN,CROSS_AND_OUT,INSIDE_UP_AND_IN,INSIDE_UP_AND_OUT,INSIDE_DOWN_AND_IN,INSIDE_DOWN_AND_OUT"); //	checking of option type

	 if ((a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CALLPUT]!=K_CALL) &&
		  (a[ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CALLPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLPUT should be K_CALL or K_PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::Nb_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}


//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_BS_SingleBarrier_2Asset_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_BS_SingleBarrier_2Asset_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BS_SingleBarrier_2Asset_Formula::Nb_Parameters)
	 {
		 throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_SingleBarrier_2Asset_Formula::value  : bad argsize");
	 }
	 double val=TwoAsset_Single_Barrier_Option(
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::FORWARD1],	//FORWARD1
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::FORWARD2],		//FORWARD2
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::STRIKE1],	//STRIKE1
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::STRIKE2],		//STRIKE2
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::MATURITY],	//MATURITY
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::VOLATILITY1],	//VOLATILITY1
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::VOLATILITY2],	//VOLATILITY2
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::CORRELATION],	//CORRELATION
			 0.0,												//discountRate 
			 0.0,												//dividend1
			 0.0,												//dividend2
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::CALLPUT],	//CALLPUT
			 a[ARM_CF_BS_SingleBarrier_2Asset_Formula::OPTIONTYPE]	//OPTIONTYPE
			);
		return val;
	}

double ARM_CF_BS_SingleBarrier_2Asset_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_BS_SingleBarrier_2Asset_Formula::value,i,a,s);

}


double ARM_CF_BS_SingleBarrier_2Asset_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_BS_SingleBarrier_2Asset_Formula::value,i,j,a,s1,s2);

}





 double ARM_CF_BS_SingleBarrier_2Asset_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::FORWARD1 :
			return 0.00001;		// FORWARD1
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::STRIKE1 :
			return 0.00001;		// STRIKE1
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::FORWARD2 :
			return 0.00001;		// FORWARD2
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::STRIKE2 :
			return 0.00001;		// STRIKE2
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::VOLATILITY1 :
			return 0.00001;		// VOLATILITY1
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::VOLATILITY2 :
			return 0.00001;		// VOLATILITY2
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::CORRELATION :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_BS_SingleBarrier_2Asset_Formula::MATURITY :
			return 0.00001;		// MATURITY
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_SingleBarrier_2Asset_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_BS_SingleBarrier_2Asset_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_BS_SingleBarrier_2Asset_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 
	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::FORWARD1]<=0) return ArgumentList_Checking_Result(false," FORWARD1 not positive"); //		positivity of FORWARD

	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::STRIKE1]<=0) return ArgumentList_Checking_Result(false," STRIKE1 not positive"); //			positivity of STRIKE
	 
	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::FORWARD2]<=0) return ArgumentList_Checking_Result(false," FORWARD2 not positive"); //			positivity of BARRIERE
	
	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::STRIKE2 ]<=0) return ArgumentList_Checking_Result(false," STRIKE2 not positive"); //		positivity of REBATE
	 	
	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::VOLATILITY1]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY
		
	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of MATURITY
	  
	 if (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," MATURITY not positive"); //		positivity of BARRIERSTARTTIME

	 if (abs(a[ARM_CF_BS_SingleBarrier_2Asset_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION) not <=1"); //		Abs(CORRELATION)<=1


	 if ((a[ARM_CF_BS_SingleBarrier_2Asset_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarrier_2Asset_Formula::DOWN_AND_IN) &&
		 (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarrier_2Asset_Formula::UP_AND_IN) &&
		 (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarrier_2Asset_Formula::DOWN_AND_OUT) &&
		 (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::OPTIONTYPE]!=ARM_CF_BS_SingleBarrier_2Asset_Formula::UP_AND_OUT) )
		  return ArgumentList_Checking_Result(false," OPTIONTYPE should be DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"); //	checking of option type

	 if ((a[ARM_CF_BS_SingleBarrier_2Asset_Formula::CALLPUT]!=K_CALL) &&
		  (a[ARM_CF_BS_SingleBarrier_2Asset_Formula::CALLPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLPUT should be K_CALL or K_PUT"); //	

 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_BS_SingleBarrier_2Asset_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>ARM_CF_BS_SingleBarrier_2Asset_Formula::Nb_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




CC_END_NAMESPACE()
 
#undef ARM_CF_SQRT_2_PI

