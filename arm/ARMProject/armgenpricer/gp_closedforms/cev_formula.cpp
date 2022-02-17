/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file SABR_Analytics.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "firsttoinc.h"
#include "gpbase/port.h"
#include <cmath>

#include "gpclosedforms/cev_formula.h"
#include "gpclosedforms/cev.h"

#include "expt.h"



CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_CEV_VanillaOption_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_CEV_VanillaOption_Formula::value(const ArgumentList& a)
{
	int argsize=a.size();
	if (argsize!=ARM_CF_CEV_VanillaOption_Formula::Nb_Parameters) 
	{
		throw ("ARM_CF_CEV_VanillaOption_Formula::value  : bad argsize");
	}
	int callput=(int) a[ARM_CF_CEV_VanillaOption_Formula::CALLORPUT];
	double beta=a[ARM_CF_CEV_VanillaOption_Formula::BETA];
	double val;
	switch (callput)
	{
	case K_CALL :
		{
			
			val= CEV_Call(
				a[ARM_CF_CEV_VanillaOption_Formula::INDEX],	//INDEX
				a[ARM_CF_CEV_VanillaOption_Formula::STRIKE],	//STRIKE
				a[ARM_CF_CEV_VanillaOption_Formula::TIMETORESET],	//TIMETOMATURITY
				a[ARM_CF_CEV_VanillaOption_Formula::MU],	//MU
				a[ARM_CF_CEV_VanillaOption_Formula::VOLATILITY],	//VOLATILITY
				-1/(2*(beta-1)),										//NU
				a[ARM_CF_CEV_VanillaOption_Formula::NBTERMS]
				);
			break;
		}
	case K_PUT :
		{
			val= CEV_Call(
				a[ARM_CF_CEV_VanillaOption_Formula::INDEX],	//INDEX
				a[ARM_CF_CEV_VanillaOption_Formula::STRIKE],	//STRIKE
				a[ARM_CF_CEV_VanillaOption_Formula::TIMETORESET],	//TIMETOMATURITY
				a[ARM_CF_CEV_VanillaOption_Formula::MU],	//MU
				a[ARM_CF_CEV_VanillaOption_Formula::VOLATILITY],	//VOLATILITY
				-1/(2*(beta-1)),	//BETA
				a[ARM_CF_CEV_VanillaOption_Formula::NBTERMS]
				)-(a[ARM_CF_CEV_VanillaOption_Formula::INDEX]-a[ARM_CF_CEV_VanillaOption_Formula::STRIKE]);
			
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_VanillaOption_Formula : callput , bad input :");
			break;
		}
	}
	return val;
}




double ARM_CF_CEV_VanillaOption_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_CEV_VanillaOption_Formula::value,i,a,s);
}


double ARM_CF_CEV_VanillaOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_CEV_VanillaOption_Formula::value,i,j,a,s1,s2);

}



 double ARM_CF_CEV_VanillaOption_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_CEV_VanillaOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_CEV_VanillaOption_Formula::INDEX :
			return 0.00001;		// Forward index 
		case ARM_CF_CEV_VanillaOption_Formula::STRIKE :
			return 0.00001;		// Strike
		case ARM_CF_CEV_VanillaOption_Formula::VOLATILITY :
			return 0.00001;		// Volatility
		case ARM_CF_CEV_VanillaOption_Formula::BETA :
			return 0.00001;		// Beta
		case ARM_CF_CEV_VanillaOption_Formula::MU :
			return 0.00001;		// Mu
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_VanillaOption_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_CEV_VanillaOption_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_CEV_VanillaOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_CEV_VanillaOption_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET
	 if (a[ARM_CF_CEV_VanillaOption_Formula::INDEX]<=0) return ArgumentList_Checking_Result(false," INDEX not positive"); //			positivity of INDEX
	 if (a[ARM_CF_CEV_VanillaOption_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	
	 if (a[ARM_CF_CEV_VanillaOption_Formula::VOLATILITY ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY	
///	 if (a[ARM_CF_CEV_VanillaOption_Formula::BETA]>1) return ArgumentList_Checking_Result(false," BETA not <1"); //		negativity of BETA

	 if (abs(a[ARM_CF_CEV_VanillaOption_Formula::NBTERMS])<=1) return ArgumentList_Checking_Result(false," NBTERMS should be >= 1"); //		NBTERMS should >=1
	
	if ((a[ARM_CF_CEV_VanillaOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_CEV_VanillaOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_CEV_VanillaOption_Formula:CALLORPUT should be CALL or PUT"); //			positivity of STRIKE
 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_CEV_VanillaOption_Formula::check_dimension(int rank)
{
	 if ((rank<0)||(rank>=ARM_CF_CEV_VanillaOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_CEV_DoubleBarrierOption_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////

double ARM_CF_CEV_DoubleBarrierOption_Formula::value(const ArgumentList& a)
{
	int argsize=a.size();
	if (argsize!=ARM_CF_CEV_DoubleBarrierOption_Formula::Nb_Parameters) 
	{
		throw ("ARM_CF_CEV_DoubleBarrierOption_Formula::value  : bad argsize");
	}
	int callput=(int) a[ARM_CF_CEV_DoubleBarrierOption_Formula::CALLORPUT];
	double val;
	switch (callput)
	{
	case K_CALL :
		{
			
			val= CEV_DoubleBarrierCall(
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::INDEX],	//INDEX
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::STRIKE],	//STRIKE
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
				0,
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::BETA]-1.,	//BETA
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::MU],	//MU
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::VOLATILITY],	//VOLATILITY
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::BARRIERDOWN],	//BARRIERDOWN
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::BARRIERUP],	//BARRIERUP
				a[ARM_CF_CEV_DoubleBarrierOption_Formula::NBTERMS]
				);
			break;
		}
	case K_PUT :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_DoubleBarrierOption_Formula : callput , Put: Not Yet Implemented :");
			
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_DoubleBarrierOption_Formula : callput , bad input :");
			break;
		}
	}
	return val;
}




double ARM_CF_CEV_DoubleBarrierOption_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_CEV_DoubleBarrierOption_Formula::value,i,a,s);
}


double ARM_CF_CEV_DoubleBarrierOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_CEV_DoubleBarrierOption_Formula::value,i,j,a,s1,s2);

}



 double ARM_CF_CEV_DoubleBarrierOption_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_CEV_DoubleBarrierOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_CEV_DoubleBarrierOption_Formula::INDEX :
			return 0.00001;		// Forward index 
		case ARM_CF_CEV_DoubleBarrierOption_Formula::STRIKE :
			return 0.00001;		// Strike
		case ARM_CF_CEV_DoubleBarrierOption_Formula::VOLATILITY :
			return 0.00001;		// Volatility
		case ARM_CF_CEV_DoubleBarrierOption_Formula::BETA :
			return 0.00001;		// Beta
		case ARM_CF_CEV_DoubleBarrierOption_Formula::MU :
			return 0.00001;		// Mu
		case ARM_CF_CEV_DoubleBarrierOption_Formula::BARRIERDOWN :
			return 0.00001;		// Barrier down
		case ARM_CF_CEV_DoubleBarrierOption_Formula::BARRIERUP :
			return 0.00001;		// Barrier up
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_DoubleBarrierOption_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_CEV_DoubleBarrierOption_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_CEV_DoubleBarrierOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET
	 if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::INDEX]<=0) return ArgumentList_Checking_Result(false," INDEX not positive"); //			positivity of INDEX
	 if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	
	 if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::VOLATILITY ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY	
	 if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::BETA]<0) return ArgumentList_Checking_Result(false," BETA not positive"); //		negativity of BETA

	  if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::BARRIERDOWN]<=0) return ArgumentList_Checking_Result(false," BARRIERDOWN not positive"); //		positivity of BARRIERDOWN
	   if (a[ARM_CF_CEV_DoubleBarrierOption_Formula::BARRIERUP]<=0) return ArgumentList_Checking_Result(false," BARRIERUP not positive"); //		positivity of BARRIERUP
	 
	 if (abs(a[ARM_CF_CEV_DoubleBarrierOption_Formula::NBTERMS])<=1) return ArgumentList_Checking_Result(false," NBTERMS should be >= 1"); //		NBTERMS should >=1
	
	if ((a[ARM_CF_CEV_DoubleBarrierOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_CEV_DoubleBarrierOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_CEV_DoubleBarrierOption_Formula:CALLORPUT should be CALL or PUT"); //			positivity of STRIKE
 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_CEV_DoubleBarrierOption_Formula::check_dimension(int rank)
{
	 if ((rank<0)||(rank>=ARM_CF_CEV_DoubleBarrierOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

 //////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_CEV_BarrierOption_Formula
///
//////////////////////////////////////////////////////////////////////////////////////////
double CEV_SingleBarrierUpandOutCall(double S,double K,double T,double r,double beta,double mu,double delta,double U,int  n);

double ARM_CF_CEV_BarrierOption_Formula::value(const ArgumentList& a)
{
	int argsize=a.size();
	if (argsize!=ARM_CF_CEV_BarrierOption_Formula::Nb_Parameters) 
	{
		throw ("ARM_CF_CEV_BarrierOption_Formula::value  : bad argsize");
	}
	double val;
	int callput=(int) a[ARM_CF_CEV_BarrierOption_Formula::CALLORPUT];
	switch (callput)
	{
	case K_CALL :
		{
			if(a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]==ARM_CF_CEV_BarrierOption_Formula::UP_AND_OUT)
			{
				
				val= CEV_SingleBarrierUpandOutCall(
					a[ARM_CF_CEV_BarrierOption_Formula::INDEX],	//INDEX
					a[ARM_CF_CEV_BarrierOption_Formula::STRIKE],	//STRIKE
					a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
					0,
					a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.,	//BETA
					a[ARM_CF_CEV_BarrierOption_Formula::MU],	//MU
					a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY],	//VOLATILITY
					a[ARM_CF_CEV_BarrierOption_Formula::BARRIER],	//BARRIER
					a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS]
					);
			}
			else
				
				if(a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]==ARM_CF_CEV_BarrierOption_Formula::UP_AND_IN)
				{
					
					double val0=CEV_Call(
						a[ARM_CF_CEV_BarrierOption_Formula::INDEX],	//INDEX
						a[ARM_CF_CEV_BarrierOption_Formula::STRIKE],	//STRIKE
						a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
						a[ARM_CF_CEV_BarrierOption_Formula::MU],	//MU
						a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY]*
						pow(	a[ARM_CF_CEV_BarrierOption_Formula::INDEX],
								a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.),	//VOLATILITY
						-1./(2.*(a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.)),		//NU
						a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS]
						);
					
					val=val0 - CEV_SingleBarrierUpandOutCall(
						a[ARM_CF_CEV_BarrierOption_Formula::INDEX],	//INDEX
						a[ARM_CF_CEV_BarrierOption_Formula::STRIKE],	//STRIKE
						a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
						0,
						a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.,	//BETA
						a[ARM_CF_CEV_BarrierOption_Formula::MU],	//MU
						a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY],	//VOLATILITY
						a[ARM_CF_CEV_BarrierOption_Formula::BARRIER],	//BARRIER
						a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS]
						);
				}
				
				else
				{	
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_BarrierOption_Formula : optiontype , Not Yet Implemented or bad input :");
					break;
				}
			break;
			
		}
	case K_PUT :
		{
		if(a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]==ARM_CF_CEV_BarrierOption_Formula::DOWN_AND_OUT)
			{
				
				val= CEV_SingleBarrierDownandOutPut(
					a[ARM_CF_CEV_BarrierOption_Formula::INDEX],	//INDEX
					a[ARM_CF_CEV_BarrierOption_Formula::STRIKE],	//STRIKE
					a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
					0,
					a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.,	//BETA
					a[ARM_CF_CEV_BarrierOption_Formula::MU],	//MU
					a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY],	//VOLATILITY
					a[ARM_CF_CEV_BarrierOption_Formula::BARRIER],	//BARRIER
					a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS]
					);
			}
			else
				
				if(a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]==ARM_CF_CEV_BarrierOption_Formula::DOWN_AND_IN)
				{
					
					double val0=CEV_Call(
						a[ARM_CF_CEV_BarrierOption_Formula::INDEX],	//INDEX
						a[ARM_CF_CEV_BarrierOption_Formula::STRIKE],	//STRIKE
						a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
						a[ARM_CF_CEV_BarrierOption_Formula::MU],	//MU
						a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY]*
						pow(	a[ARM_CF_CEV_BarrierOption_Formula::INDEX],
								a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.),	//VOLATILITY
						-1./(2.*(a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.)),										//NU
						a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS]
						)-(a[ARM_CF_CEV_VanillaOption_Formula::INDEX]-a[ARM_CF_CEV_VanillaOption_Formula::STRIKE]);
					
					val=val0 - CEV_SingleBarrierDownandOutPut(
						a[ARM_CF_CEV_BarrierOption_Formula::INDEX],	//INDEX
						a[ARM_CF_CEV_BarrierOption_Formula::STRIKE],	//STRIKE
						a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET],	//TIMETOMATURITY
						0,
						a[ARM_CF_CEV_BarrierOption_Formula::BETA]-1.,	//BETA
						a[ARM_CF_CEV_BarrierOption_Formula::MU],	//MU
						a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY],	//VOLATILITY
						a[ARM_CF_CEV_BarrierOption_Formula::BARRIER],	//BARRIER
						a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS]
						);
				}
				
				else
				{	
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_BarrierOption_Formula : optiontype , Not Yet Implemented or bad input :");
					break;
				}
			break;
			
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_BarrierOption_Formula : callput , bad input :");
			break;
		}
	}
	return val;
}




double ARM_CF_CEV_BarrierOption_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_CEV_BarrierOption_Formula::value,i,a,s);
}


double ARM_CF_CEV_BarrierOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_CEV_BarrierOption_Formula::value,i,j,a,s1,s2);

}



 double ARM_CF_CEV_BarrierOption_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_CEV_BarrierOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_CEV_BarrierOption_Formula::INDEX :
			return 0.00001;		// Forward index 
		case ARM_CF_CEV_BarrierOption_Formula::STRIKE :
			return 0.00001;		// Strike
		case ARM_CF_CEV_BarrierOption_Formula::VOLATILITY :
			return 0.00001;		// Volatility
		case ARM_CF_CEV_BarrierOption_Formula::BETA :
			return 0.00001;		// Beta
		case ARM_CF_CEV_BarrierOption_Formula::MU :
			return 0.00001;		// Mu
		case ARM_CF_CEV_BarrierOption_Formula::BARRIER :
			return 0.00001;		// Barrier 
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_CEV_BarrierOption_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_CEV_BarrierOption_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_CEV_BarrierOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_CEV_BarrierOption_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET
	 if (a[ARM_CF_CEV_BarrierOption_Formula::INDEX]<=0) return ArgumentList_Checking_Result(false," INDEX not positive"); //			positivity of INDEX
	 if (a[ARM_CF_CEV_BarrierOption_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //			positivity of STRIKE
	
	 if (a[ARM_CF_CEV_BarrierOption_Formula::VOLATILITY ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY not positive"); //		positivity of VOLATILITY	
	 if (a[ARM_CF_CEV_BarrierOption_Formula::BETA]<0) return ArgumentList_Checking_Result(false," BETA not positive"); //		negativity of BETA

	  if (a[ARM_CF_CEV_BarrierOption_Formula::BARRIER]<=0) return ArgumentList_Checking_Result(false," BARRIER not positive"); //		positivity of BARRIER
	 
	 if (abs(a[ARM_CF_CEV_BarrierOption_Formula::NBTERMS])<=1) return ArgumentList_Checking_Result(false," NBTERMS should be >= 1"); //		NBTERMS should >=1
	
	if ((a[ARM_CF_CEV_BarrierOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_CEV_BarrierOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_CEV_BarrierOption_Formula:CALLORPUT should be CALL or PUT"); //			Call or Put

	 if ((a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]!=ARM_CF_CEV_BarrierOption_Formula::DOWN_AND_IN) &&
		 (a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]!=ARM_CF_CEV_BarrierOption_Formula::UP_AND_IN) &&
		 (a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]!=ARM_CF_CEV_BarrierOption_Formula::DOWN_AND_OUT) &&
		 (a[ARM_CF_CEV_BarrierOption_Formula::OPTIONTYPE]!=ARM_CF_CEV_BarrierOption_Formula::UP_AND_OUT))
		  return ArgumentList_Checking_Result(false," OPTIONTYPE should be DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT or UP_AND_OUT"); //	checking of option type
							//		Up or Down
 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_CEV_BarrierOption_Formula::check_dimension(int rank)
{
	 if ((rank<0)||(rank>=ARM_CF_CEV_BarrierOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}


CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/