/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file shifted_lognormal_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/vanilla_shifted_lognormal.h"
#include "gpclosedforms/shifted_lognormal_formula.h"

#include <glob/expt.h>

CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////
///
///      Value of the formula 
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_Shifted_Lognormal_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_Shifted_Lognormal_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Shifted_Lognormal_Formula::value  : bad argsize");
	 }
	 int callput=a[ARM_CF_Shifted_Lognormal_Formula::CALLORPUT];

switch (callput)
	{
	case K_CALL :
		{
			return shifted_lognormal_vanilla_call(
				a[ARM_CF_Shifted_Lognormal_Formula::INDEX],				//INDEX			
				a[ARM_CF_Shifted_Lognormal_Formula::STRIKE],				//STRIKE			
				a[ARM_CF_Shifted_Lognormal_Formula::TIMETOMATURITY],		//TIMETOMATURITY	
				a[ARM_CF_Shifted_Lognormal_Formula::VOLATILITY],			//VOLATILITY		
				a[ARM_CF_Shifted_Lognormal_Formula::ALPHA]	//ALPHA				
				);
			break;
		}
	case K_PUT :
		{
			return shifted_lognormal_vanilla_call(
				a[ARM_CF_Shifted_Lognormal_Formula::INDEX],				//INDEX			
				a[ARM_CF_Shifted_Lognormal_Formula::STRIKE],				//STRIKE			
				a[ARM_CF_Shifted_Lognormal_Formula::TIMETOMATURITY],		//TIMETOMATURITY	
				a[ARM_CF_Shifted_Lognormal_Formula::VOLATILITY],			//VOLATILITY		
				a[ARM_CF_Shifted_Lognormal_Formula::ALPHA]	//ALPHA			
				)-(a[ARM_CF_Shifted_Lognormal_Formula::INDEX]-a[ARM_CF_Shifted_Lognormal_Formula::STRIKE]);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Shifted_Lognormal_Formula : callput , bad input :");
			break;
		}
	}
}



////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Shifted_Lognormal_Formula::value(int i,const ArgumentList& a, double s)
{
	 		return standard_first_derivative(ARM_CF_Shifted_Lognormal_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Shifted_Lognormal_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_Shifted_Lognormal_Formula::value,i,j,a,s1,s2);
}


////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Shifted_Lognormal_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_Shifted_Lognormal_Formula::INDEX :
		return 0.00001;		// INDEX
	case ARM_CF_Shifted_Lognormal_Formula::STRIKE :
		return 0.00001;		// STRIKE
	case ARM_CF_Shifted_Lognormal_Formula::TIMETOMATURITY :
		return 0.00001;		// TIMETOMATURITY
	case ARM_CF_Shifted_Lognormal_Formula::VOLATILITY :
		return 0.00001;		// VOLATILITY
	case ARM_CF_Shifted_Lognormal_Formula::ALPHA :  
		return 0.00001;		// ALPHA
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Shifted_Lognormal_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_Shifted_Lognormal_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_Shifted_Lognormal_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	}
	
	if(a[ARM_CF_Shifted_Lognormal_Formula::INDEX]<0) return ArgumentList_Checking_Result(false," INDEX  not positive");					//			positivity of INDEX
	if(a[ARM_CF_Shifted_Lognormal_Formula::STRIKE]<0) return ArgumentList_Checking_Result(false," STRIKE  not positive");					//			positivity of STRIKE
	if(a[ARM_CF_Shifted_Lognormal_Formula::TIMETOMATURITY]<0) return ArgumentList_Checking_Result(false," TIMETOMATURITY   not positive");		//			positivity of TIMETOMATURITY
	if(a[ARM_CF_Shifted_Lognormal_Formula::VOLATILITY]<0) return ArgumentList_Checking_Result(false," VOLATILITY   not positive or zero");	//			positivity of VOLATILITY	
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_Shifted_Lognormal_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///		so the number of parameters is 19
	if ((rank<0)||(rank>=ARM_CF_Shifted_Lognormal_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
