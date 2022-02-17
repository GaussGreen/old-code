/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_spreadoption_formula.cpp
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
#include "gpclosedforms/spreadoption_bisabr_interface.h"
#include "gpclosedforms/bisabr_spreadoption_formula.h"
#include "gpclosedforms/bisabr_spreadoption.h"


#include <glob/expt.h>

CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
///      Value of the formula 
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_SpreadOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BiSABR_SpreadOption_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BiSABR_SpreadOption_Formula::value  : bad argsize");
	 }
	 int callput=a[ARM_CF_BiSABR_SpreadOption_Formula::CALLORPUT];
	 
	 
	 return Packaged_BiSABR_SpreadOption(
		 a[ARM_CF_BiSABR_SpreadOption_Formula::INDEX1],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::ALPHA1],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::BETA1],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::RHO1],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::NU1],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::INDEX2],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::ALPHA2],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::BETA2],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::RHO2],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::NU2],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::K],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::T],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::CALLORPUT],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::RHOS],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::RHOV],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::RHOC12],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::RHOC21],
		 a[ARM_CF_BiSABR_SpreadOption_Formula::FLAG]
		 );
	 
	 
	 
	 
}



////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_SpreadOption_Formula::value(int i,const ArgumentList& a, double s)
{
	 		return centred_standard_first_derivative(ARM_CF_BiSABR_SpreadOption_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_SpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_BiSABR_SpreadOption_Formula::value,i,j,a,s1,s2);
}


////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_SpreadOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_BiSABR_SpreadOption_Formula::INDEX1 :
		return 0.00001;		// INDEX1
	case ARM_CF_BiSABR_SpreadOption_Formula::ALPHA1 :
		return 0.00001;		// ALPHA1
	case ARM_CF_BiSABR_SpreadOption_Formula::BETA1 :
		return 0.00001;		// BETA1
	case ARM_CF_BiSABR_SpreadOption_Formula::RHO1 :
		return 0.00001;		// RHO1
	case ARM_CF_BiSABR_SpreadOption_Formula::NU1 :
		return 0.00001;		// NU1
	case ARM_CF_BiSABR_SpreadOption_Formula::INDEX2 :
		return 0.00001;		// INDEX2
	case ARM_CF_BiSABR_SpreadOption_Formula::ALPHA2 :
		return 0.00001;		// ALPHA2
	case ARM_CF_BiSABR_SpreadOption_Formula::BETA2 :
		return 0.00001;		// BETA2
	case ARM_CF_BiSABR_SpreadOption_Formula::RHO2 :
		return 0.00001;		// RHO2
	case ARM_CF_BiSABR_SpreadOption_Formula::NU2 :
		return 0.00001;		// NU2
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOS :
		return 0.00001;		// RHOS
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOV :
		return 0.00001;		// RHOV
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOC12 :
		return 0.00001;		// RHOC12
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOC21 :
		return 0.00001;		// RHOC21
	case ARM_CF_BiSABR_SpreadOption_Formula::K :
		return 0.00001;		// K
	case ARM_CF_BiSABR_SpreadOption_Formula::T :
		return 0.00001;		// T

	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BiSABR_SpreadOption_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_BiSABR_SpreadOption_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_BiSABR_SpreadOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	}
	
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1  not positive"); //			positivity of INDEX1
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::INDEX2]<0) return ArgumentList_Checking_Result(false," INDEX2  not positive"); //			positivity of INDEX2
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::ALPHA1]<=0) return ArgumentList_Checking_Result(false," ALPHA1   not positive"); //			positivity of ALPHA1
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::BETA1]<=0) return ArgumentList_Checking_Result(false," BETA1   not positive"); //			positivity of BETA1
	if(abs(a[ARM_CF_BiSABR_SpreadOption_Formula::RHO1])>1.) return ArgumentList_Checking_Result(false," RHO1     not <1 and >-1"); //	Abs(RHO1)<=1
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::NU1]<=0) return ArgumentList_Checking_Result(false," NU1      not positive"); //			positivity of NU1
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::ALPHA2]<0) return ArgumentList_Checking_Result(false," ALPHA2    not positive"); //			positivity of ALPHA2
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::BETA2]<=0) return ArgumentList_Checking_Result(false," BETA2   not positive"); //			positivity of BETA2
	if(abs(a[ARM_CF_BiSABR_SpreadOption_Formula::RHO2])>1.) return ArgumentList_Checking_Result(false," RHO2    not <1 and >-1"); //		Abs(RHO2)<=1 
	if(a[ARM_CF_BiSABR_SpreadOption_Formula::NU2]<=0) return ArgumentList_Checking_Result(false," NU2      not positive"); //			positivity of NU2
	if(abs(a[ARM_CF_BiSABR_SpreadOption_Formula::RHOS])>1.) return ArgumentList_Checking_Result(false," RHOS    not <1 and >-1"); //		Abs(RHOS)<=1 
	if(abs(a[ARM_CF_BiSABR_SpreadOption_Formula::RHOV])>1.) return ArgumentList_Checking_Result(false," RHOV    not <1 and >-1"); //		Abs(RHOV)<=1 
	if(abs(a[ARM_CF_BiSABR_SpreadOption_Formula::RHOC12])>1.) return ArgumentList_Checking_Result(false," RHOC12    not <1 and >-1"); //		Abs(RHOC12)<=1 
	if(abs(a[ARM_CF_BiSABR_SpreadOption_Formula::RHOC21])>1.) return ArgumentList_Checking_Result(false," RHOC21    not <1 and >-1"); //		Abs(RHOC21)<=1 


	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_BiSABR_SpreadOption_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///		so the number of parameters is 19
	if ((rank<0)||(rank>=ARM_CF_BiSABR_SpreadOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
