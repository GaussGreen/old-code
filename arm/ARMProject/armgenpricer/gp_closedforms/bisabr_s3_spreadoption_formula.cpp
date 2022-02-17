/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_s3_spreadoption_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date November 2006
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/spreadoption_bisabr_interface.h"
#include "gpclosedforms/bisabr_s3_spreadoption_formula.h"
#include "gpclosedforms/powerspreadoption.h"
#include "gpclosedforms/sabrvanilla.h"

#include "expt.h"

CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
///      Value of the formula 
///
//////////////////////////////////////////////////////////////////////////////////////////




double ARM_CF_BiSABR_S3_SpreadOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_BiSABR_S3_SpreadOption_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value  : bad argsize");
	 }

	 ///		Here the limitation of the dependence in the smile and the copula is that they 
	 ///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	 ///


	 double val=structure::PowerSpreadOption_GaussianPricing(
		 /// ****************  spread **************************
		  ArgumentList(
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX1],	//INDEX1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA1],	//ALPHA1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA1],	//BETA1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO1],	//RHO1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NU1], 	//NU1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX2],	//INDEX2
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA2],	//ALPHA2
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA2],	//BETA2
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO2],	//RHO2
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NU2],	    //NU2
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOS],	//RHOS
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOV],	//RHOV
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOC12],	//RHOC12
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOC21],	//RHOC21
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::FLAG]	    //FLAG
		 ),
		 ///   *************  index ****************************
		 ArgumentList(
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX3],	//INDEX3
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA3],	//ALPHA3
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA3],	//BETA3
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO3],	//RHO3
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NU3],	    //NU3
		 ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL,		// a choice 
		a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NBSTEPS]	//NBSTEPS											// a choice, but it is not necessary
		 ),
		/// ****************** correlation ******************* 
		
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::CORRELATION]	//CORRELATION
		 ,
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::T],			//TIMETOEXPIRATION

		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NBSTEPS],		//NBSTEPS
		/// *******************  descriptionof the payoff *****
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::A1],	//A1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::B1],	//B1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::K1],	//K1
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::A2],	//A2
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::B2],	//B2 
		 a[ARM_CF_BiSABR_S3_SpreadOption_Formula::K2]	//K2
		 );
		 
	 return val;
	 
}


////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_S3_SpreadOption_Formula::value(int i,const ArgumentList& a, double s)
{
	 		return centred_standard_first_derivative(ARM_CF_BiSABR_S3_SpreadOption_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_S3_SpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_BiSABR_S3_SpreadOption_Formula::value,i,j,a,s1,s2);
}


////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_BiSABR_S3_SpreadOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX1 :
		return 0.00001;		// INDEX1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA1 :
		return 0.00001;		// ALPHA1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA1 :
		return 0.00001;		// BETA1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO1 :
		return 0.00001;		// RHO1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::NU1 :
		return 0.00001;		// NU1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX2 :
		return 0.00001;		// INDEX2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA2 :
		return 0.00001;		// ALPHA2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA2 :
		return 0.00001;		// BETA2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO2 :
		return 0.00001;		// RHO2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::NU2 :
		return 0.00001;		// NU2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOS :
		return 0.00001;		// RHOS
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOV :
		return 0.00001;		// RHOV
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOC12 :
		return 0.00001;		// RHOC12
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOC21 :
		return 0.00001;		// RHOC21
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX3 :
		return 0.00001;		// INDEX3
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA3 :
		return 0.00001;		// ALPHA3
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA3 :
		return 0.00001;		// BETA3
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO3 :
		return 0.00001;		// RHO3
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::NU3 :
		return 0.00001;		// NU3
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::CORRELATION :
		return 0.00001;		// CORRELATION
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::A1 :
		return 0.00001;		// A1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::B1 :
		return 0.00001;		// B1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::K1 :
		return 0.00001;		// K1
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::A2 :
		return 0.00001;		// A2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::B2 :
		return 0.00001;		// B2
	case ARM_CF_BiSABR_S3_SpreadOption_Formula::K2 :
		return 0.00001;		// K2

	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BiSABR_S3_SpreadOption_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_BiSABR_S3_SpreadOption_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_BiSABR_S3_SpreadOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	}
	
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1  not positive"); //			positivity of INDEX1
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX2]<0) return ArgumentList_Checking_Result(false," INDEX2  not positive"); //			positivity of INDEX2
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA1]<=0) return ArgumentList_Checking_Result(false," ALPHA1   not positive"); //			positivity of ALPHA1
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA1]<=0) return ArgumentList_Checking_Result(false," BETA1   not positive"); //			positivity of BETA1
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO1])>1.) return ArgumentList_Checking_Result(false," RHO1     not <1 and >-1"); //	Abs(RHO1)<=1
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NU1]<=0) return ArgumentList_Checking_Result(false," NU1      not positive"); //			positivity of NU1
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA2]<0) return ArgumentList_Checking_Result(false," ALPHA2    not positive"); //			positivity of ALPHA2
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA2]<=0) return ArgumentList_Checking_Result(false," BETA2   not positive"); //			positivity of BETA2
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO2])>1.) return ArgumentList_Checking_Result(false," RHO2    not <1 and >-1"); //		Abs(RHO2)<=1 
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NU2]<=0) return ArgumentList_Checking_Result(false," NU2      not positive"); //			positivity of NU2
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOS])>1.) return ArgumentList_Checking_Result(false," RHOS    not <1 and >-1"); //		Abs(RHOS)<=1 
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOV])>1.) return ArgumentList_Checking_Result(false," RHOV    not <1 and >-1"); //		Abs(RHOV)<=1 
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOC12])>1.) return ArgumentList_Checking_Result(false," RHOC12    not <1 and >-1"); //		Abs(RHOC12)<=1 
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHOC21])>1.) return ArgumentList_Checking_Result(false," RHOC21    not <1 and >-1"); //		Abs(RHOC21)<=1 
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::INDEX3]<=0) return ArgumentList_Checking_Result(false," INDEX3  not positive"); //			positivity of INDEX3
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::ALPHA3]<=0) return ArgumentList_Checking_Result(false," ALPHA3   not positive"); //			positivity of ALPHA3
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::BETA3]<=0) return ArgumentList_Checking_Result(false," BETA3   not positive"); //			positivity of BETA3
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::RHO3])>1.) return ArgumentList_Checking_Result(false," RHO3     not <1 and >-1"); //	Abs(RHO3)<=1
	if(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::NU3]<=0) return ArgumentList_Checking_Result(false," NU3      not positive"); //			positivity of NU3
	if(abs(a[ARM_CF_BiSABR_S3_SpreadOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," CORRELATION     not <1 and >-1"); //	Abs(CORRELATION)<=1



	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_BiSABR_S3_SpreadOption_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///		so the number of parameters is 19
	if ((rank<0)||(rank>=ARM_CF_BiSABR_S3_SpreadOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
