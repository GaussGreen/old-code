/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_shiftedlognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_shiftedlognormal_formula.h"
#include "gpclosedforms/spreadoption_sabr_formula.h"

#include "expt.h"

/// uses the template PowerSpreadOption_Pricing_With_Limits which declared in powerspreadoption.h



CC_BEGIN_NAMESPACE(ARM)



//////////////////////////////////////////////////////////////////////////////////////////
///
///      Value of the formula 
///
//////////////////////////////////////////////////////////////////////////////////////////


double ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value(const ArgumentList& a)
	{
	 int nbsteps;
	 int argsize=a.size();
	 if (argsize!=ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value  : bad argsize");
	 }

	 ///		Here the limitation of the dependence in the smile and the copula is that they 
	 ///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	 ///


	 double val=structure::PowerSpreadOption_Pricing_With_Limits(
		 ArgumentList(
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX1],	//INDEX1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA1],	//SIGMA1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA1]	//ALPHA1
		 
		 ),
		 ArgumentList(
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX2],	//INDEX2
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA2],	//SIGMA2
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA2]	//ALPHA2
		 
		 ),
		 ArgumentList(
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::CORRELATION]	//CORRELATION
		 ),
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION],	//TIMETOEXPIRATION
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::NBSTEPS],	//NBSTEPS
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::A1],	//A1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::B1],	//B1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::K1],	//K1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::A2],	//A2
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::B2],	//B2 
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::K2],	//K2
		 GENERIC_SPREADOPTION
		 );
		 
	 return val;
}


// cette fonction utilise PowerSpreadOption_pricing qui supose que les payoff sont differentiable ce qui est rarement le cas, 
/// dans le cas non differentiable, une erreur plus ou moins importante se manifeste
double ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula_value2(const ArgumentList& a)
	{
	 int nbsteps;
	 int argsize=a.size();
	 if (argsize<14)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value  : bad argsize");
	 }
	 if (argsize>14)
	 {
		 nbsteps=a[14];
	 }
	 else
	 {
		 nbsteps=64;
	 }					
	 ///		Here the limitation of the dependence in the smile and the copula is that they 
	 ///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	 ///


	 double val=ARM_CF_PowerSpreadOption_Formula<ShiftedLogNormal_Smile,GaussianCopula>::PowerSpreadOption_Pricing(
		 ArgumentList(
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX1],	//INDEX1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA1],	//SIGMA1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA1]	//ALPHA1
		 
		 ),
		 ArgumentList(
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX2],	//INDEX2
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA2],	//SIGMA2
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA2]	//ALPHA2
		 
		 ),
		 ArgumentList(
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::CORRELATION]	//CORRELATION
		 ),
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION],	//TIMETOEXPIRATION
		 nbsteps,
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::A1],	//A1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::B1],	//B1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::K1],	//K1
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::A2],	//A2
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::B2],	//B2 
		 a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::K2],	//K2
		 GENERIC_SPREADOPTION
		 );
		 
	 return val;
}

////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value(int i,const ArgumentList& a, double s)
{
	 		return standard_first_derivative(ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::value,i,j,a,s1,s2);
}


////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX1 :
		return 0.00001;		// INDEX1
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX2 :
		return 0.00001;		// INDEX2
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA1 :
		return 0.00001;		// SIGMA1
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA1 :
		return 0.00001;		// ALPHA1
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA2 :  
		return 0.00001;		// SIGMA2
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA2 :
		return 0.00001;		// ALPHA2
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::CORRELATION :
		return 0.00001;		// CORRELATION
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION :
		return 0.00001;		// TIMETOEXPIRATION
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::A1 :
		return 0.00001;		// A1
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::B1 :
		return 0.00001;		// B1
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::K1 :
		return 0.00001;		// K1
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::A2 :
		return 0.00001;		// A2
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::B2  :
		return 0.00001;		// B2
	case ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::K2 :
		return 0.00001;		// K2
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()<ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments : <7 ");
	}
	
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1  not positive"); //			positivity of INDEX1
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2  not positive"); //			positivity of INDEX2
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA1]<=0) return ArgumentList_Checking_Result(false," SIGMA1   not positive"); //			positivity of ALPHA1
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA1]<0) return ArgumentList_Checking_Result(false," ALPHA1   not positive or zero"); //			positivity of BETA1
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::SIGMA2]<=0) return ArgumentList_Checking_Result(false," SIGMA2    not positive"); //			positivity of ALPHA2
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::ALPHA2]<0) return ArgumentList_Checking_Result(false," ALPHA2   not positive or 0"); //			positivity of BETA2
	if(abs(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," CORRELATION   not within [-1,+1]"); //	Abs(CORRELATION)<=1 
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION]<=0) return ArgumentList_Checking_Result(false," TIMETOEXPIRATION not positive"); //positivity of TIMETOEXPIRATION
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::NBSTEPS]<2) return ArgumentList_Checking_Result(false," NBSTEPS   not>2"); //			NBSTEPS>=2
	if(a[ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::NBSTEPS]>=178) return ArgumentList_Checking_Result(false," NBSTEPS   not<178"); //			NBSTEPS<2

	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///		so the number of parameters is 19
	 if ((rank<0)||(rank>=ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
