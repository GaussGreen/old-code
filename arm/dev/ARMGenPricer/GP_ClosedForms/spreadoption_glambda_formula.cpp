/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_glambda_student_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/glambda_student_formula.h"
#include "gpclosedforms/smile_glambda.h"
#include "gpclosedforms/digital_templated_pricing.h"


#include <glob/expt.h>




CC_BEGIN_NAMESPACE(ARM)


////////////////////////////////////////////////////////////////////////
///
///      Value of the formula (introduces the checking and the ability to automatically derive
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_GLambda_Student_PowerSpreadOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_GLambda_Student_PowerSpreadOption_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_GLambda_Student_PowerSpreadOption_Formula::value  : bad argsize");
	 }
	 ///		Here the limitation of the dependence in the smile and the copula is that they 
	 ///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	 ///
	 ///		we use here only PowerSpreadOption_Pricing and not Basic_Pricing_With_Limits which is more efficient but which involve the computation of 
	 ///		DistributionZaInverseLimit<smile,copula> that is not simple and prone to error in the case of GLambda 
	 ///
	 ArgumentList* CopulaArg=new ArgumentList(
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::CORRELATION],	//CORRELATION
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::DEGRE]		//DEGRE
		 );
	 StudentCopula cop(CopulaArg);
	 double val=structure::PowerSpreadOption_Pricing(
		 ArgumentList(
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L1],	//FIRST_UNDERLYING_L1
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L2],	//FIRST_UNDERLYING_L2
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L3],	//FIRST_UNDERLYING_L3
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L4],	//FIRST_UNDERLYING_L4
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L5]	//FIRST_UNDERLYING_L5
		 ),
		 ArgumentList(
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L1],	//SECOND_UNDERLYING_L1
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L2],	//SECOND_UNDERLYING_L2
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L3],	//SECOND_UNDERLYING_L3
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L4],	//SECOND_UNDERLYING_L4
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L5]	//SECOND_UNDERLYING_L5
		 ),
		 cop,
		 0,
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::NBSTEPS],
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::A1],	//A1
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::B1],	//B1
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::K1],	//K1
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::A2],	//A2
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::B2],	//B2 
		 a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::K2],	//K2
		 GENERIC_SPREADOPTION	//Algorithm
		 );
		 
	 return val;
}


////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_GLambda_Student_PowerSpreadOption_Formula::value(int i,const ArgumentList& a, double s)

{
		return standard_first_derivative(ARM_CF_GLambda_Student_PowerSpreadOption_Formula::value,i,a,s);
}


////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_GLambda_Student_PowerSpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_GLambda_Student_PowerSpreadOption_Formula::value,i,j,a,s1,s2);
}



////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_GLambda_Student_PowerSpreadOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L1 :
		return 0.00001;		// FIRST_UNDERLYING_L1
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L2 :
		return 0.00001;		// FIRST_UNDERLYING_L2
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L3 :
		return 0.00001;		// FIRST_UNDERLYING_L3
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L4 :
		return 0.00001;		// FIRST_UNDERLYING_L4
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L5 :
		return 0.00001;		// FIRST_UNDERLYING_L5
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L1 :
		return 0.00001;		// SECOND_UNDERLYING_L1
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L2 :  
		return 0.00001;		// SECOND_UNDERLYING_L2
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L3 :
		return 0.00001;		// SECOND_UNDERLYING_L3
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L4  :
		return 0.00001;		// SECOND_UNDERLYING_L4
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L5 :
		return 0.00001;		// SECOND_UNDERLYING_L5
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::CORRELATION :
		return 0.00001;		// CORRELATION
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::DEGRE :
		return 0.00001;		// DEGRE
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::A1 :
		return 0.00001;		// A1
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::B1 :
		return 0.00001;		// B1
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::K1 :
		return 0.00001;		// K1
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::A2 :
		return 0.00001;		// A2
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::B2  :
		return 0.00001;		// B2
	case ARM_CF_GLambda_Student_PowerSpreadOption_Formula::K2 :
		return 0.00001;		// K2
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_GLambda_Student_PowerSpreadOption_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_GLambda_Student_PowerSpreadOption_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_GLambda_Student_PowerSpreadOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	}
	
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L2]<=0) return ArgumentList_Checking_Result(false," FIRST_UNDERLYING_L2  not positive"); //			positivity of INDEX1
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L3]<=0) return ArgumentList_Checking_Result(false," FIRST_UNDERLYING_L3  not positive"); //			positivity of INDEX1
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L4]<=0) return ArgumentList_Checking_Result(false," FIRST_UNDERLYING_L4  not positive"); //			positivity of INDEX1
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::FIRST_UNDERLYING_L5]<=0) return ArgumentList_Checking_Result(false," FIRST_UNDERLYING_L5  not positive"); //			positivity of INDEX1

//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L2]<=0) return ArgumentList_Checking_Result(false," SECOND_UNDERLYING_L2  not positive"); //			positivity of INDEX1
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L3]<=0) return ArgumentList_Checking_Result(false," SECOND_UNDERLYING_L3  not positive"); //			positivity of INDEX1
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L4]<=0) return ArgumentList_Checking_Result(false," SECOND_UNDERLYING_L4  not positive"); //			positivity of INDEX1
//	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::SECOND_UNDERLYING_L5]<=0) return ArgumentList_Checking_Result(false," SECOND_UNDERLYING_L5  not positive"); //			positivity of INDEX1


	if(abs(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," CORRELATION   not positive"); //	Abs(CORRELATION)<=1 
	if(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::DEGRE]<=0) return ArgumentList_Checking_Result(false," DEGRE    not positive"); //			positivity of DEGRE
 
	if( fabs( a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::NBSTEPS]-
		floor(a[ARM_CF_GLambda_Student_PowerSpreadOption_Formula::NBSTEPS]) ) > K_NEW_DOUBLE_TOL)
			return ArgumentList_Checking_Result(false," NBSteps not an integer"); //					Integer for NbSteps

	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_GLambda_Student_PowerSpreadOption_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///		so the number of parameters is 19
	 if ((rank<0)||(rank>=ARM_CF_GLambda_Student_PowerSpreadOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_SQRT2 
#undef ARM_CF_INVSQRTPI 

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/