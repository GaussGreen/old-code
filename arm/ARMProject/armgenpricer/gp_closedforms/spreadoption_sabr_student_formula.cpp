/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file speradoption_sabr_student_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_sabr_student_formula.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/digital_templated_pricing.h"


#include "expt.h"




CC_BEGIN_NAMESPACE(ARM)


////////////////////////////////////////////////////////////////////////
///
///      Value of the formula (introduces the checking and the ability to automatically derive
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_Student_PowerSpreadOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_SABR_Student_PowerSpreadOption_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_Student_PowerSpreadOption_Formula::value  : bad argsize");
	 }
	 ///		Here the limitation of the dependence in the smile and the copula is that they 
	 ///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	 ///
	 ///		we use here only PowerSpreadOption_Pricing and not Basic_Pricing_With_Limits which is more efficient but which involve the computation of 
	 ///		DistributionZaInverseLimit<smile,copula> that is not simple and prone to error in the case of SABR 
	 ///
	 ArgumentList* CopulaArg=new ArgumentList(
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::CORRELATION],	//CORRELATION
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::DEGRE]		//DEGRE
		 );
	 StudentCopula cop(CopulaArg);
	 double val=structure::PowerSpreadOption_Pricing(
		 ArgumentList(
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX1],	//INDEX1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA1],	//ALPHA1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1],	//BETA1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO1],	//RHO1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU1],		//NU1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG],	//FLAG
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS],	//NBSTEPS
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA_EXP],	//ALPHA_EXP
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA_TANH],	//ALPHA_TANH
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::KB_TANH]		//KB_TANH
		 ),
		 ArgumentList(
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX2],	//INDEX2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA2],	//ALPHA2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2],	//BETA2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO2],	//RHO2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU2],		//NU2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG],	//FLAG
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS],	//NBSTEPS
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA_EXP],	//ALPHA_EXP
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA_TANH],	//ALPHA_TANH
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::KB_TANH]		//KB_TANH
		 ),
		 cop,
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::TIMETOEXPIRATION],	//TIMETOEXPIRATION
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS],
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::A1],	//A1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::B1],	//B1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::K1],	//K1
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::A2],	//A2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::B2],	//B2 
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::K2],	//K2
		 a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALGORITHM]	//Algorithm
		 );
		 
	 return val;
}


////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_Student_PowerSpreadOption_Formula::value(int i,const ArgumentList& a, double s)
///     i vaut:
///		INDEX1,
///		INDEX2,
///		ALPHA1,
///		BETA1,
///		RHO1,
///		NU1,
///		ALPHA2,
///		BETA2,
///		RHO2,
///		NU2,
///		CORRELATION,
///		TIMETOEXPIRATION,
///		A1,
///		B1,
///		K1,
///		A2,
///		B2,
///		K2,

{
	int algorithm=a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALGORITHM];
	int der_i;
	 		
	switch (algorithm)
	{
	case GENERIC_SPREADOPTION :
		{
			return standard_first_derivative(ARM_CF_SABR_Student_PowerSpreadOption_Formula::value,i,a,s);
		}
	case DIGITAL_SPREADOPTION :
		{
			ArgumentList Underlying1(
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX1],	//INDEX1
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA1],	//ALPHA1
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1],	//BETA1
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO1],		//RHO1
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU1],		//NU1
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG],		//FLAG
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS]	//NBSTEPS
				);
			ArgumentList Underlying2 (
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX2],	//INDEX2
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA2],	//ALPHA2
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2],	//BETA2
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO2],		//RHO2
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU2],		//NU2
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG],		//FLAG
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS]	//NBSTEPS
				);
			ArgumentList copula_arg (
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::CORRELATION],	//CORRELATION
				a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::DEGRE]		//DEGRE
				);
			switch(i)
			{
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX1 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA1 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO1 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU1 :
				{
					if(i==ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX1) der_i=ARM_CF_SABR_VanillaOption_Formula::FORWARD;
					else der_i=ARM_CF_SABR_VanillaOption_Formula::ALPHA+(i-ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA1);
					
					return  Templated_Pricing<SABR_smile,StudentCopula>::Digital_Generic_Pricing_First_Derivative_1(
						der_i,
						Underlying1,
						Underlying2,
						copula_arg,
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::TIMETOEXPIRATION],	//TIMETOEXPIRATION
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS],
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::K1],	//K1
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::A2],	//A2
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::B2],	//B2 
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::K2]	//K2
						);
				}
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX2 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA2 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO2 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU2 :
				{
					if(i==ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX2) der_i=ARM_CF_SABR_VanillaOption_Formula::FORWARD;
					else der_i=ARM_CF_SABR_VanillaOption_Formula::ALPHA+(i-ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA2);
					
					return  Templated_Pricing<SABR_smile,StudentCopula>::Digital_Generic_Pricing_First_Derivative_2(
						der_i,
						Underlying1,
						Underlying2,
						copula_arg,
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::TIMETOEXPIRATION],	//TIMETOEXPIRATION
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS],
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::K1],	//K1
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::A2],	//A2
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::B2],	//B2 
						a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::K2]	//K2
						);
				}
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::CORRELATION :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::DEGRE :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::TIMETOEXPIRATION :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::K1 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::A2 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::B2 :
			case ARM_CF_SABR_Student_PowerSpreadOption_Formula::K2 :
				{
					return standard_first_derivative(ARM_CF_SABR_Student_PowerSpreadOption_Formula::value,i,a,s);
				}
			default :
				{	
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_Student_PowerSpreadOption_Formula::value: derivative index  : bad value ");
					break;
				}
			}
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_Student_PowerSpreadOption_Formula::value: ALGORITHM : bad value ");
			break;
		}
		
	}
}


////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_Student_PowerSpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_SABR_Student_PowerSpreadOption_Formula::value,i,j,a,s1,s2);
}



////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_Student_PowerSpreadOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX1 :
		return 0.00001;		// INDEX1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX2 :
		return 0.00001;		// INDEX2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA1 :
		return 0.00001;		// ALPHA1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1 :
		return 0.00001;		// BETA1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO1 :
		return 0.00001;		// RHO1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU1 :
		return 0.00001;		// NU1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA2 :  
		return 0.00001;		// ALPHA2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2 :
		return 0.00001;		// BETA2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO2  :
		return 0.00001;		// RHO2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU2 :
		return 0.00001;		// NU2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::CORRELATION :
		return 0.00001;		// CORRELATION
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::DEGRE :
		return 0.00001;		// DEGRE
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::TIMETOEXPIRATION :
		return 0.00001;		// TIMETOEXPIRATION
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::A1 :
		return 0.00001;		// A1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::B1 :
		return 0.00001;		// B1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::K1 :
		return 0.00001;		// K1
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::A2 :
		return 0.00001;		// A2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::B2  :
		return 0.00001;		// B2
	case ARM_CF_SABR_Student_PowerSpreadOption_Formula::K2 :
		return 0.00001;		// K2
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_Student_PowerSpreadOption_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_SABR_Student_PowerSpreadOption_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_SABR_Student_PowerSpreadOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments : <7 ");
	}
	
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1  not positive"); //			positivity of INDEX1
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2  not positive"); //			positivity of INDEX2
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA1]<=0) return ArgumentList_Checking_Result(false," ALPHA1   not positive"); //			positivity of ALPHA1
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1]<=0) return ArgumentList_Checking_Result(false," BETA1   not positive"); //			positivity of BETA1
	if(abs(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO1])>1.) return ArgumentList_Checking_Result(false," RHO1     not positive"); //	Abs(RHO1)<=1
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU1]<=0) return ArgumentList_Checking_Result(false," NU1      not positive"); //			positivity of NU1
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::ALPHA2]<=0) return ArgumentList_Checking_Result(false," ALPHA2    not positive"); //			positivity of ALPHA2
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2]<=0) return ArgumentList_Checking_Result(false," BETA2   not positive"); //			positivity of BETA2
	if(abs(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::RHO2])>1.) return ArgumentList_Checking_Result(false," RHO2    not positive"); //		Abs(RHO2)<=1 
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NU2]<=0) return ArgumentList_Checking_Result(false," NU2      not positive"); //			positivity of NU2
	if(abs(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," CORRELATION   not positive"); //	Abs(CORRELATION)<=1 
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::DEGRE]<=0) return ArgumentList_Checking_Result(false," DEGRE    not positive"); //			positivity of DEGRE
	if(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::TIMETOEXPIRATION]<=0) return ArgumentList_Checking_Result(false," TIMETOEXPIRATION not positive"); //positivity of TIMETOEXPIRATION
    if(
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::ANALYTIC)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NINTEGRATION)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::AUTOMATIC)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTGEOMETRIC)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTARITHMETIC)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NORMALEXACT)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NORMALGEOMETRIC)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NORMALARITHMETIC)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP2)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_A)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_G)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC1)&&
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC2)
		)	return ArgumentList_Checking_Result(false," Unkown Extended_Flag"); //				LEGAL Extended_Flag
	if( fabs( a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS]-
		floor(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::NBSTEPS]) ) > K_NEW_DOUBLE_TOL)
			return ArgumentList_Checking_Result(false," NBSteps not an integer"); //					Integer for NbSteps

	if ((a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]==ARM_CF_SABR_ImplicitVol_Formula::SABR_A)&&					//		compatibility of the conventions with the Kernel
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1]!=1.0))
	{
		return ArgumentList_Checking_Result(false," ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA1(SABR_A) should be 1"); //
	}
	if ((a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::FLAG]==ARM_CF_SABR_ImplicitVol_Formula::SABR_A)&&					//		compatibility of the conventions with the Kernel
		(a[ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2]!=1.0))
	{
		return ArgumentList_Checking_Result(false," ARM_CF_SABR_Student_PowerSpreadOption_Formula::BETA2(SABR_A) should be 1"); //
	}
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_SABR_Student_PowerSpreadOption_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///		so the number of parameters is 19
	 if ((rank<0)||(rank>=ARM_CF_SABR_Student_PowerSpreadOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_SQRT2 
#undef ARM_CF_INVSQRTPI 

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/