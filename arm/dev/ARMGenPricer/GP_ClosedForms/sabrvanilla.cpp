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


#include "gpclosedforms/sabrvanilla.h"

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/sabr_newversion.h"

#include <glob/expt.h>


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////////////
///
///				Interface of the formula
///
///////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
///
///      Value of the implicit vol formula 
///  (introduces the checking and the ability to automatically derive
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_ImplicitVol_Formula::value(const ArgumentList& a)
{
	double val;
	int argsize=a.size();
	if (argsize!=ARM_CF_SABR_ImplicitVol_Formula::Nb_Parameters)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_ImplicitVol_Formula::value  : bad argsize");
	}
	
	int flag=a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG];
	switch(flag)
	{
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC1 :
		{
		val=SABR_StrikeCuterExp_ImplicitVol(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU],	//NU
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA_EXP]	//ALPHA_EXP
				);
			break;
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC2 :
		{
		val=SABR_StrikeCuterTanh_ImplicitVol(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU],	//NU
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA_TANH],	//ALPHA_TANH
				a[ARM_CF_SABR_ImplicitVol_Formula::KB_TANH]	//KB_TANH
				);
			break;
		}
	case ARM_CF_SABR_ImplicitVol_Formula::NINTEGRATION :
	case ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0 :
		{
			val=SABR_NumericalIntegration_ImplicitVol(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU],	//NU
				a[ARM_CF_SABR_ImplicitVol_Formula::NBSTEPS]	//NBSTEPS
				);
			break;
		}
	case ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP2 :
		{
			val=SABR_NumericalIntegration_ImplicitVol2(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU],	//NU
				a[ARM_CF_SABR_ImplicitVol_Formula::NBSTEPS]	//NBSTEPS
				);
			break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTGEOMETRIC :
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTARITHMETIC :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
		{
			val=SABR_implicit_vol_direct(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU],	//NU
				flag
				);
			break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::NORMALEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_G :
	case ARM_CF_SABR_ImplicitVol_Formula::NORMALGEOMETRIC :
	case ARM_CF_SABR_ImplicitVol_Formula::NORMALARITHMETIC :
		{
			val=SABR_implicit_vol_normal(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU],	//NU
				flag
				);
			break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_A :
		{
			val=SABR_BetEqOne_CompatibleKernel_ImplVol(
				a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_ImplicitVol_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_ImplicitVol_Formula::NU]	//NU
				);
			break;
			
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_ImplicitVol_Formula::value: methodflag : bad input :");
			break;
		}
	}
	return val;
}




////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_ImplicitVol_Formula::value(int i,const ArgumentList& a, double s)
{
	int flag=a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG];
	if (flag==ARM_CF_SABR_ImplicitVol_Formula::ANALYTIC)
		
	{
		switch(i)
		{
		case ARM_CF_SABR_ImplicitVol_Formula::FORWARD :
			{
				return SABR_ImpVol_DerForward(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE], 
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA], 
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO],
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}
		case ARM_CF_SABR_ImplicitVol_Formula::STRIKE :
			{
				return SABR_ImpVol_DerStrike(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE], 
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA],
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO],
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}

		case ARM_CF_SABR_ImplicitVol_Formula::MATURITY :
			{
				return SABR_ImpVol_DerMaturity(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD], 
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY], 
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA], 
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA],
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO], 
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}
		case ARM_CF_SABR_ImplicitVol_Formula::ALPHA :
			{
				return SABR_ImpVol_DerAlpha(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA], 
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA],
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO],
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}
		case ARM_CF_SABR_ImplicitVol_Formula::BETA :
			{
				return SABR_ImpVol_DerBeta(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA], 
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO], 
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}
		case ARM_CF_SABR_ImplicitVol_Formula::RHO :
			{
				return SABR_ImpVol_DerRho(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE], 
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY],
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA],
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA], 
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO],
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}
		case ARM_CF_SABR_ImplicitVol_Formula::NU :
			{
				return SABR_ImpVol_DerNu(
					a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD],
					a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE],
					a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY], 
					a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA], 
					a[ARM_CF_SABR_ImplicitVol_Formula::BETA],
					a[ARM_CF_SABR_ImplicitVol_Formula::RHO], 
					a[ARM_CF_SABR_ImplicitVol_Formula::NU]);
			}
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_ImplicitVol_Formula::value : incorrect dimension");
			
		}


	}
	else
	{
	 		return standard_first_derivative(ARM_CF_SABR_ImplicitVol_Formula::value,i,a,s);
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_ImplicitVol_Formula::value(int i,
        int j,
        const ArgumentList& a,
        double s1,
        double s2)
{
	return standard_second_derivative(ARM_CF_SABR_ImplicitVol_Formula::value,i,j,a,s1,s2);
}



////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_ImplicitVol_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_SABR_ImplicitVol_Formula::FORWARD :
		return 0.00001;		// FORWARD
	case ARM_CF_SABR_ImplicitVol_Formula::STRIKE :
		return 0.000001;		// STRIKE
	case ARM_CF_SABR_ImplicitVol_Formula::MATURITY :
		return 0.00001;		// MATURITY
	case ARM_CF_SABR_ImplicitVol_Formula::ALPHA :
		return 0.00001;		// ALPHA
	case ARM_CF_SABR_ImplicitVol_Formula::BETA :
		return 0.00001;		// BETA
	case ARM_CF_SABR_ImplicitVol_Formula::RHO :
		return 0.00001;		// RHO
	case ARM_CF_SABR_ImplicitVol_Formula::NU :  
		return 0.00001;		// NU
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_ImplicitVol_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_SABR_ImplicitVol_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()<ARM_CF_SABR_ImplicitVol_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments : <7 ");
	}
	
	if(a[ARM_CF_SABR_ImplicitVol_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::FORWARD  not positive"); //			positivity of FORWARD
	if(a[ARM_CF_SABR_ImplicitVol_Formula::STRIKE]<0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::STRIKE  not positive"); //			positivity of STRIKE
	if(a[ARM_CF_SABR_ImplicitVol_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::MATURITY   not positive"); //		positivity of MATURITY
	if(a[ARM_CF_SABR_ImplicitVol_Formula::ALPHA]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::ALPHA   not positive"); //			positivity of ALPHA
	if(a[ARM_CF_SABR_ImplicitVol_Formula::BETA]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::BETA     not positive"); //			positivity of BETA
	if(abs(a[ARM_CF_SABR_ImplicitVol_Formula::RHO])>1.) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::RHO   not positive"); //		Abs(RHO)<=1 
	if(a[ARM_CF_SABR_ImplicitVol_Formula::NU]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::NU not positive"); //				positivity of NU
	if(
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::ANALYTIC)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NINTEGRATION)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::AUTOMATIC)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTGEOMETRIC)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::DIRECTARITHMETIC)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NORMALEXACT)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NORMALGEOMETRIC)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::NORMALARITHMETIC)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP2)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_A)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_G)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC1)&&
		(a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC2)
		)	return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::Unkown Extended_Flag"); //				LEGAL Extended_Flag
	if( fabs( a[ARM_CF_SABR_ImplicitVol_Formula::NBSTEPS]-floor(a[ARM_CF_SABR_ImplicitVol_Formula::NBSTEPS]) ) > K_NEW_DOUBLE_TOL)
			return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::NBSteps not an integer"); //			Integer for NbSteps
	if ((a[ARM_CF_SABR_ImplicitVol_Formula::EXTENDEDFLAG]==ARM_CF_SABR_ImplicitVol_Formula::SABR_A)&&					//		compatibility of the conventions with the Kernel
		(1.0-a[ARM_CF_SABR_ImplicitVol_Formula::BETA])>1e-8)
	{
		return ArgumentList_Checking_Result(false," ARM_CF_SABR_ImplicitVol_Formula::BETA(SABR_A) should be 1"); //
	}
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_SABR_ImplicitVol_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///		are described by alpha,beta,rho,nu + flag, steps  
		///		so the number of parameters is 9
	 if ((rank<0)||(rank>=ARM_CF_SABR_ImplicitVol_Formula::Nb_Derivable_Parameters)) 
		 return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

////////////////////////////////////////////////////////////////////////
///
///      Value of the VanillaOption formula 
///  (introduces the checking and the ability to automatically derive
///
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_VanillaOption_Formula::value(const ArgumentList& a)
	{
	double val;
	int argsize=a.size();
	if (argsize!=ARM_CF_SABR_VanillaOption_Formula::Nb_Parameters)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_VanillaOption_Formula::value  : bad argsize");
	}
	
	int flag=a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG];
	switch(flag)
	{
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC1 :
		{
		val=SABR_StrikeCuterExp_VanillaOption(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU],	//NU
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA_EXP],	//ALPHA_EXP
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC2 :
		{
		val=SABR_StrikeCuterTanh_VanillaOption(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU],	//NU
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA_TANH],	//ALPHA_TANH
				a[ARM_CF_SABR_VanillaOption_Formula::KB_TANH],	//KB_TANH
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
		}

	case ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0 :
	case ARM_CF_SABR_ImplicitVol_Formula::NINTEGRATION :
		{
			val=VanillaOptionFromCall(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				
				SABR_NumericalIntegration_Call(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU],	//NU
				a[ARM_CF_SABR_VanillaOption_Formula::NBSTEPS]	//NBSTEPS
				),
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
		}
	case ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP2 :
		{
			val=VanillaOptionFromCall(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				
				SABR_NumericalIntegration_Call2(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU],	//NU
				a[ARM_CF_SABR_VanillaOption_Formula::NBSTEPS]	//NBSTEPS
				),
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTGEOMETRIC :
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTARITHMETIC :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
		{
			val=VanillaOptionFromCall(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				
				BS(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				
				SABR_implicit_vol_direct(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU],	//NU
				flag
				)),
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
			
		}
		
	case ARM_CF_SABR_ImplicitVol_Formula::NORMALEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_G :
	case ARM_CF_SABR_ImplicitVol_Formula::NORMALGEOMETRIC :
	case ARM_CF_SABR_ImplicitVol_Formula::NORMALARITHMETIC :
		{
			val=VanillaOptionFromCall(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				
				BS(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				
				SABR_implicit_vol_normal(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::BETA],	//BETA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU],	//NU
				flag
				)),
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_A :
		{
			val=VanillaOptionFromCall(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				
				BS(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				
				SABR_BetEqOne_CompatibleKernel_ImplVol(
				a[ARM_CF_SABR_VanillaOption_Formula::FORWARD],	//FORWARD,
				a[ARM_CF_SABR_VanillaOption_Formula::STRIKE],	//STRIKE,
				a[ARM_CF_SABR_VanillaOption_Formula::MATURITY],	//MATURITY,
				a[ARM_CF_SABR_VanillaOption_Formula::ALPHA],	//ALPHA,
				a[ARM_CF_SABR_VanillaOption_Formula::RHO],	//RHO,
				a[ARM_CF_SABR_VanillaOption_Formula::NU]	//NU
				)),
				a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]	//CALLORPUT
				);
			break;
			
		}
		
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_VanillaOption_Formula::value: methodflag : bad input :");
			break;
		}
	 }
	 return val;
}




////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_VanillaOption_Formula::value(int i,const ArgumentList& a, double s)
{	
		return standard_first_derivative(ARM_CF_SABR_VanillaOption_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_VanillaOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_SABR_VanillaOption_Formula::value,i,j,a,s1,s2);
}



////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_SABR_VanillaOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_SABR_VanillaOption_Formula::FORWARD :
		return 0.0001;		// FORWARD
	case ARM_CF_SABR_VanillaOption_Formula::STRIKE :
		return 0.00001;		// STRIKE
	case ARM_CF_SABR_VanillaOption_Formula::MATURITY :
		return 0.001;		// MATURITY
	case ARM_CF_SABR_VanillaOption_Formula::ALPHA :
		return 0.001;		// ALPHA
	case ARM_CF_SABR_VanillaOption_Formula::BETA :
		return 0.0001;		// BETA
	case ARM_CF_SABR_VanillaOption_Formula::RHO :
		return 0.001;		// RHO
	case ARM_CF_SABR_VanillaOption_Formula::NU :  
		return 0.001;		// NU
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_SABR_ImplicitVol_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_SABR_VanillaOption_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_SABR_VanillaOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	}
	if(a[ARM_CF_SABR_VanillaOption_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::FORWARD  not positive"); //			positivity of FORWARD
	if(a[ARM_CF_SABR_VanillaOption_Formula::STRIKE]<0) 
	{
		double vv=a[ARM_CF_SABR_VanillaOption_Formula::STRIKE];
		return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::STRIKE  not positive");															//	positivity of STRIKE
	}
	if(a[ARM_CF_SABR_VanillaOption_Formula::MATURITY]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::MATURITY   not positive");		//	positivity of MATURITY
	if(a[ARM_CF_SABR_VanillaOption_Formula::ALPHA]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::ALPHA   not positive");			//	positivity of ALPHA
	if(a[ARM_CF_SABR_VanillaOption_Formula::BETA]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::BETA     not positive");			//	positivity of BETA
	if(abs(a[ARM_CF_SABR_VanillaOption_Formula::RHO])>1.) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::RHO   not positive");			//	Abs(RHO)<=1 
	if(a[ARM_CF_SABR_VanillaOption_Formula::NU]<=0) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::NU not positive");					//	positivity of NU
	 if ((a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]!=K_CALL) &&(a[ARM_CF_SABR_VanillaOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::CALLORPUT should be K_CALL or K_PUT");											//	Should be Call or Put

	if(
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::ANALYTIC)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::NINTEGRATION)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::AUTOMATIC)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::DIRECTEXACT)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::DIRECTEXACTSTRIKE)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::DIRECTGEOMETRIC)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::DIRECTARITHMETIC)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::NORMALEXACT)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::NORMALGEOMETRIC)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::NORMALARITHMETIC)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::ANALYTICZP0)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::ANALYTICZP2)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::SABR_IMPLNVOL)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::SABR_IMPLNVOL2)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::SABR_A)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::SABR_G)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::SABR_DIRECT_SC1)&&
		(a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]!=ARM_CF_SABR_VanillaOption_Formula::SABR_DIRECT_SC2)
		)	return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::Unkown Extended_Flag"); //				LEGAL Extended_Flag
	if( fabs( a[ARM_CF_SABR_VanillaOption_Formula::NBSTEPS]-floor(a[ARM_CF_SABR_VanillaOption_Formula::NBSTEPS]) ) > K_NEW_DOUBLE_TOL)
			return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::NBSteps not an integer"); //					Integer for NbSteps

	if ((a[ARM_CF_SABR_VanillaOption_Formula::EXTENDEDFLAG]==ARM_CF_SABR_VanillaOption_Formula::SABR_A)&&					//		compatibility of the conventions with the Kernel
		(1.0-a[ARM_CF_SABR_VanillaOption_Formula::BETA])>1e-8)
	{
		return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::BETA(SABR_A) should be 1"); //
	}
	
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_SABR_VanillaOption_Formula::check_dimension(int rank)
{
	 	///		Here the limitation of the dependence in the smile and the copula is that they 
		///		are described by alpha,beta,rho,nu + flag, steps  
		///		so the number of parameters is 9
	 if ((rank<0)||(rank>=ARM_CF_SABR_VanillaOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," ARM_CF_SABR_VanillaOption_Formula::check_dimension Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/