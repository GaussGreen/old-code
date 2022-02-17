
/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Functions for   BlackSholes 
 *
 *	\file stochasticvol_ln_formula.cpp
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
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/stochasticvol_ln.h"
#include "gpclosedforms/stochasticvol_ln_formula.h"
#include "gpinfra/gramfunctorstochvolcf.h"
#include "gpinfra/argconvdefault.h"

#include <glob/expt.h>



CC_BEGIN_NAMESPACE( ARM )


double ARM_CF_StochasticVol_LN_Formula::value(const ArgumentList& a)
{
	double val;
	int selector=a[ARM_CF_StochasticVol_LN_Formula::AVERAGINGTYPE];
	switch (selector) 
	{
	case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST :
		{
			val= StochasticVol_LN_Geometric_VanillaOption_with_Reset(
				a[ARM_CF_StochasticVol_LN_Formula::FORWARD],
				a[ARM_CF_StochasticVol_LN_Formula::STRIKE],
				a[ARM_CF_StochasticVol_LN_Formula::MATURITY],
				a[ARM_CF_StochasticVol_LN_Formula::DRIFT],
				a[ARM_CF_StochasticVol_LN_Formula::VOLATILITY],
				a[ARM_CF_StochasticVol_LN_Formula::VOLDRIFT],
				a[ARM_CF_StochasticVol_LN_Formula::VOLVOLATILITY],
				a[ARM_CF_StochasticVol_LN_Formula::AVERAGINGPERIOD],
				a[ARM_CF_StochasticVol_LN_Formula::RESET],
				a[ARM_CF_StochasticVol_LN_Formula::CALLORPUT],
				a[ARM_CF_StochasticVol_LN_Formula::NBTERMS]
				);
			break;
		}
	case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST :
		{
			val= StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(
				a[ARM_CF_StochasticVol_LN_Formula::FORWARD],
				a[ARM_CF_StochasticVol_LN_Formula::STRIKE],
				a[ARM_CF_StochasticVol_LN_Formula::MATURITY],
				a[ARM_CF_StochasticVol_LN_Formula::DRIFT],
				a[ARM_CF_StochasticVol_LN_Formula::VOLATILITY],
				a[ARM_CF_StochasticVol_LN_Formula::VOLDRIFT],
				a[ARM_CF_StochasticVol_LN_Formula::VOLVOLATILITY],
				a[ARM_CF_StochasticVol_LN_Formula::AVERAGINGPERIOD],
				a[ARM_CF_StochasticVol_LN_Formula::RESET],
				a[ARM_CF_StochasticVol_LN_Formula::CALLORPUT],
				a[ARM_CF_StochasticVol_LN_Formula::NBTERMS]
				);
			break;
		}
		
		
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : flag  bad input :");
			break;
		}
	}
	return val;	
}


double ARM_CF_StochasticVol_LN_Formula::value(int i,const ArgumentList& a, double s)
{
	return standard_first_derivative(ARM_CF_StochasticVol_LN_Formula::value,i,a,s);
}


double ARM_CF_StochasticVol_LN_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	return standard_second_derivative(ARM_CF_StochasticVol_LN_Formula::value,i,j,a,s1,s2);
}


ArgumentList_Checking_Result ARM_CF_StochasticVol_LN_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_StochasticVol_LN_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 if (a[ARM_CF_StochasticVol_LN_Formula::FORWARD				]<=0) return ArgumentList_Checking_Result(false," FORWARD				not positive"); //		positivity of FORWARD		
	 if (a[ARM_CF_StochasticVol_LN_Formula::STRIKE				]<=0) return ArgumentList_Checking_Result(false," STRIKE		 		not positive"); //		positivity of STRIKE		
	 if (a[ARM_CF_StochasticVol_LN_Formula::MATURITY			]<0) return ArgumentList_Checking_Result(false," MATURITY		 		not positive"); //		positivity of MATURITY			
     if (a[ARM_CF_StochasticVol_LN_Formula::VOLATILITY			]<=0) return ArgumentList_Checking_Result(false," VOLATILITY	 		not positive"); //		positivity of VOLATILITY		
	 if (a[ARM_CF_StochasticVol_LN_Formula::VOLVOLATILITY		]<0) return ArgumentList_Checking_Result(false," VOLVOLATILITY 		not positive"); //		positivity of VOLVOLATILITY
	 if (a[ARM_CF_StochasticVol_LN_Formula::AVERAGINGPERIOD		]<0) return ArgumentList_Checking_Result(false," AVERAGINGPERIOD 		not positive"); //		positivity of AVERAGINGPERIOD
	 if (a[ARM_CF_StochasticVol_LN_Formula::RESET		]<0)		return ArgumentList_Checking_Result(false," RESET 		not positive"); //		positivity of RESET

	 if ((a[ARM_CF_StochasticVol_LN_Formula::CALLORPUT]!=K_CALL) && (a[ARM_CF_StochasticVol_LN_Formula::CALLORPUT]!=K_PUT))	  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); 

 	return ArgumentList_Checking_Result(true,string(""));
 }


	ArgumentList_Checking_Result ARM_CF_StochasticVol_LN_Formula::check_dimension(int rank)
	{
		if ((rank<0)||(rank>=ARM_CF_StochasticVol_LN_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




 double ARM_CF_StochasticVol_LN_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_StochasticVol_LN_Formula::FORWARD			:return 0.001;		// FORWARD		
		case ARM_CF_StochasticVol_LN_Formula::STRIKE			:return 0.001;		// STRIKE		
		case ARM_CF_StochasticVol_LN_Formula::MATURITY			:return 0.001;		// MATURITY		
		case ARM_CF_StochasticVol_LN_Formula::DRIFT				:return 0.001;		// DRIFT		
		case ARM_CF_StochasticVol_LN_Formula::VOLATILITY		:return 0.001;		// VOLATILITY	
		case ARM_CF_StochasticVol_LN_Formula::VOLDRIFT			:return 0.001;		// VOLDRIFT		
		case ARM_CF_StochasticVol_LN_Formula::VOLVOLATILITY		:return 0.001;		// VOLVOLATILITY
		case ARM_CF_StochasticVol_LN_Formula::AVERAGINGPERIOD	:return 0.001;		// AVERAGINGPERIOD
		case ARM_CF_StochasticVol_LN_Formula::RESET				:return 0.001;		// RESET
		
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_StochasticVol_LN_Formula::specific_shift : incorrect input");
	
		}
	}







CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/