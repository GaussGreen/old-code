
/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Functions for   BlackSholes 
 *
 *	\file vanille_bs_formula.cpp
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
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanille_bs_formula.h"
//#include "gpinfra/argconvdefault.h"

#include "expt.h"



CC_BEGIN_NAMESPACE( ARM )


double ARM_CF_BS_Formula::value(const ArgumentList& a)
	{
		double val= BlackSholes_Formula(
			a[ARM_CF_BS_Formula::FORWARD],
			a[ARM_CF_BS_Formula::TOTALVOLATILITY],
			a[ARM_CF_BS_Formula::DISCOUNT],
			a[ARM_CF_BS_Formula::STRIKE],
			a[ARM_CF_BS_Formula::CALLORPUT]
			);
		return val;
	}


double ARM_CF_BS_Formula::value(int i,const ArgumentList& a, double s)
{
	switch (i)
	{
	case ARM_CF_BS_Formula::FORWARD :
		return BlackSholes_Derivative_1(
			a[ARM_CF_BS_Formula::FORWARD],
			a[ARM_CF_BS_Formula::TOTALVOLATILITY],
			a[ARM_CF_BS_Formula::DISCOUNT],
			a[ARM_CF_BS_Formula::STRIKE],
			a[ARM_CF_BS_Formula::CALLORPUT]
			);
	case ARM_CF_BS_Formula::TOTALVOLATILITY :
		return BlackSholes_Derivative_2(
			a[ARM_CF_BS_Formula::FORWARD],
			a[ARM_CF_BS_Formula::TOTALVOLATILITY],
			a[ARM_CF_BS_Formula::DISCOUNT],
			a[ARM_CF_BS_Formula::STRIKE],
			a[ARM_CF_BS_Formula::CALLORPUT]
			);
			
	default :
		return standard_first_derivative(ARM_CF_BS_Formula::value,i,a,s);
	}
}


double ARM_CF_BS_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	switch (i)
	{
	case ARM_CF_BS_Formula::FORWARD :
		{
			switch (j)
			{
			case ARM_CF_BS_Formula::FORWARD :
				return BlackSholes_Derivative_1_1(
					a[ARM_CF_BS_Formula::FORWARD],
			a[ARM_CF_BS_Formula::TOTALVOLATILITY],
			a[ARM_CF_BS_Formula::DISCOUNT],
			a[ARM_CF_BS_Formula::STRIKE],
			a[ARM_CF_BS_Formula::CALLORPUT]
					);
			default :
				return standard_second_derivative(ARM_CF_BS_Formula::value,i,j,a,s1,s2);
			}
		}
	default :
		return standard_second_derivative(ARM_CF_BS_Formula::value,i,j,a,s1,s2);
	}
}


ArgumentList_Checking_Result ARM_CF_BS_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=5) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	 }
	 if (a[ARM_CF_BS_Formula::FORWARD]<=0) return ArgumentList_Checking_Result(false," FORWARD not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_BS_Formula::TOTALVOLATILITY]<0) return ArgumentList_Checking_Result(false," TOTALVOLATILITY not positive"); //			positivity of INDEX1
	 
	 if (a[ARM_CF_BS_Formula::DISCOUNT]<=0) return ArgumentList_Checking_Result(false," DISCOUNT not positive"); //			positivity of INDEX2
	
	 if (a[ARM_CF_BS_Formula::STRIKE ]<0) return ArgumentList_Checking_Result(false," STRIKE not positive"); //		positivity of VOLATILITY1

	 if ((a[ARM_CF_BS_Formula::CALLORPUT]!=K_CALL) && (a[ARM_CF_BS_Formula::CALLORPUT]!=K_PUT))
	 {
		  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); 
	 }

 	return ArgumentList_Checking_Result(true,string(""));
 }


	ArgumentList_Checking_Result ARM_CF_BS_Formula::check_dimension(int rank)
	{
		if ((rank<0)||(rank>=ARM_CF_BS_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




 double ARM_CF_BS_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_BS_Formula::FORWARD :
			return 0.001;		// Forward
		case ARM_CF_BS_Formula::TOTALVOLATILITY :
			return 0.001;		// total volatility
		case ARM_CF_BS_Formula::DISCOUNT :
			return 0.001;		// discount
		case ARM_CF_BS_Formula::STRIKE :
			return 0.001;		//  strike
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_Formula::specific_shift : incorrect input");
	
		}
	}







CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/