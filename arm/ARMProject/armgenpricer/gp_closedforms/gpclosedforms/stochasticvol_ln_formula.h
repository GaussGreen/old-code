/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Header for the   BlackSholes templates and related
 *
 *	\file stochasticvol_ln_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _GP_CF_STOCHASTICVOL_LN_FORMULA_H
#define _GP_CF_STOCHASTICVOL_LN_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <gpbase/argconv.h>
#include "basic_distributions.h" /// for ArgumentList_Checking_Result

CC_BEGIN_NAMESPACE( ARM )

struct ARM_CF_StochasticVol_LN_Formula
{
	enum ArgumentType
	{
		FORWARD,
		STRIKE,
		MATURITY,
		DRIFT,
		VOLATILITY,
		VOLDRIFT,
		VOLVOLATILITY,
		AVERAGINGPERIOD,
		RESET,
		AVERAGINGTYPE,
		CALLORPUT,
		NBTERMS
	};
	enum
	{ 
		Nb_Parameters =12
	};
	enum
	{ 
		Nb_Derivable_Parameters =9
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);
	static int nb_arguments();
	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);

	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);

};

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
