/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Header for the   BlackSholes templates and related
 *
 *	\file vanille_bs_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _GP_CF_VANILLE_BS_FORMULA_H
#define _GP_CF_VANILLE_BS_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <gpbase/argconv.h>
#include "basic_distributions.h" /// for ArgumentList_Checking_Result

CC_BEGIN_NAMESPACE( ARM )

struct ARM_CF_BS_Formula
{
	enum ArgumentType
	{
		FORWARD,
		TOTALVOLATILITY,
        DISCOUNT,
        STRIKE,
		CALLORPUT
	};
	enum
	{ 
		Nb_Parameters =5
	};
	enum
	{ 
		Nb_Derivable_Parameters =4
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);
	static int nb_arguments();
	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);

	static double callimplicit_totalvolatility(double forward, double discount, double strike,
		double callprice,double opt, double accuracy);
	static double callimplicit_totalvolatility1(double forward, double discount, double strike,
		double callprice,double opt, double accuracy);
	static double callimplicit_totalvolatility2(double forward, double discount, double strike,
		double callprice,double opt, double accuracy);
	
	static double callimplicit_totalvolatility_DerOpt(double forward, double discount,
		double strike, double callprice,double opt,double accuracy);
	static double callimplicit_totalvolatility_DerForward(double forward, double discount,
		double strike, double CallPut,double opt,double accuracy);
	static double ARM_CF_BS_Formula::callimplicit_totalvolatility_DerStrike(double forward, double discount,
		double strike, double CallPut,double opt,double accuracy);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);
	static double BSphi(double m, double u);
	static double BSphiInverse(double param,double objective, double accuracy);
	static double BSphiInverse1(double param,double objective, double accuracy);
	static double BSphiInverse2(double param,double objective, double accuracy);

	static double BSphiDerivative(double m, double u);
	static double BSphiDerivative2(double m, double u);


private:

	static void funcd(const double x,const double param,const double objective,double &fn, double &df);
	static void funcd2(const double x,const double param,const double objective,double &fn, double &df);
	static double findroot(void funcd(const double, const double, const double,double &, double &), const double x1, const double x2,
			const double xacc,const double param,const double objective);




};

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
