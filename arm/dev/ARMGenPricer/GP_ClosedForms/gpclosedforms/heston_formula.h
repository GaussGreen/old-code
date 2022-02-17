/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=theta*(lgtvol^2-V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///			V(0)=sig^2
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
//////
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_HESTON_FORMULA_H
#define _GP_CF_HESTON_FORMULA_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result


CC_BEGIN_NAMESPACE(ARM)

/// forward declaration
class ArgumentList;

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Heston_JumpDiffusion_Formula
///			Purpose : Evaluation of avanilla option in the Heston with jumps model
///			Assumptions: lognormal hypothesis with stochastic volatility and normal jumps of the log
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Heston_JumpDiffusion_Formula
{
	enum ArgumentType
	{
		INDEX = 0,					// i=0
		STRIKE,						// i=1
		INITIALVOL,			// i=2
		TIMETOMATURITY,				// i=3
		LONGTERMVOL,				// i=4  omega
		VOLSQUAREREVERTINGSPEED,	// i=5  theta
		VOLSQUAREVOLATILITY,		// i=6  ksi
		CORRELATION,				// i=7  rho
		JUMPPROBABILITY,			// i=8	lambda	
		JUMPMEAN,					// i=9  muJ
		JUMPVOLATILITY,				// i=10  sigmaJ
		CALLORPUT,					// i=11
		NBTERMS						// i=12
	};
	enum Vector_Interpolation_Method
	{
		LINEARINTERPOLATION=0,
		PARABOLICINTERPOLATION
	};
	enum 
	{ 
		Nb_Parameters =13
	};
	enum
	{ 
		Nb_Derivable_Parameters =11
	};
	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Heston_Formula
///			Purpose : Evaluation of avanilla option in the Heston with jumps model
///			Assumptions: lognormal hypothesis with stochastic volatility 
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Heston_Formula
{
	enum ArgumentType
	{
		INDEX = 0,					// i=0
		STRIKE,						// i=1
		INITIALVOLSQUARE,					// i=2
		TIMETOMATURITY,				// i=3
		LONGTERMVOLSQUARE,				// i=4
		VOLSQUAREREVERTINGSPEED,	// i=5
		VOLSQUAREVOLATILITY,		// i=6
		CORRELATION,				// i=7
		CALLORPUT,					// i=8
		NBTERMS						// i=9
	};
	enum Vector_Interpolation_Method
	{
		LINEARINTERPOLATION=0,
		PARABOLICINTERPOLATION
	};
	enum 
	{ 
		Nb_Parameters =10
	};
	enum
	{ 
		Nb_Derivable_Parameters =8
	};
	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
