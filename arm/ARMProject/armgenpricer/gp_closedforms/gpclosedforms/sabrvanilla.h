/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabrvanilla.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_EXTENDED_SABR_FORMULA_H
#define _GP_CF_EXTENDED_SABR_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h"
#include "armdef.h"
#include "gpclosedforms/sabrbdiff1.h"


CC_BEGIN_NAMESPACE(ARM)




struct ARM_CF_SABR_ImplicitVol_Formula
{
	enum ArgumentType
	{
		FORWARD,
		STRIKE,
		MATURITY,
		ALPHA,
		BETA,
		RHO,
		NU,
		EXTENDEDFLAG,
		NBSTEPS,
		ALPHA_EXP,
		ALPHA_TANH,
		KB_TANH
	};
	enum 
	{ 
		Nb_Parameters =12
	};
	enum 
	{ 
		Nb_Derivable_Parameters =7
	};
	
	
	enum Extended_Flag
	{
		ANALYTIC=			K_SABR_ANALYTIC,
		NINTEGRATION=		K_SABR_NINTEGRATION,
		AUTOMATIC=			K_SABR_AUTOMATIC,
		DIRECTEXACT=		K_SABR_DIRECTEXACT,
		DIRECTEXACTSTRIKE=	K_SABR_DIRECTEXACTSTRIKE,
		DIRECTGEOMETRIC=	K_SABR_DIRECTGEOMETRIC,
		DIRECTARITHMETIC=	K_SABR_DIRECTARITHMETIC,
		NORMALEXACT=		K_SABR_NORMALEXACT,
		NORMALGEOMETRIC=	K_SABR_NORMALGEOMETRIC,
		NORMALARITHMETIC=	K_SABR_NORMALARITHMETIC,
		ANALYTICZP0=		K_SABR_ANALYTICZP0,
		ANALYTICZP2=		K_SABR_ANALYTICZP2,
		SABR_IMPLNVOL=		K_SABR_IMPLNVOL,			// pour la compatibilité avec le Kernel, equivalent a DIRECTEXACT  (directe/exact/exact/arithmetic)
		SABR_IMPLNVOL2=		K_SABR_IMPLNVOL2,			// pour la compatibilité avec le Kernel, equivalent a DIRECTEXACTSTRIKE  (directe/exact/exact/strike)
		SABR_A=				K_SABR_ARITH,				// pour la compatibilité avec le Kernel, 
		SABR_G=				K_SABR_GEO,					// pour la compatibilité avec le Kernel, 
		SABR_DIRECT_SC1=	K_SABR_DIRECT_SC1,					// Directe avec strike cutter= Exponentiel
		SABR_DIRECT_SC2=	K_SABR_DIRECT_SC2					// Directe avec strike cutter= Tanh, eqquivalent au strike cutter de Natexis pour les K faible


	};


	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};



struct ARM_CF_SABR_VanillaOption_Formula
{
	enum ArgumentType
	{
		FORWARD,
		STRIKE,
		MATURITY,
		ALPHA,
		BETA,
		RHO,
		NU,
		CALLORPUT,
		EXTENDEDFLAG,
		NBSTEPS,
		ALPHA_EXP,
		ALPHA_TANH,
		KB_TANH
	};
	enum 
	{ 
		Nb_Parameters =13
	};
	enum 
	{ 
		Nb_Derivable_Parameters =7
	};
	enum Extended_Flag
	{
		ANALYTIC=			K_SABR_ANALYTIC,	
		NINTEGRATION=		K_SABR_NINTEGRATION,
		AUTOMATIC=			K_SABR_AUTOMATIC,
		DIRECTEXACT=		K_SABR_DIRECTEXACT,
		DIRECTEXACTSTRIKE=	K_SABR_DIRECTEXACTSTRIKE,
		DIRECTGEOMETRIC=	K_SABR_DIRECTGEOMETRIC,
		DIRECTARITHMETIC=	K_SABR_DIRECTARITHMETIC,
		NORMALEXACT=		K_SABR_NORMALEXACT,
		NORMALGEOMETRIC=	K_SABR_NORMALGEOMETRIC,
		NORMALARITHMETIC=	K_SABR_NORMALARITHMETIC,
		ANALYTICZP0=		K_SABR_ANALYTICZP0,
		ANALYTICZP2=		K_SABR_ANALYTICZP2,
		SABR_IMPLNVOL=		K_SABR_IMPLNVOL,								// pour la compatibilité avec le Kernel, equivalent a DIRECTEXACT  (directe/exact/exact/arithmetic)
		SABR_IMPLNVOL2=		K_SABR_IMPLNVOL2,								// pour la compatibilité avec le Kernel, equivalent a DIRECTEXACTSTRIKE  (directe/exact/exact/strike)
		SABR_A=				K_SABR_ARITH,									// pour la compatibilité avec le Kernel, 
		SABR_G=				K_SABR_GEO,										// pour la compatibilité avec le Kernel,
		SABR_DIRECT_SC1=	K_SABR_DIRECT_SC1,								// Directe avec strike cutter= Exponentiel
		SABR_DIRECT_SC2=	K_SABR_DIRECT_SC2								// Directe avec strike cutter= Tanh, eqquivalent au strike cutter de Natexis pour les K faible

	};


	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};

double SABR_ComputeImpliedVol( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       int SABRflag );

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

