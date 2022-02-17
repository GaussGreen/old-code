/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_s3_spreadoption_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date November 2006
 */
 
#ifndef _GP_CF_BISABR_S3_SPREADOPTION_FORMULA_H
#define _GP_CF_BISABR_S3_SPREADOPTION_FORMULA_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gpclosedforms/powerspreadoption_nonsymmetric.h"


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_BiSABR_S3_SpreadOption_Formula
///			Purpose : Evaluation of spread option under the BiSABR model
///			
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_BiSABR_S3_SpreadOption_Formula
{
	typedef ARM_CF_PowerSpreadOption_GuaussianNonSymmetric_Formula<BiSABR_smile,SABR_smile> structure;
	enum ArgumentType
	{
		INDEX1,				// i=0
		ALPHA1,				// i=1
		BETA1,				// i=2
		RHO1,				// i=3
		NU1,				// i=4
		INDEX2,				// i=5
		ALPHA2,				// i=6
		BETA2,				// i=7
		RHO2,				// i=8
		NU2,				// i=9
		INDEX3,				// i=10
		ALPHA3,				// i=11
		BETA3,				// i=12
		RHO3,				// i=13
		NU3,				// i=14
		RHOS,				// i=15
		RHOV,				// i=16
		RHOC12,				// i=17
		RHOC21,				// i=18
		CORRELATION,		// i=19
		T,					// i=20
		A1,					// i=21
		B1,					// i=22
		K1,					// i=23
		A2,					// i=24
		B2,					// i=25
		K2,					// i=26
		FLAG,				// i=27
		NBSTEPS				// i=28
	};
	enum 
	{ 
		Nb_Parameters =29
	};

	enum
	{ 
		Nb_Derivable_Parameters =27
	};
	static double value(const ArgumentList& a);
		static double value2(const ArgumentList& a);
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
