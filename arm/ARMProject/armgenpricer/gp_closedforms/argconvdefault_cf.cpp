/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  argconv table : methods and functions 
 *
 *	\file argconvdefault_CF.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/argconvdefault_CF.h"

#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/spreadoption_normal_formula.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_sabr.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/barriere_bs_formula.h"
#include "gpclosedforms/heston_formula.h"
#include "gpclosedforms/optimization1.h"
#include "gpclosedforms/spreadoption_nonparametric_formula.h"




using namespace std;


CC_BEGIN_NAMESPACE( ARM )

ARGConvTable BS_Formula_ArgType_Table[] = 
{
	/// name			/// number
	"FORWARD",			ARM_CF_BS_Formula::FORWARD,
	"TOTALVOLATILITY",	ARM_CF_BS_Formula::TOTALVOLATILITY,
	"DISCOUNT",			ARM_CF_BS_Formula::DISCOUNT,
	"STRIKE",			ARM_CF_BS_Formula::STRIKE,
	"CALLORPUT",		ARM_CF_BS_Formula::CALLORPUT,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};


const ARM_ArgConv ARM_ArgConv_BS_Formula_ArgType( BS_Formula_ArgType_Table, "BS_Formula_ArgType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_BS_Formula_ArgType( BS_Formula_ArgType_Table, "BS_Formula_ArgType" );


ARGConvTable SpreadDigitalOption_Formula_DerivedStruct_Table[] = 
{
	/// name			/// number
	"DIGITALOPTION",			ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION,
	"SPREADOPTION",				ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,
	"PAYFIRST",					ARM_CF_SpreadDigitalOption_Formula::PAYFIRST,
	"PAYSECOND",				ARM_CF_SpreadDigitalOption_Formula::PAYSECOND,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///													SABR Vanilla option
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


const ARM_ArgConv ARM_ArgConv_SpreadDigitalOption_Formula_DerivedStruct( SpreadDigitalOption_Formula_DerivedStruct_Table, 
																		"SpreadDigitalOption_Formula_DerivedStruct" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SpreadDigitalOption_Formula_DerivedStruct( SpreadDigitalOption_Formula_DerivedStruct_Table, 
																					  "SpreadDigitalOption_Formula_DerivedStruct" );

ARGConvTable SABR_ImplicitVol_Formula_Extended_Flag_Table[] = 
{
	/// name			/// number
	"ANALYTIC",			ARM_CF_SABR_ImplicitVol_Formula::ANALYTIC,
	"NINTEGRATION",		ARM_CF_SABR_ImplicitVol_Formula::NINTEGRATION,
	"AUTOMATIC",		ARM_CF_SABR_ImplicitVol_Formula::AUTOMATIC,
	"DIRECTEXACT",		ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT,
	"DIRECTEXACTSTRIKE",ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE,
	"DIRECTGEOMETRIC",	ARM_CF_SABR_ImplicitVol_Formula::DIRECTGEOMETRIC,
	"DIRECTARITHMETIC",	ARM_CF_SABR_ImplicitVol_Formula::DIRECTARITHMETIC,
	"NORMALEXACT",		ARM_CF_SABR_ImplicitVol_Formula::NORMALEXACT,
	"NORMALGEOMETRIC",	ARM_CF_SABR_ImplicitVol_Formula::NORMALGEOMETRIC,
	"NORMALARITHMETIC",	ARM_CF_SABR_ImplicitVol_Formula::NORMALARITHMETIC,
	"ANALYTICZP0",		ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0,
	"ANALYTICZP2",		ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP2,
	"SABR_IMPLNVOL",	ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL,
	"SABR_IMPLNVOL2",	ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2,
	"SABR_A",			ARM_CF_SABR_ImplicitVol_Formula::SABR_A,
	"SABR_G",			ARM_CF_SABR_ImplicitVol_Formula::SABR_G,
	"SABR_DIRECT_SC1",	ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC1,
	"SABR_DIRECT_SC2",	ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC2,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};


const ARM_ArgConv ARM_ArgConv_SABR_ImplicitVol_Formula_Extended_Flag( SABR_ImplicitVol_Formula_Extended_Flag_Table, 
																		"SABR_ImplicitVol_Formula_Extended_Flag" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SABR_ImplicitVol_Formula_Extended_Flag( SABR_ImplicitVol_Formula_Extended_Flag_Table, 
																					  "SABR_ImplicitVol_Formula_Extended_Flag" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///													BS Single barrier option
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ARGConvTable BS_EuropeanBarriere_Formula_OptionType_Table[] = 
{
	/// name			/// number
	"DOWN_AND_IN",			ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN,
	"UP_AND_IN",			ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN,
	"DOWN_AND_OUT",			ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT,
	"UP_AND_OUT",			ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag( BS_EuropeanBarriere_Formula_OptionType_Table, 
																		"BS_EuropeanBarriere_Formula_OptionType_Flag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag(BS_EuropeanBarriere_Formula_OptionType_Table, 
																					  "BS_EuropeanBarriere_Formula_OptionType_Flag" );




ARGConvTable BS_EuropeanBarriere_Formula_InOutFlag_Table[] = 
{
	/// name			/// number
	"IN",		ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN,
	"OUT",		ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};


const ARM_ArgConv ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_InOut_Flag( BS_EuropeanBarriere_Formula_InOutFlag_Table, 
																		"BS_EuropeanBarriere_Formula_InOutFlag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_EuropeanBarriere_Formula_InOut_Flag( BS_EuropeanBarriere_Formula_InOutFlag_Table, 
																					  "BS_EuropeanBarriere_Formula_InOutFlag" );
ARGConvTable BS_EuropeanBarriere_Formula_UpDownFlag_Table[] = 
{
	/// name			/// number
	"UP",		ARM_CF_BS_SingleBarriere_Formula::UP,
	"DOWN",		ARM_CF_BS_SingleBarriere_Formula::DOWN,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};


const ARM_ArgConv ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag( BS_EuropeanBarriere_Formula_UpDownFlag_Table, 
																		"BS_EuropeanBarriere_Formula_UpDownFlag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag( BS_EuropeanBarriere_Formula_UpDownFlag_Table, 
																					  "BS_EuropeanBarriere_Formula_UpDownFlag" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///													Partial time BS Single barrier option
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ARGConvTable BS_PartialTime_Barriere_End_Formula_OptionType_Table[] = 
{
	/// name			/// number
	"CROSS_AND_IN",				ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_IN,
	"CROSS_AND_OUT",			ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_OUT,
	"INSIDE_UP_AND_IN",			ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_IN,
	"INSIDE_UP_AND_OUT",		ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_OUT,
	"INSIDE_DOWN_AND_IN",		ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_IN,
	"INSIDE_DOWN_AND_OUT",		ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_OUT,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConvGP_CFBS_PartialTime_Barriere_End_Formula_OptionType_Flag( BS_PartialTime_Barriere_End_Formula_OptionType_Table, 
																		"BS_PartialTime_Barriere_End_Formula_OptionType_Flag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_PartialTime_End_Barriere_Formula_OptionType_Flag(BS_PartialTime_Barriere_End_Formula_OptionType_Table, 
																					  "BS_PartialTime_Barriere_End_Formula_OptionType_Flag" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///													Heston Vector Formula: Interpolation Method
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ARGConvTable Heston_Vector_InterpolationMethod_Table[] = 
{
	/// name			/// number
	"LINEAR",				ARM_CF_Heston_JumpDiffusion_Formula::LINEARINTERPOLATION,
	"PARABOLIC",			ARM_CF_Heston_JumpDiffusion_Formula::PARABOLICINTERPOLATION,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConvGP_CFHeston_Vector_InterpolationMethod_Flag( Heston_Vector_InterpolationMethod_Table, 
																		"Heston_Vector_InterpolationMethod_Flag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFHeston_Vector_InterpolationMethod_Flag(Heston_Vector_InterpolationMethod_Table, 
																					  "Heston_Vector_InterpolationMethod_Flag" );



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///													NonParametric Methods: Extrapolations
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ARGConvTable NonParametricExtrapolationMethod_Table[] = 
{
	/// name			/// number
	"VOLTAILINDEX",				ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX,
	"VOLCONSTANT",			ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConvGP_CFNonParametricExtrapolationMethod_Flag( Heston_Vector_InterpolationMethod_Table, 
																		"Heston_Vector_InterpolationMethod_Flag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFNonParametricExtrapolationMethod_Flag(Heston_Vector_InterpolationMethod_Table, 
																					  "Heston_Vector_InterpolationMethod_Flag" );


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///													Optimization : Available Algorithm
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ARGConvTable Optimization_ObjectiveFuntion_Algorithm_Table[] = 
{
	/// name			/// number
	"NLIN_LSQ",				Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ,
	"LSQ_DERIV",			Optimization_ObjectiveFuntion::NAG_OPT_LSQ_DERIV,
	"NLIN_LSQ_CONSTRAINT",	Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ_CONSTRAINT,
	"LSQ_CHECK_DERIV",		Optimization_ObjectiveFuntion::NAG_OPT_LSQ_CHECK_DERIV,
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag( Optimization_ObjectiveFuntion_Algorithm_Table, 
																		"Optimization_ObjectiveFuntion_Algorithm_Flag" );
const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag(Optimization_ObjectiveFuntion_Algorithm_Table, 
																					  "Optimization_ObjectiveFuntion_Algorithm_Flag" );



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

