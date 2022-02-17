
/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  argconv table : classes and templates 
 *
 *	\file argconvdefault_CF.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */



#ifndef _GP_CF_ARGCONVDEFAULT_H
#define _GP_CF_ARGCONVDEFAULT_H



#include "firsttoinc.h"
#include "gpbase/port.h"

#include <gpbase/argconv.h>

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_BS_Formula_ArgType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_BS_Formula_ArgType;

extern const ARM_ArgConv ARM_ArgConv_SpreadDigitalOption_Formula_DerivedStruct;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SpreadDigitalOption_Formula_DerivedStruct;

extern const ARM_ArgConv ARM_ArgConv_SABR_ImplicitVol_Formula_Extended_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SABR_ImplicitVol_Formula_Extended_Flag;

extern const ARM_ArgConv ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag;

extern const ARM_ArgConv ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_InOut_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_EuropeanBarriere_Formula_InOut_Flag;

extern const ARM_ArgConv ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag;

extern const ARM_ArgConv ARM_ArgConvGP_CFBS_PartialTime_Barriere_End_Formula_OptionType_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFBS_PartialTime_End_Barriere_Formula_OptionType_Flag;

extern const ARM_ArgConv ARM_ArgConvGP_CFHeston_Vector_InterpolationMethod_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFHeston_Vector_InterpolationMethod_Flag;

extern const ARM_ArgConv ARM_ArgConvGP_CFNonParametricExtrapolationMethod_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFNonParametricExtrapolationMethod_Flag;


extern const ARM_ArgConv ARM_ArgConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag;
extern const ARM_ArgConvReverse ARM_ArgConvReverseGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
