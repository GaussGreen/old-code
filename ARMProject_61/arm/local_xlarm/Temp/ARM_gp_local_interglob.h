/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_gp_local_interglob.h,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */


////////////////////////////////////////////////
/// This file is for global interpolation and
/// various conversion routine for the generic
/// pricer... Everything specific to the generic
/// pricer should be included here and not in ARM_local_interglob.cpp
/// to avoid exporting the generic pricer in the 
/// local dllarm project (activex project)
/// since the interglob.h file is included in the 
/// local dllarm project (activex project)
////////////////////////////////////////////////

#ifndef _ARM_GP_LOCAL_INTERGLOB_H
#define _ARM_GP_LOCAL_INTERGLOB_H

#include <ARM\libarm\ARM_result.h>
#include <libCCtools++\CCString_STL.h>

extern long ARM_ConvGPNumeraire( const CCString& numeraireType, ARM_result& result);
extern long ARM_ConvGPModelParam( const CCString& modelParamType, ARM_result& result);
extern long ARM_ConvGPModelParamDataType( const CCString& modelParamType, ARM_result& result);
extern long ARM_ConvGPCalibType( const CCString& calibType, ARM_result& result);
extern long ARM_ConvGPCalib2DDirection( const CCString& calibDirection, ARM_result& result);
extern long ARM_ConvGPTargetFuncType( const CCString& calibType, ARM_result& result);
extern long ARM_ConvGP_CFSpreadDigitalOption_Formula_DerivedStruct( const CCString& nType, ARM_result& result);
extern long ARM_ConvGP_CFSABR_ImplicitVol_Formula_Extended_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ConvGP_CFBS_EuropeanBarriere_Formula_InOut_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ConvGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ArgConvGP_CFBS_PartialTime_Barriere_End_Formula_OptionType_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ArgConvGP_CFHeston_Vector_InterpolationMethod_Flag( const CCString& nType, ARM_result& result);
extern long ARM_ConvGPModelParamDataType( const CCString& modelParamDataType, ARM_result& result);
extern long ARM_ConvGPInfSwoptComputationMethod( const CCString& type, ARM_result& result);
extern long ARM_ConvGPSamplerType( const CCString& samplerType, ARM_result& result);
extern long ARM_ConvGPTruncatorType( const CCString& truncatorType, ARM_result& result);
extern long ARM_ConvGPSchedulerType( const CCString& schedulerType, ARM_result& result);
extern long ARM_ConvGPReconnectorType( const CCString& reconnectorType, ARM_result& result);
extern long ARM_ConvGPSmootherType( const CCString& smootherType, ARM_result& result);
extern long ARM_ConvGPImpSamplerType( const CCString& impSamplerType, ARM_result& result);
extern long ARM_ConvGPPathSchemeType( const CCString& pathSchemeType, ARM_result& result);
extern long ARM_ConvGPODESolverType( const CCString& solverType, ARM_result& result);
extern long ARM_ConvHWSVFormulaType( const CCString& formulaType, ARM_result& result);

#endif	

