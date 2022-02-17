/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_nummethod_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_nummethod_local.h
 *
 *  \brief file for the numerical method part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#ifndef ARM_XL_GP_NUMMETHOD_LOCAL_H
#define ARM_XL_GP_NUMMETHOD_LOCAL_H

#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"


//////////////////////////////////
/// 
/// In order to factorise code
/// many function are calling
/// common function
/// these functions are not
/// declared here 
/// 
/// only exported functions are 
/// included here to facilitate 
/// the reading of the header
/// 
//////////////////////////////////



///////////////////////////////////
/// Creates a backward induction num method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy );

///////////////////////////////////
/// Creates a backward induction num method
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy );


///////////////////////////////////
/// Creates a forward induction num method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy );

///////////////////////////////////
/// Creates a forward induction num method
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy );


///////////////////////////////////
/// Set Numerical method to model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetNumMethodtoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numMethodId );

	
///////////////////////////////////
/// Set Numerical method to model
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetNumMethodtoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numMethodId );



///////////////////////////////////
/// Set Numerical method to model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetNumerairetoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numeraireId );


///////////////////////////////////
/// Set Numerical method to model
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetNumerairetoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numeraireId );


///////////////////////////////////
/// Create a Numeraire
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Numeraire_Create(
	LPXLOPER XL_NumeraireType,
	LPXLOPER XL_NumeraireTimes);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Numeraire_Create(
	LPXLOPER XL_NumeraireType,
	LPXLOPER XL_NumeraireTimes);


///////////////////////////////////
/// Create a Tree Method or a Tree ND depending of internal flag
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TreeMethod_Create(
    LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbStdDev,
    LPXLOPER XL_MinStdDev,
    LPXLOPER XL_NbMinSteps);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TreeMethod_Create(
    LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbStdDev,
    LPXLOPER XL_MinStdDev,
    LPXLOPER XL_NbMinSteps);


///////////////////////////////////
/// Create a MC Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MCMethod_Create(
	LPXLOPER XL_ItersNb,
	LPXLOPER XL_FixStep,
	LPXLOPER XL_RandGenIds,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
	LPXLOPER XL_SchedulerType,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_ExercBoundCalcId,
	LPXLOPER XL_MaxBucketSize,
	LPXLOPER XL_ImpSamplerType,
	LPXLOPER XL_ImpSamplerDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MCMethod_Create(
	LPXLOPER XL_ItersNb,
	LPXLOPER XL_FixStep,
	LPXLOPER XL_RandGenIds,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
	LPXLOPER XL_SchedulerType,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_ExercBoundCalcId, 
	LPXLOPER XL_MaxBucketSize,
	LPXLOPER XL_ImpSamplerType,
	LPXLOPER XL_ImpSamplerDatas);

///////////////////////////////////
/// Create a PDE Method
///////////////////////////////////

//Initial version
__declspec(dllexport) LPXLOPER WINAPI Local_PDEMethod_Create(
	LPXLOPER XL_MethodName,
	LPXLOPER XL_SchedulerData,
	LPXLOPER XL_NX,
	LPXLOPER XL_GridData,
	LPXLOPER XL_NY,
	LPXLOPER XL_NZ,
	LPXLOPER XL_Theta1,
	LPXLOPER XL_Theta2,
	LPXLOPER XL_Theta3,
	LPXLOPER XL_BoundaryConditionName,
	LPXLOPER XL_Lambda,
	LPXLOPER XL_GridType
	);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PDEMethod_Create(
	LPXLOPER XL_MethodName,
	LPXLOPER XL_SchedulerData,
	LPXLOPER XL_NX,
	LPXLOPER XL_GridData,
	LPXLOPER XL_NY,
	LPXLOPER XL_NZ,
	LPXLOPER XL_Theta1,
	LPXLOPER XL_Theta2,
	LPXLOPER XL_Theta3,
	LPXLOPER XL_BoundaryConditionName,
	LPXLOPER XL_Lambda,
	LPXLOPER XL_GridType);

//ND version
__declspec(dllexport) LPXLOPER WINAPI Local_PdeND_Common(
		LPXLOPER XL_MethodName,
        LPXLOPER XL_SchulerDataId,
		LPXLOPER XL_SchedulerType,
		LPXLOPER XL_SpaceDataId,
		LPXLOPER XL_SchemeDataId,
		LPXLOPER XL_BoundCondName
		);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PdeND_Common(
		LPXLOPER XL_MethodName,
		LPXLOPER XL_SchedulerType,
        LPXLOPER XL_SchulerDataId,
		LPXLOPER XL_SpaceDataId,
		LPXLOPER XL_SchemeDataId,
		LPXLOPER XL_BoundCondName,
		LPXLOPER XL_Lambda
		);

///////////////////////////////////
/// Create a Random nb generator
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_RandGen_Create(
	LPXLOPER XL_genType,
	LPXLOPER XL_algo,
	LPXLOPER XL_baseGen1Id,
	LPXLOPER XL_seed,
	LPXLOPER XL_dim,
	LPXLOPER XL_factorDim,
	LPXLOPER XL_nbOfPoints,
	LPXLOPER XL_nbStdDevs,
	LPXLOPER XL_baseGen2Id,
	LPXLOPER XL_nbFirstDims,
	LPXLOPER XL_order,
	LPXLOPER XL_firstSimulations,
	LPXLOPER XL_firstNbDims);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_RandGen_Create(
	LPXLOPER XL_genType,
	LPXLOPER XL_algo,
	LPXLOPER XL_baseGenId,
	LPXLOPER XL_seed,
	LPXLOPER XL_dim,
	LPXLOPER XL_factorDim,
	LPXLOPER XL_nbOfPoints,
	LPXLOPER XL_nbStdDevs,
	LPXLOPER XL_baseGen2Id,
	LPXLOPER XL_nbFirstDims,
	LPXLOPER XL_order,
	LPXLOPER XL_firstSimulations,
	LPXLOPER XL_firstNbDims);


///////////////////////////////////
/// Create a Random nb generator
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SimpleRandGen_Create(
	LPXLOPER XL_genMode,
	LPXLOPER XL_firstTimes,
	LPXLOPER XL_firstNbDims,
	LPXLOPER XL_skipNbStdDevs
	);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SimpleRandGen_Create(
	LPXLOPER XL_genMode,
	LPXLOPER XL_firstTimes,
	LPXLOPER XL_firstNbDims,
	LPXLOPER XL_skipNbStdDevs
	);

__declspec(dllexport) LPXLOPER Local_RandomGen_DrawVector(
	LPXLOPER XL_RandGenId,
	LPXLOPER XL_size );

/////////////////////////////////////////////////////////////
//// Get the numerical method of a GP model
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetNumMethodFromModel(
	LPXLOPER XL_modelId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetNumMethodFromModel(
	LPXLOPER XL_modelId);

///////////////////////////////////
/// Create an Andersen Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_AMCAndersen_Create(
	LPXLOPER XL_ItersNb, 
	LPXLOPER XL_sortedMaximisation );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AMCAndersen_Create(
	LPXLOPER XL_ItersNb, 
	LPXLOPER XL_sortedMaximisation );


///////////////////////////////////
/// Create an LongstaffSchwartz Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_AMCLongstaffSchwartz_Create(
	LPXLOPER XL_ItersNb );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AMCLongstaffSchwartz_Create(
	LPXLOPER XL_ItersNb );


///////////////////////////////////
/// Create an Tree ND Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TreeND_Create(
	LPXLOPER XL_NbDims,
	LPXLOPER XL_SchedulerType,
    LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
    LPXLOPER XL_TruncatorType,
	LPXLOPER XL_TruncatorDatas,
    LPXLOPER XL_ProbasFlag,
	LPXLOPER XL_ReconnectorType,
	LPXLOPER XL_SmootherType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TreeND_Create(
	LPXLOPER XL_NbDims,
	LPXLOPER XL_SchedulerType,
    LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
    LPXLOPER XL_TruncatorType,
	LPXLOPER XL_TruncatorDatas,
    LPXLOPER XL_ProbasFlag,
	LPXLOPER XL_ReconnectorType,
	LPXLOPER XL_SmootherType);

extern LPXLOPER Local_SetSomethingToModelCommon(
	LPXLOPER XL_modelId,
	LPXLOPER XL_obj2Id,
	const string& Obj2Name,
	long (*Function)(long,long,ARM_result&, long ),
	bool PersistentInXL );

__declspec(dllexport) LPXLOPER WINAPI Local_TreeND_Create(
	LPXLOPER XL_GenSecId,
	LPXLOPER XL_ModelId,
    LPXLOPER XL_InitGuess,
	LPXLOPER XL_LowerBound,
    LPXLOPER XL_UpperBound);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ImpSampler_Optimize(
	LPXLOPER XL_GenSecId,
	LPXLOPER XL_ModelId,
    LPXLOPER XL_InitGuess,
	LPXLOPER XL_LowerBound,
    LPXLOPER XL_UpperBound);


#endif
