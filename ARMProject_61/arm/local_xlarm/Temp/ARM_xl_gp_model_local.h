/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_xl_gp_model_local.h
 *
 *  \brief file for the model part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#ifndef ARM_XL_GP_MODEL_LOCAL_H
#define ARM_XL_GP_MODEL_LOCAL_H

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
/// Creates an IR Fwd Model
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_IRFwdMod_Create(
	LPXLOPER XL_ZeroCurveId );
///////////////////////////////////
/// Creates an IR Fwd Model
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRFwdMod_Create(
	LPXLOPER XL_ZeroCurveId );


///////////////////////////////////
/// Creates an INF Fwd Model
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfFwdMod_Create(
	LPXLOPER XL_InfCurveId);
///////////////////////////////////
/// Creates an INF Fwd Model
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfFwdMod_Create(
	LPXLOPER XL_InfCurveId);

///////////////////////////////////
/// Creates an Inflation Equity Model
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InflationEquityModel_Create(
	LPXLOPER XL_InfCurveId,LPXLOPER XL_Param1Id,LPXLOPER XL_Param2Id);

///////////////////////////////////
/// Creates an Inflation Equity Model
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InflationEquityModel_Create(
	LPXLOPER XL_InfCurveId,LPXLOPER XL_Param1Id,LPXLOPER XL_Param2Id);

///----------------------------------------------
///----------------------------------------------
///             SABR  Equity Model
/// Inputs :
///     ZcCurve
///     spot( double)
///		model params Id
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SABREquityModel_Create(
	LPXLOPER XL_ZcCurveId,
	LPXLOPER XL_Spot,
	LPXLOPER XL_modelParamsId );

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABREquityModel_Create(
	LPXLOPER XL_ZcCurveId,
	LPXLOPER XL_Spot,
	LPXLOPER XL_modelParamsId );




///////////////////////////////////
/// Create a Model Parameter
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamTimes,
	LPXLOPER XL_ParamValues,
	LPXLOPER XL_ModelParamName,
    LPXLOPER XL_LowerBoundary,
	LPXLOPER XL_UpperBoundary,
    LPXLOPER XL_InterpolMethod,
	LPXLOPER XL_AdviseBreakPointTimes);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamTimes,
	LPXLOPER XL_ParamValues,
    LPXLOPER XL_ModelParamName,
    LPXLOPER XL_LowerBoundary,
	LPXLOPER XL_UpperBoundary,
    LPXLOPER XL_InterpolMethod,
	LPXLOPER XL_AdviseBreakPointTimes);


///////////////////////////////////
/// Create a Cst Model Parameter
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CstModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamValue,
	LPXLOPER XL_AdviseBreakPointTimes );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CstModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamValue,
	LPXLOPER XL_AdviseBreakPointTimes  );

///////////////////////////////////
/// Create a Trigo Correl Parameter
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TrigoCorrelParam_Create(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DateStripId,
	LPXLOPER XL_Theta,
	LPXLOPER XL_interpolatorName,
    LPXLOPER XL_lowerBound,
    LPXLOPER XL_upperBound );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TrigoCorrelParam_Create(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DateStripId,
	LPXLOPER XL_Theta,
	LPXLOPER XL_interpolatorName,
    LPXLOPER XL_lowerBound,
    LPXLOPER XL_upperBound );


///////////////////////////////////
/// Get a model param from a pricing model
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_PricingModel_GetModelParam(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ModelParamDataType,
	LPXLOPER XL_ModelParamIndex,
	LPXLOPER XL_FactorNb );


///////////////////////////////////
/// Create a HW1F model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HW1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER C_MeanReversionParamId,
	LPXLOPER XL_Flags);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HW1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER C_MeanReversionParamId,
	LPXLOPER XL_Flags);



///////////////////////////////////
/// Create a HW2F model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HW2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SigmaRatioParamId,
	LPXLOPER XL_MeanReversionSpreadParamId,
	LPXLOPER XL_CorrelationParamId,
	LPXLOPER XL_Flags);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HW2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SigmaRatioParamId,
	LPXLOPER XL_MeanReversionSpreadParamId,
	LPXLOPER XL_CorrelationParamId,
	LPXLOPER XL_Flags);


///////////////////////////////////
/// Create a MF model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MarkovFunctionalModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MarkovFunctionalModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId);


///////////////////////////////////
/// Create a SFRM model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_volType,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_ShiftConvPortId,
	LPXLOPER XL_NonParamDrift);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_volType,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_ShiftConvPortId,
	LPXLOPER XL_NonParamDrift);

__declspec(dllexport) LPXLOPER Local_SFRMModelVolSwapVolFRADump(
	LPXLOPER XL_SwaptionId,
	LPXLOPER XL_SFRMId);



///////////////////////////////////
/// Set Fix Scheduler SFRM model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetSFRMFixScheduler(
	LPXLOPER XL_SFRMModId,
	LPXLOPER XL_StartDate);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetSFRMFixScheduler(
	LPXLOPER XL_SFRMModId,
	LPXLOPER XL_StartDate);

///////////////////////////////////
/// Create a Hybrid Basis Fwd IR model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HybridBasisFwdIRModel_Create(
	LPXLOPER XL_refIRModelId,
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_basisZcCurveId,
	LPXLOPER XL_forexId,
    LPXLOPER XL_modelNames);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HybridBasisFwdIRModel_Create(
	LPXLOPER XL_refIRModelId,
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_basisZcCurveId,
	LPXLOPER XL_forexId,
    LPXLOPER XL_modelNames);

///////////////////////////////////
/// Set Curve in IR model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Model_ZCCurveSet(
	LPXLOPER XL_IRModelId,
	LPXLOPER XL_zeroCurveId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Model_ZCCurveSet(
	LPXLOPER XL_IRModelId,
	LPXLOPER XL_zeroCurveId);

///////////////////////////////////
/// Create a QGM1F model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QGM1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER C_MeanReversionParamId,
	LPXLOPER C_SkewParamId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QGM1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER C_MeanReversionParamId,
	LPXLOPER C_SkewParamId);


///////////////////////////////////
/// Create a QGM2F model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QGM2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsVecFactor1Id,
	LPXLOPER XL_ParamsVecFactor2Id);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QGM2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsVecFactor1Id,
	LPXLOPER XL_ParamsVecFactor2Id);

///////////////////////////////////
/// Create a Q1F model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QModel1F_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_DegenerateInHW);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QModel1F_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_DegenerateInHW);


///////////////////////////////////
/// Create a Q1F model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QModel1FAna_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_QParamId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QModel1FAna_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_QParamId);



//////////////////////////////////////////////////
//// Function to create a surface model param
//////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SurfaceModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_ModelParamName,
	LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
	LPXLOPER XL_AdviseBreakPointTimes );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SurfaceModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_ModelParamName,
	LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
	LPXLOPER XL_AdviseBreakPointTimes );



//////////////////////////////////////////////////
//// Function to create a surface List model param
//////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SurfaceListModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_Index,
	LPXLOPER XL_SurfaceListId,
	LPXLOPER XL_ModelParamName );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SurfaceListModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_Index,
	LPXLOPER XL_SurfaceListId,
	LPXLOPER XL_ModelParamName );



//////////////////////////////////////////////////
//// Function to create an Heston Model
//////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Heston_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Heston_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

//////////////////////////////////////////////////
//// Function to create a BS Model
//////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BS_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BS_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );


//////////////////////////////////////////////////
//// Function to create a CEV Model
//////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CEV_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

//////////////////////////////////////////////////
//// Function to create a Merton Model
//////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Merton_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Merton_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

//////////////////////////////////////////////////
//// Function to create a Normal Model
//////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Normal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

//////////////////////////////////////////////////
//// Function to create a SABR Model
//////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ImpliedVolType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABR_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ImpliedVolType);

//////////////////////////////////////////////////
//// Function to create a SLN Model
//////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SLN_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SLN_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_GetVarianceSqueeze(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_detail);


///----------------------------------------------
///----------------------------------------------
///             HW1F Model Param Create
/// Inputs :
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HW1FModelParam_Create( LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HW1FModelParam_Create( LPXLOPER XL_ModelParamsId );


///----------------------------------------------
///----------------------------------------------
///             QNF Model Param Create
/// Inputs :
///     QParam Id
///     Vector of model params Id
///     CorrelMat Id
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_QNFFModelParam_Create(
	LPXLOPER XL_QParamId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_CorrelMatId,
	LPXLOPER XL_DegenerateInHW);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QNFFModelParam_Create(
	LPXLOPER XL_QParamId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_CorrelMatId,
	LPXLOPER XL_DegenerateInHW);

///----------------------------------------------
///----------------------------------------------
///             QNF Model Create
/// Inputs :
///     Zc curve Id
///     QNF Model param Id
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_QNFFModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QNFFModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId );



///----------------------------------------------
///----------------------------------------------
///             Q1F FX Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		For curve Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_Q1FModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Q1FModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

///----------------------------------------------
///----------------------------------------------
///             Heston Eq Model Create
/// Inputs :
///     Zc curve Id
///     Heston Model param Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HestonModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HestonModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_Spot );

///----------------------------------------------
///----------------------------------------------
///             Q1F Eq Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_Q1FModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Q1FModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_Spot );



///----------------------------------------------
///----------------------------------------------
///             BS  FX Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Fwd curve Id
///		For curve Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_BSModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

///----------------------------------------------
///----------------------------------------------
///             Heston  FX Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Fwd curve Id
///		For curve Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HestonModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HestonModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

///----------------------------------------------
///----------------------------------------------
///             Heston  FX Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Fwd curve Id
///		For curve Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_SABRModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABRModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_Spot );


///----------------------------------------------
///----------------------------------------------
///             BS Eq Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Fwd curve Id
///		Spot
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_BSModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_Spot );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_Spot );


///----------------------------------------------
///----------------------------------------------
///             Model Name Map Create
/// Inputs :
///     vector of names
///     vector of models id
///		vector of vector of string
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER Local_ModelNameMap_Create(
	LPXLOPER XL_Names,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_OtherModelNames );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ModelNameMap_Create(
	LPXLOPER XL_Names,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_OtherModelNames );



///----------------------------------------------
///----------------------------------------------
///             Multi-Assets Model Create
/// Inputs :
///     Model Name Map
///     Correlation matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER Local_MultiAssetsModel_Create(
	LPXLOPER XL_ModelNameMapId,
	LPXLOPER XL_CorrelationMatrixId,
	LPXLOPER XL_MultiAssetsModelName);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MultiAssetsModel_Create(
	LPXLOPER XL_ModelNameMapId,
	LPXLOPER XL_CorrelationMatrixId,
	LPXLOPER XL_MultiAssetsModelName);

///----------------------------------------------
///----------------------------------------------
///             Get Model Map from Hybrid Model
/// Inputs :
///     Model
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER Local_PricingModel_GetModelMap(
	LPXLOPER XL_ModelId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PricingModel_GetModelMap(
	LPXLOPER XL_ModelId);

///----------------------------------------------
///----------------------------------------------
///             Set Model Map from Hybrid Model
/// Inputs :
///     Model Name Map
///     Model
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PricingModel_SetModelMap(
	LPXLOPER XL_modelId,
	LPXLOPER XL_ModelMapId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PricingModel_SetModelMap(
	LPXLOPER XL_modelId,
	LPXLOPER XL_ModelMapId );

///----------------------------------------------
///----------------------------------------------
///             Create a 2 IR + FX Hybrid Model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_Create2IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Create2IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId );


///----------------------------------------------
///----------------------------------------------
///             Create a 1 IR + FX Hybrid Model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_Create1IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Create1IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId );

///----------------------------------------------
///----------------------------------------------
///             Create a HW HW Quanto Model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HWHWQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWHWQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag);

///----------------------------------------------
///----------------------------------------------
///             Create a HW HW2F Quanto Model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HWHW2FQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWHW2FQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag);


///----------------------------------------------
///----------------------------------------------
///             Create a Fwd Margin Model
/// Inputs :
///		basis curve
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_CreateFwdMarginModel(
	LPXLOPER XL_basisZcCurveId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateFwdMarginModel(
	LPXLOPER XL_basisZcCurveId );


///----------------------------------------------
///----------------------------------------------
///             Local_SetRefModelNameToMultiAsset
///		function to set the reference model to a multi-asset
///			using its name
/// Inputs :
///			Multi-asset Model
///			Name
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_SetRefModelNameToMultiAsset(
	LPXLOPER XL_Name,
	LPXLOPER XL_MultiAssetModel );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetRefModelNameToMultiAsset(
	LPXLOPER XL_Name,
	LPXLOPER XL_MultiAssetModel );

///----------------------------------------------
///----------------------------------------------
///             Local Normal Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_LocalNormal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalNormal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId);

///----------------------------------------------
///----------------------------------------------
///             Local Shifted LogNormal Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_LocalSLN_Model_Create(
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalSLN_Model_Create(
	LPXLOPER XL_ModelParamsId);

///----------------------------------------------
///----------------------------------------------
///             Local Model Calibration
/// Inputs :
///     MultiAssets model Id
///     Name of Local Model embedded in MultiAssets
///		Portfolio Id
///		Eval Dates
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_Local_Model_Calibrate(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_EvalDates);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Local_Model_Calibrate(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_EvalDates);

__declspec(dllexport) LPXLOPER WINAPI Local_Local_Model_VarSqueeze(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName);

__declspec(dllexport) LPXLOPER WINAPI Local_LocalModel_CalibrateFunctional(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_SecIds,
	LPXLOPER XL_DensIds,
	LPXLOPER XL_Rescaling);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalModel_CalibrateFunctional(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_SecIds,
	LPXLOPER XL_DensIds,
	LPXLOPER XL_Rescaling);


///----------------------------------------------
///----------------------------------------------
///             Market IR Model
/// Inputs :
///     Market Datas (Mkt data Manager + Keys)
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_MarketIRModel_Create(
	LPXLOPER XL_MarketDatas, LPXLOPER XL_VnsPricingMethod);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MarketIRModel_Create(
	LPXLOPER XL_MarketDatas, LPXLOPER XL_VnsPricingMethod);


///----------------------------------------------
///----------------------------------------------
///             Smiled FRM Model
/// Inputs :
///     ZcCurve, Correl, Hump, NbFactors
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_SmiledFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_SwaptionApprox);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_SwaptionApprox);

///----------------------------------------------
///----------------------------------------------
///             HWSV 1F Model
/// Inputs :
///     ZcCurve, ModelParams, SolverType, SolverParams
///		FormulaType, FormulaParams
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HWSV1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverType,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSV1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverType,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO);


///----------------------------------------------
///----------------------------------------------
///             HWSV 2F Model
/// Inputs :
///     ZcCurve, ModelParams, SolverParams
///		FormulaType, FormulaParams
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_HWSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO);


///----------------------------------------------
///----------------------------------------------
///             EQHWSV Param Model
/// Inputs :
///		ModelParams
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_EQHWSV_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EQHWSV_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId);


///----------------------------------------------
///----------------------------------------------
///             EQHWSV Numeric Method
/// Inputs :
///     IntStep		Gauss legendre Nb Steps
///		ImAxis		Imaginary axis integration
///		MaxDecay	Critera of discretization of time depending functions
///
///----------------------------------------------
///---------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_EQHWSV_NumMethodsCreate(
	LPXLOPER XL_IntStep,
	LPXLOPER XL_ImAxis,
	LPXLOPER XL_MaxDecay);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EQHWSV_NumMethodsCreate(
	LPXLOPER XL_IntStep,
	LPXLOPER XL_ImAxis,
	LPXLOPER XL_MaxDecay);

///----------------------------------------------
///----------------------------------------------
///             EQHWSV
/// Inputs :
///     EQHWSV_ModelParams
///		EQHWSV_NumMethods
///		EQHWSV_MktDatas
///
///----------------------------------------------
///---------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_EQHWSV_Create(
	LPXLOPER XL_Currency,
	LPXLOPER XL_ModelParams,
	LPXLOPER XL_NumMethods,
	LPXLOPER XL_MktDatas,
	LPXLOPER XL_Dilatation);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EQHWSV_Create(
	LPXLOPER XL_Currency,
	LPXLOPER XL_ModelParams,
	LPXLOPER XL_NumMethods,
	LPXLOPER XL_MktDatas,
	LPXLOPER XL_Dilatation);

///----------------------------------------------
///----------------------------------------------
///             Smiled Market Model
/// Inputs :
///     Pattern, ZcCurve, Correl, Hump, NbFactors
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_SmiledMarketModel_Create(
	LPXLOPER XL_Pattern,
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_CalibProxy,
	LPXLOPER XL_Recorrel);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledMarketModel_Create(
	LPXLOPER XL_Pattern,
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_CalibProxy,
	LPXLOPER XL_Recorrel);

__declspec(dllexport) LPXLOPER WINAPI Local_SmiledMarketModelDS_Create(
	LPXLOPER XL_CalibPattern,
	LPXLOPER XL_StartDate,
	LPXLOPER XL_EndDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_IndexFreq,
	LPXLOPER XL_IndexType,	
	LPXLOPER XL_ResetTiming,		
	LPXLOPER XL_DayCount,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_FwdRule,
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledMarketModelDS_Create(
	LPXLOPER XL_CalibPattern,
	LPXLOPER XL_StartDate,
	LPXLOPER XL_EndDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_IndexFreq,
	LPXLOPER XL_IndexType,	
	LPXLOPER XL_ResetTiming,		
	LPXLOPER XL_DayCount,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_FwdRule,
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap);

///----------------------------------------------
///----------------------------------------------
///             SVBGM Model
/// Inputs :
///     ZcCurve, shift, alpha, nu, rho, correl(taux,taux), correl(taux,vol), correl(vol,vol), NbFactors
///
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_SVBGMModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_alphaId,
	LPXLOPER XL_nuId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_rvcorrelParamId,
	LPXLOPER XL_vvcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_Proxy);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SVBGMModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_alphaId,
	LPXLOPER XL_nuId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_rvcorrelParamId,
	LPXLOPER XL_vvcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_Proxy);


__declspec(dllexport) LPXLOPER WINAPI Local_BGMSV1FModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRevId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_localRhoCalib,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BGMSV1FModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRevId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_LocalCalibration,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_localRhoCalib,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy);

__declspec(dllexport) LPXLOPER WINAPI Local_BGMSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_BetaCorrelId,
	LPXLOPER XL_V01,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_V02,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_LocalRho1Calib,
	LPXLOPER XL_LocalRho2Calib,
	LPXLOPER XL_Shift,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BGMSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_BetaCorrelId,
	LPXLOPER XL_V01,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_V02,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_LocalRho1Calib,
	LPXLOPER XL_LocalRho2Calib,
	LPXLOPER XL_Shift,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy);

__declspec(dllexport) LPXLOPER WINAPI Local_SVMMSpreadModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRevId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SVMMSpreadModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRevId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio);

__declspec(dllexport) LPXLOPER WINAPI Local_BiSVMM_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BiSVMM_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId);

__declspec(dllexport) LPXLOPER WINAPI Local_HWxSVMMSpread_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_hw2fId,
	LPXLOPER XL_corrIndexEndTimes,
	LPXLOPER XL_constantCrossCorrel);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWxSVMMSpread_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_hw2fId,
	LPXLOPER XL_corrIndexEndTimes,
	LPXLOPER XL_constantCrossCorrel);

__declspec(dllexport) LPXLOPER WINAPI Local_HWSBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_correlationMatrix);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_correlationMatrix);

__declspec(dllexport) LPXLOPER WINAPI Local_HWSVBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_correlationMatrix);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSVBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_correlationMatrix);

#endif


