/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_calib_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou,emezzine
 * Initial version
 *
 */

/*! \file ARM_xl_gp_calib_local.h
 *
 *  \brief file for the generic calibration part in the generic pricer
 *
 *	\author  E Benhamou, Ezzine Mostafa
 *	\version 1.0
 *	\date December 2003
 */

#ifndef ARM_XL_GP_CALIB_LOCAL_H
#define ARM_XL_GP_CALIB_LOCAL_H

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
/// Create a Calibrate Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Calibrate(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_CalibMethodId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Calibrate(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_CalibMethodId);

///////////////////////////////////
/// Create a Calib Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethod_Create(
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_Type,
    LPXLOPER XL_Max_iter,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_Validate,
	LPXLOPER XL_advise);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethod_Create(
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_Type,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_Validate,
	LPXLOPER XL_advise);
//////////////////////////////////////////
/// Create a Calib Method with description
///////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethodWithDesc_Create(
	LPXLOPER XL_MethodType,
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_ModelFitterDesId,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_NbIteration,
	LPXLOPER XL_Validate,
	LPXLOPER XL_advise);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethodWithDesc_Create(
    LPXLOPER XL_MethodType,
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_ModelFitterDesId,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_NbIteration,
	LPXLOPER XL_Validate,
	LPXLOPER XL_advise);

//////////////////////////////////////////
/// Create a Calib Method 2D with description
///////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethod2D_Create(
    LPXLOPER XL_Portfolio1Id,
	LPXLOPER XL_Portfolio2Id,
	LPXLOPER XL_CalibParams1Ids,
	LPXLOPER XL_CalibParams2Ids,
    LPXLOPER XL_ModelFitterDes1Id,
	LPXLOPER XL_ModelFitterDes2Id,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethod2D_Create(
    LPXLOPER XL_Portfolio1Id,
	LPXLOPER XL_Portfolio2Id,
	LPXLOPER XL_CalibParams1Ids,
	LPXLOPER XL_CalibParams2Ids,
    LPXLOPER XL_ModelFitterDes1Id,
	LPXLOPER XL_ModelFitterDes2Id,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb);


///////////////////////////////////
/// Set the detail flag in the calib method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetDetailFlagToCalibMethod(
	LPXLOPER XL_CalibMethodId,
	LPXLOPER XL_DetailFlag );


///////////////////////////////////
//// Get duration from a calib method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetDurationFromCalibMethod(
	LPXLOPER XL_CalibMethodId );


///////////////////////////////////
//// Create a model fitter by type
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Optimizer_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_tol,
    LPXLOPER XL_stepMax,
	LPXLOPER XL_localSearch,
    LPXLOPER XL_printLevel);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Optimizer_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_tol,
    LPXLOPER XL_stepMax,
	LPXLOPER XL_localSearch,
    LPXLOPER XL_printLevel);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Solver_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_xTol,
	LPXLOPER XL_fxTol,
	LPXLOPER XL_gradTol,
    LPXLOPER XL_printLevel);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_Solver_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_xTol,
	LPXLOPER XL_fxTol,
	LPXLOPER XL_gradTol,
    LPXLOPER XL_printLevel);

///////////////////////////////////
/// Create a Numerical Calib Method
///////////////////////////////////
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_NumericalCalibMethod_Create(
    LPXLOPER XL_CalibDateStripId,
    LPXLOPER XL_VanillaSecDensitiesId,
	LPXLOPER XL_PortfolioId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NumericalCalibMethod_Create(
    LPXLOPER XL_CalibDateStripId,
    LPXLOPER XL_VanillaSecDensitiesId,
	LPXLOPER XL_PortfolioId);
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SLNDensityFunctor(
    LPXLOPER XL_Volatility,
    LPXLOPER XL_Shift );
__declspec(dllexport) LPXLOPER WINAPI Local_SLNDensityFunctor(
    LPXLOPER XL_Volatility,
    LPXLOPER XL_Shift );
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABRDensityFunctor(
    LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize);
__declspec(dllexport) LPXLOPER WINAPI Local_SABRDensityFunctor(
    LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize);
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BiSABRDensityFunctor(
    LPXLOPER XL_Alpha1,
	LPXLOPER XL_Beta1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Nu1,
    LPXLOPER XL_Alpha2,
	LPXLOPER XL_Beta2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_RhoS1S2,
	LPXLOPER XL_RhoS1V2,
	LPXLOPER XL_RhoS2V1,
	LPXLOPER XL_RhoV1V2,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize);
__declspec(dllexport) LPXLOPER WINAPI Local_BiSABRDensityFunctor(
    LPXLOPER XL_Alpha1,
	LPXLOPER XL_Beta1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Nu1,
    LPXLOPER XL_Alpha2,
	LPXLOPER XL_Beta2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_RhoS1S2,
	LPXLOPER XL_RhoS1V2,
	LPXLOPER XL_RhoS2V1,
	LPXLOPER XL_RhoV1V2,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize);
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SplineDensityFunctor(
    LPXLOPER XL_Money,
	LPXLOPER XL_Vol,
	LPXLOPER XL_VolType,
	LPXLOPER XL_Smile);
__declspec(dllexport) LPXLOPER WINAPI Local_SplineDensityFunctor(
    LPXLOPER XL_Money,
	LPXLOPER XL_Vol,
	LPXLOPER XL_VolType,
	LPXLOPER XL_Smile);

__declspec(dllexport) LPXLOPER WINAPI Local_NormalHestonDensityFunctor(
	LPXLOPER C_Fwd,
	LPXLOPER C_V0,
	LPXLOPER C_Kappa,
	LPXLOPER C_Theta,
	LPXLOPER C_VVol,
	LPXLOPER C_Rho,
	LPXLOPER C_Level);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NormalHestonDensityFunctor(
	LPXLOPER C_Fwd,
	LPXLOPER C_V0,
	LPXLOPER C_Kappa,
	LPXLOPER C_Theta,
	LPXLOPER C_VVol,
	LPXLOPER C_Rho,
	LPXLOPER C_Level);

__declspec(dllexport) LPXLOPER WINAPI Local_DensityFunctor_CallOption(
	LPXLOPER XL_DensityFunctorId,
    LPXLOPER XL_Forward,
    LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity);
__declspec(dllexport) LPXLOPER WINAPI Local_DensityFunctor_Quantile(
	LPXLOPER XL_DensityFunctorId,
    LPXLOPER XL_Forward,
    LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity);
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IrFwdDensityFunctor();
__declspec(dllexport) LPXLOPER WINAPI Local_IrFwdDensityFunctor();
//--------------------
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VanillaSecurityDensity(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate,
    LPXLOPER XL_EndDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency,
	LPXLOPER XL_DayCount,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight,
	LPXLOPER XL_AdjFwdAdd,
	LPXLOPER XL_AdjFwdMult);
__declspec(dllexport) LPXLOPER WINAPI Local_VanillaSecurityDensity(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate,
    LPXLOPER XL_EndDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency,
	LPXLOPER XL_DayCount,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight,
	LPXLOPER XL_AdjFwdAdd,
	LPXLOPER XL_AdjFwdMult);
//--------------------

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VanillaSecurityDensitySpread(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency1,
	LPXLOPER XL_DayCount1,
	LPXLOPER XL_Frequency2,
	LPXLOPER XL_DayCount2,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight);
__declspec(dllexport) LPXLOPER WINAPI Local_VanillaSecurityDensitySpread(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency1,
	LPXLOPER XL_DayCount1,
	LPXLOPER XL_Frequency2,
	LPXLOPER XL_DayCount2,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight);

///////////////////////////////////
/// Create a Calib Method fo HW2F
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethodHW2F_Create(
	LPXLOPER XL_PortfolioId1,
    LPXLOPER XL_ParamId1,
    LPXLOPER XL_PortfolioId2,
	LPXLOPER XL_ParamId2,
    LPXLOPER XL_PortfolioId3,
	LPXLOPER XL_ParamId3,
	LPXLOPER XL_FlagOptim);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethodHW2F_Create(
	LPXLOPER XL_PortfolioId1,
    LPXLOPER XL_ParamId1,
    LPXLOPER XL_PortfolioId2,
	LPXLOPER XL_ParamId2,
    LPXLOPER XL_PortfolioId3,
	LPXLOPER XL_ParamId3,
	LPXLOPER XL_FlagOptim);

///////////////////////////////////
/// Basket decomposition of a security
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BasketDecomp_Create(
	LPXLOPER XL_SecuritiesId,
    LPXLOPER XL_ModelsId,
    LPXLOPER XL_Weights,
	LPXLOPER XL_NotionalId,
	LPXLOPER XL_ExerDateStripId,
	LPXLOPER XL_ExerFeesId,
	LPXLOPER XL_Side,
	LPXLOPER XL_MkmoId,
	LPXLOPER XL_Method,
	LPXLOPER XL_Strike);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasketDecomp_Create(
	LPXLOPER XL_SecuritiesId,
    LPXLOPER XL_ModelsId,
    LPXLOPER XL_Weights,
	LPXLOPER XL_NotionalId,
	LPXLOPER XL_ExerDateStripId,
	LPXLOPER XL_ExerFeesId,
	LPXLOPER XL_Side,
	LPXLOPER XL_MkmoId,
	LPXLOPER XL_Method,
	LPXLOPER XL_Strike);

__declspec(dllexport) LPXLOPER WINAPI Local_GetData_FromBasket(
	LPXLOPER XL_BasketId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetData_FromBasket(
	LPXLOPER XL_BasketId,
    LPXLOPER XL_KeyId );


///////////////////////////////////
/// SmileViewer of a security
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SmileViewer_Create(
	LPXLOPER XL_SecurityId,
    LPXLOPER XL_ModelId,
	LPXLOPER XL_Moneyness,
	LPXLOPER XL_MoneyType,
	LPXLOPER XL_Strike);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmileViewer_Create(
	LPXLOPER XL_SecurityId,
    LPXLOPER XL_ModelId,
	LPXLOPER XL_Moneyness,
	LPXLOPER XL_MoneyType,
	LPXLOPER XL_Strike);

///////////////////////////////////
/// GenSecurityDensity
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SecurityDensityGen_Create(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_DecStrike,
	LPXLOPER XL_IsDirect,
    LPXLOPER XL_MinProba,
    LPXLOPER XL_MaxProba);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SecurityDensityGen_Create(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_DecStrike,
	LPXLOPER XL_IsDirect,
    LPXLOPER XL_MinProba,
    LPXLOPER XL_MaxProba);


#endif