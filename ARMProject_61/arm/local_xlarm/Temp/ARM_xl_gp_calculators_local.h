/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_calculators_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARM_XL_GP_CALCULATORS_LOCAL_H
#define ARM_XL_GP_CALCULATORS_LOCAL_H


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

__declspec(dllexport) LPXLOPER WINAPI Local_DateStripCombiner_Create(
	LPXLOPER XL_dateStripIds,
	LPXLOPER XL_mergeFuncName );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DateStripCombiner_Create(
	LPXLOPER XL_dateStripIds,
	LPXLOPER XL_mergeFuncName );

__declspec(dllexport) LPXLOPER WINAPI Local_DataStripCombiner_GetData(
    LPXLOPER XL_dateStripCombinerId,
    LPXLOPER XL_dataType,
	LPXLOPER XL_dateStripNb );

__declspec(dllexport) LPXLOPER WINAPI Local_DataStripCombiner_GetMergeData(
    LPXLOPER XL_dateStripCombinerId );

__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_GetPricingData(
	LPXLOPER XL_GenCalculatorId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_GetData(
	LPXLOPER XL_GenCalculatorId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenCalculator_GetData(
	LPXLOPER XL_GenCalculatorId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_SetData(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenCalculator_SetData(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_fixEndDate,
    LPXLOPER XL_fixDayCount,
    LPXLOPER XL_cpnResetGap,
    LPXLOPER XL_leverage,
    LPXLOPER XL_cpnMin,
    LPXLOPER XL_cpnMax,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
    LPXLOPER XL_exerGap,
    LPXLOPER XL_nbNonCall,
    LPXLOPER XL_exerFee,
    LPXLOPER XL_flags,
    LPXLOPER XL_FundNominal);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_fixEndDate,
    LPXLOPER XL_fixDayCount,
    LPXLOPER XL_cpnResetGap,
    LPXLOPER XL_leverage,
    LPXLOPER XL_cpnMin,
    LPXLOPER XL_cpnMax,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
    LPXLOPER XL_exerGap,
    LPXLOPER XL_nbNonCall,
    LPXLOPER XL_exerFee,
    LPXLOPER XL_flags,
    LPXLOPER XL_FundNominal);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCcyFromGenCalculator(
          LPXLOPER XL_calcId,
          LPXLOPER XL_CcyType);

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_Initialize(
	LPXLOPER XL_crfId,
	LPXLOPER XL_mktmanagerId,
    LPXLOPER XL_toCalSigma,
    LPXLOPER XL_toCalMrs,
	LPXLOPER XL_strikeTypeToCalMRS,
    LPXLOPER XL_toAdjKcap,
    LPXLOPER XL_toAdjKfloor,
    LPXLOPER XL_pricingModel,
    LPXLOPER XL_toCalSkew,
	LPXLOPER XL_Kshift,
	LPXLOPER XL_KFrontier
    );

__declspec(dllexport) LPXLOPER WINAPI PXL_Local_CRFCalculator_Initialize(
	LPXLOPER XL_crfId,
	LPXLOPER XL_mktmanagerId,
    LPXLOPER XL_toCalSigma,
    LPXLOPER XL_toCalMrs,
	LPXLOPER XL_strikeTypeToCalMRS,
    LPXLOPER XL_toAdjKcap,
    LPXLOPER XL_toAdjKfloor,
    LPXLOPER XL_pricingModel,
    LPXLOPER XL_toCalSkew,
    LPXLOPER XL_Kshift,
	LPXLOPER XL_KFrontier
	);


__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_GetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_GetMRS(LPXLOPER XL_crfId);

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_GetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_GetMRS(LPXLOPER XL_crfId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_SetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetAutoCalFlags(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags );

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetOneCalFlag(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flag ,
	LPXLOPER XL_type);

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetAutoCalFlagsAndClone(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_SetAutoCalFlagsAndClone(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags );

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetProductToPrice(
	LPXLOPER XL_crfId,
	LPXLOPER XL_ProdToPrice );

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_Update(
	LPXLOPER XL_gcId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType);

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_leverage,
	LPXLOPER XL_lifeTimeCapTarget,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundNominal,
	LPXLOPER XL_fees,
    LPXLOPER XL_asOf);

__declspec(dllexport) LPXLOPER WINAPI Local_TarnSetOutputFlags(
    LPXLOPER XL_SecurityId,
    LPXLOPER XL_ProductFlags);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_leverage,
	LPXLOPER XL_lifeTimeCapTarget,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundNominal,
	LPXLOPER XL_fees,
    LPXLOPER XL_asOf);

__declspec(dllexport) LPXLOPER WINAPI Local_TARNSnowBallCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_coupon0,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_leverage,
	LPXLOPER XL_levPrev,
	LPXLOPER XL_lifeTimeCapTarget,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundingCcyName);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNSnowBallCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_coupon0,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_leverage,
	LPXLOPER XL_levPrev,
	LPXLOPER XL_lifeTimeCapTarget,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundingCcyName);

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_GetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_SetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNCalculator_GetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNCalculator_SetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_Update(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_GetPricingData(
	LPXLOPER XL_TARNCalculatorId,
    LPXLOPER XL_KeyId );

////////////////////////////////////////////////////////////////////////////////
///////////////////caption calculator
////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_Create(
   LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_RequiredCpnDatas,
	LPXLOPER XL_coupon,
	LPXLOPER XL_FundIdxTerm,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_exerStyle,
	LPXLOPER XL_Notional,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_cpnSpread,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_fundSpread,
	LPXLOPER XL_factorNb,
	LPXLOPER XL_CalibMod,
	LPXLOPER XL_flags);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CaptionCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_RequiredCpnDatas,
	LPXLOPER XL_coupon,
	LPXLOPER XL_FundIdxTerm,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_exerStyle,
	LPXLOPER XL_Notional,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_cpnSpread,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_fundSpread,
	LPXLOPER XL_factorNb,
	LPXLOPER XL_CalibMod,
	LPXLOPER XL_flags);

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_GetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_SetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CaptionCalculator_GetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CaptionCalculator_SetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_Update(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_GetPricingData(
	LPXLOPER XL_CaptionCalculatorId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_PRDCCalculator_Create(
    LPXLOPER XL_prdcId,
	LPXLOPER XL_modelId,
    LPXLOPER XL_otherMktDataIds,
    LPXLOPER XL_schedulerDatas,
    LPXLOPER XL_truncatorDatas,
    LPXLOPER XL_columnsToPrice,
    LPXLOPER XL_markovianDriftSamplerFlag,
    LPXLOPER XL_fxLocalModelFlag,
    LPXLOPER XL_calibType,
    LPXLOPER XL_calibDatas,
    LPXLOPER XL_basisIRCalibFlag,
	LPXLOPER XL_marginFlag);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRDCCalculator_Create(
    LPXLOPER XL_prdcId,
	LPXLOPER XL_modelId,
    LPXLOPER XL_otherMktDataIds,
    LPXLOPER XL_schedulerDatas,
    LPXLOPER XL_truncatorDatas,
    LPXLOPER XL_columnsToPrice,
    LPXLOPER XL_markovianDriftSamplerFlag,
    LPXLOPER XL_fxLocalModelFlag,
    LPXLOPER XL_calibType,
    LPXLOPER XL_calibDatas,
    LPXLOPER XL_basisIRCalibFlag,
	LPXLOPER XL_marginFlag);

__declspec(dllexport) LPXLOPER WINAPI Local_PRCSCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_initialFx,
	LPXLOPER XL_Cpns,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRCSCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_initialFx,
	LPXLOPER XL_Cpns,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_PRDKOCalculator_Create(
	LPXLOPER XL_Dates,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_CpnCurves,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRDKOCalculator_Create(
	LPXLOPER XL_Dates,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_CpnCurves,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas);



__declspec(dllexport) LPXLOPER WINAPI Local_PRDCCalculator_GetData(
    LPXLOPER XL_prdcId,
    LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRDCCalculator_GetData(
    LPXLOPER XL_prdcId,
    LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_payOrRec,
	LPXLOPER XL_coupon,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnCurves,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CallableSBCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_payOrRec,
	LPXLOPER XL_coupon,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnCurves,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas);


__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_SetData(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_Update(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CallableSBCalculator_SetData(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_GetPricingData(
	LPXLOPER XL_CallableSBCalculatorId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_CallOnMepiVanillaArg_Create(
	LPXLOPER XL_CurveName,
	LPXLOPER XL_EquityModelName,
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_RiskFactor,
	LPXLOPER XL_Strike,
	LPXLOPER XL_MaxBorrow,
	LPXLOPER XL_ProtectionCurveStart,
	LPXLOPER XL_ProtectionCurveEnd,
	LPXLOPER XL_StartingPf,
	LPXLOPER XL_StartingCash,
	LPXLOPER XL_MinInvested,
	LPXLOPER XL_LvgCost,
	LPXLOPER XL_CashSpread,
	LPXLOPER XL_Fees );


__declspec(dllexport) LPXLOPER WINAPI Local_CSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Leverage,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Leverage,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas);

__declspec(dllexport) LPXLOPER WINAPI Local_ExtendedCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_LeverageLong,
	LPXLOPER XL_LeverageShort,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
    LPXLOPER XL_fundingNotional);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ExtendedCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_LeverageLong,
	LPXLOPER XL_LeverageShort,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
    LPXLOPER XL_fundingNotional);



__declspec(dllexport) LPXLOPER WINAPI Local_BasicCSOCalculator_Create(LPXLOPER XL_startDate,
                                                                      LPXLOPER XL_endDate,
	                                                                  LPXLOPER XL_couponStructure,
	                                                                  LPXLOPER XL_fundingStructure,
	                                                                  LPXLOPER XL_exerciseStructure,
	                                                                  LPXLOPER XL_Notional,
	                                                                  LPXLOPER XL_CpnMin,
	                                                                  LPXLOPER XL_CpnMax,
	                                                                  LPXLOPER XL_LeverageLong,
	                                                                  LPXLOPER XL_LeverageShort,
	                                                                  LPXLOPER XL_Strike,
	                                                                  LPXLOPER XL_FixCpnCurve,
	                                                                  LPXLOPER XL_FundSpread,
	                                                                  LPXLOPER XL_fees,
	                                                                  LPXLOPER XL_FundLeverage,
                                                                      LPXLOPER XL_CpnCcy,
                                                                      LPXLOPER XL_FundCcy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasicCSOCalculator_Create(LPXLOPER XL_startDate,
																		  LPXLOPER XL_endDate,
																		  LPXLOPER XL_couponStructure,
																		  LPXLOPER XL_fundingStructure,
																		  LPXLOPER XL_exerciseStructure,
																		  LPXLOPER XL_Notional,
																		  LPXLOPER XL_CpnMin,
																		  LPXLOPER XL_CpnMax,
																		  LPXLOPER XL_LeverageLong,
																		  LPXLOPER XL_LeverageShort,
																		  LPXLOPER XL_Strike,
																		  LPXLOPER XL_FixCpnCurve,
																		  LPXLOPER XL_FundSpread,
																		  LPXLOPER XL_fees,
																		  LPXLOPER XL_FundLeverage,
																		  LPXLOPER XL_CpnCcy,
																		  LPXLOPER XL_FundCcy);

	/// it is ugly FIX FIX FIX, we have to merge all interfaces.
__declspec(dllexport) LPXLOPER WINAPI Local_BasisCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Cpns,
	LPXLOPER XL_Leverages,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas
    );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasisCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Cpns,
	LPXLOPER XL_Leverages,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas
    );

__declspec(dllexport) LPXLOPER WINAPI Local_CSOCalculator_GetData(
	LPXLOPER XL_csoId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CSOCalculator_GetData(
	LPXLOPER XL_csoId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_Update(
	LPXLOPER XL_genId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CFMethod_Create();
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CFMethod_Create();


__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_GetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_SetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_GetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_SetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_RootMrs(
	LPXLOPER XL_bsId,
	LPXLOPER XL_targetPrice,
	LPXLOPER XL_fTolerance,
	LPXLOPER XL_maxIer);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_RootMrs(
	LPXLOPER XL_bsId,
	LPXLOPER XL_targetPrice,
	LPXLOPER XL_fTolerance,
	LPXLOPER XL_maxIter);

__declspec(dllexport) LPXLOPER WINAPI Local_GetWarning_OnObject( 
	LPXLOPER XL_ObjId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetWarning_OnObject(
	LPXLOPER XL_ObjId );


__declspec(dllexport) LPXLOPER WINAPI Local_CRALocalCalculator_Create(
			LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRALocalCalculator_Create(
    		LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib);

__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_Create(
			LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_Create(
    		LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib);

__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_CreateFromPf(
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_isSfrmStdCalib);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_CreateFromPf(
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_isSfrmStdCalib);

__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_LocalCRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_LocalCRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalCRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalCRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys,
	LPXLOPER XL_update);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys,
	LPXLOPER XL_update);

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_CRAVMSCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refTenorCurve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRAVMSCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refTenorCurve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas,
	LPXLOPER XL_productsToPrice);


__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_currency,
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_ModelDatas,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_refResetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_currency,
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_ModelDatas,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_refResetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_CreateFromSwaption(
			LPXLOPER XL_currency,
			LPXLOPER XL_swaption,
			LPXLOPER XL_ModelDatas,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_refResetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_CreateFromSwaption(
			LPXLOPER XL_currency,
			LPXLOPER XL_swaption,
			LPXLOPER XL_ModelDatas,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_refResetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_CreateWithoutMarketData(LPXLOPER XL_generalDatas,
																					   LPXLOPER XL_callDatas,
																					   LPXLOPER XL_fundDatas,
																					   LPXLOPER XL_cpnDatas,
																					   LPXLOPER XL_notionalCurve,
																					   LPXLOPER XL_callFeesCurve,
																					   LPXLOPER XL_fundSpreadCurve,
																					   LPXLOPER XL_boostedFixCurve,
																					   LPXLOPER XL_payIndexMultCurve,
																					   LPXLOPER XL_barrierDownCurve,
																					   LPXLOPER XL_barrierUpCurve,
																					   LPXLOPER XL_refCoeff1Curve,
																					   LPXLOPER XL_refCoeff2Curve);

__declspec(dllexport) LPXLOPER WINAPI Local_CRADoubleCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_rateBarrierDownCurve,
	LPXLOPER XL_rateBarrierUpCurve,
	LPXLOPER XL_spreadBarrierDownCurve,
	LPXLOPER XL_spreadBarrierUpCurve,
	LPXLOPER XL_spreadCoeff1,
	LPXLOPER XL_spreadCoeff2,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRADoubleCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_rateBarrierDownCurve,
	LPXLOPER XL_rateBarrierUpCurve,
	LPXLOPER XL_spreadBarrierDownCurve,
	LPXLOPER XL_spreadBarrierUpCurve,
	LPXLOPER XL_spreadCoeff1,
	LPXLOPER XL_spreadCoeff2,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_BasicCRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_currency,
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_refResetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasicCRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_currency,
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_refResetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundData,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_globalCapParams,
			LPXLOPER XL_pastFixings,
			LPXLOPER XL_capLevCurve,
			LPXLOPER XL_capFixedCurve,
			LPXLOPER XL_capStrikeCurve,
			LPXLOPER XL_capSpreadCurve,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_CreateFromPf(
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve,			
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_CreateFromPfWithoutMktData(
			LPXLOPER XL_asOf,
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve);

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_SetData(
			LPXLOPER XL_calculatorId,
			LPXLOPER XL_dataId,
			LPXLOPER XL_setPortfolioType,
			LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GlobalCapCalculator_SetData(
			LPXLOPER XL_calculatorId,
			LPXLOPER XL_dataId,
			LPXLOPER XL_setPortfolioType,
			LPXLOPER XL_mktDataKeys);

__declspec(dllexport) LPXLOPER WINAPI Local_SnowRangeCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundingData,
			LPXLOPER XL_couponData,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_strikeCurve,			
			LPXLOPER XL_ratchetCurve,			
			LPXLOPER XL_cashFlowCurve,
			LPXLOPER XL_fixedRateCurve,
			LPXLOPER XL_leverageCurve,
			LPXLOPER XL_snowRangeParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_mcParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SnowRangeCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundingData,
			LPXLOPER XL_couponData,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_strikeCurve,			
			LPXLOPER XL_ratchetCurve,			
			LPXLOPER XL_cashFlowCurve,
			LPXLOPER XL_fixedRateCurve,
			LPXLOPER XL_leverageCurve,
			LPXLOPER XL_snowRangeParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_mcParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_CMRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_exerProbas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CMRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_exerProbas,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_TarnFxCalculator_Create(LPXLOPER XL_StartDate,
																	LPXLOPER XL_EndDate,
																	LPXLOPER XL_Currencies,
																	LPXLOPER XL_PayRec,
																	LPXLOPER XL_CpnData,
																	LPXLOPER XL_CpnNominal,
																	LPXLOPER XL_DomCoupon,			
																	LPXLOPER XL_ForCoupon,			
																	LPXLOPER XL_InitialFX,
																	LPXLOPER XL_MinCouponCurve,
																	LPXLOPER XL_MaxCouponCurve,
																	LPXLOPER XL_FundingData,
																	LPXLOPER XL_FundingSpread,
																	LPXLOPER XL_FundingNominal,
																	LPXLOPER XL_RedemptionData,
																	LPXLOPER XL_Fees,
																	LPXLOPER XL_MktDataManager,
																	LPXLOPER XL_ProductsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TarnFxCalculator_Create(LPXLOPER XL_StartDate,
																		LPXLOPER XL_EndDate,
																		LPXLOPER XL_Currencies,
																		LPXLOPER XL_PayRec,
																		LPXLOPER XL_CpnData,
																		LPXLOPER XL_CpnNominal,
																		LPXLOPER XL_DomCoupon,			
																		LPXLOPER XL_ForCoupon,																		
																		LPXLOPER XL_InitialFX,
																		LPXLOPER XL_MinCouponCurve,
																		LPXLOPER XL_MaxCouponCurve,
																		LPXLOPER XL_FundingData,
																		LPXLOPER XL_FundingSpread,
																		LPXLOPER XL_FundingNominal,
																		LPXLOPER XL_RedemptionData,
																		LPXLOPER XL_Fees,
																		LPXLOPER XL_MktDataManager,
																		LPXLOPER XL_ProductsToPrice);
__declspec(dllexport) LPXLOPER WINAPI Local_ConvertToVarNotionalSwaption(
			LPXLOPER XL_yieldCurve,
			LPXLOPER XL_swaption);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ConvertToVarNotionalSwaption(
			LPXLOPER XL_yieldCurve,
			LPXLOPER XL_swaption);

__declspec(dllexport) LPXLOPER WINAPI Local_CRAQuantoCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_fixCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRAQuantoCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_fixCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_productsToPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_FXVanillaCalculator_CreateFromSecurity(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_BasketType,
	LPXLOPER XL_DigitType,
    LPXLOPER XL_VanillaType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FXVanillaCalculator_CreateFromSecurity(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_BasketType,
	LPXLOPER XL_DigitType,
    LPXLOPER XL_VanillaType);

__declspec(dllexport) LPXLOPER WINAPI Local_VolBondCalculator_Create(
	LPXLOPER XL_NominalStartEndDate,
	LPXLOPER XL_PayFreq,
	LPXLOPER XL_ResetFreq,	
	LPXLOPER XL_DayCount,
	LPXLOPER XL_Tenor,	
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap,
	LPXLOPER XL_PayCalendar,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_OdeSolvers,
	LPXLOPER XL_RKParameters,	
	LPXLOPER XL_MCParameters,
	LPXLOPER XL_RandomGenerator,
	LPXLOPER XL_PayOffType,	
	LPXLOPER XL_MarketDataManager,
	LPXLOPER XL_MarketDataManagerKeys,
	LPXLOPER XL_ProductsToPrice);

#endif
