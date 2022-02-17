/*----------------------------------------------------------------------------------*/

#ifndef _ARM_XL_DKPRCS_LOCAl_H
#define _ARM_XL_DKPRCS_LOCAl_H



__declspec(dllexport) LPXLOPER WINAPI Local_PRCS3F_LatticePricing(LPXLOPER XL_dLatticeGeometryDataIn,
																  LPXLOPER XL_dNumTimeLinesBeforeFirstNotice,
																  LPXLOPER XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
																  LPXLOPER XL_dNumTimeLines, 
																  LPXLOPER XL_evalDate,
																  /*
																  LPXLOPER XL_dBaseYieldCurveId,
																  LPXLOPER XL_dForeignYieldCurveId,
																  LPXLOPER XL_volSwopBaseId,
																  LPXLOPER XL_volSwopForeignId,
																  LPXLOPER XL_volFxId,
																  */
																  LPXLOPER XL_curves,
																  LPXLOPER XL_dNoticeDatesIn,
																  LPXLOPER XL_dStrike,
																  LPXLOPER XL_dType,
																  LPXLOPER XL_dOptionExpiry, 
																  LPXLOPER XL_dMeanReversionBase,
																  LPXLOPER XL_dMeanReversionForeign,  
																  LPXLOPER XL_dSpotFX,
																  LPXLOPER XL_dBaseForeignCorrelationId,
																  LPXLOPER XL_dBaseSpotFXCorrelationId,
																  LPXLOPER XL_dForeignSpotFXCorrelationId,
																  LPXLOPER XL_dProductModelCodeId,
																  LPXLOPER XL_paraModel,
																  LPXLOPER XL_dBoosterDataIn);



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TREE3FACT(LPXLOPER XL_asof,
														  LPXLOPER XL_xbsfx,
														  LPXLOPER XL_volSwopBaseId,
														  LPXLOPER XL_volSwopForeignId,
														  LPXLOPER XL_dBaseForeignCorrelationId,
														  LPXLOPER XL_dLatticeGeometryDataIn,
														  LPXLOPER XL_dNumTimeLinesBeforeFirstNotice,
														  LPXLOPER XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
														  LPXLOPER XL_dNumTimeLines,
														  LPXLOPER XL_dMeanReversionBase,
														  LPXLOPER XL_dMeanReversionForeign,
														  LPXLOPER XL_limitVector,
														  LPXLOPER XL_dOptimal,
														  LPXLOPER XL_dTimeBoost,
														  LPXLOPER XL_dDeltaFlag,
														  LPXLOPER XL_dSmoothingValue,
														  LPXLOPER XL_sCalcProbaSurvOrNot,
														  LPXLOPER XL_convInVolFwd);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_TREE3FACT(LPXLOPER XL_asof,
															  LPXLOPER XL_xbsfx,
															  LPXLOPER XL_volSwopBaseId,
															  LPXLOPER XL_volSwopForeignId,
															  LPXLOPER XL_dBaseForeignCorrelationId,
															  LPXLOPER XL_dLatticeGeometryDataIn,
															  LPXLOPER XL_dNumTimeLinesBeforeFirstNotice,
															  LPXLOPER XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
															  LPXLOPER XL_dNumTimeLines,
															  LPXLOPER XL_dMeanReversionBase,
															  LPXLOPER XL_dMeanReversionForeign,
															  LPXLOPER XL_limitVector,
															  LPXLOPER XL_dOptimal,
															  LPXLOPER XL_dTimeBoost,
															  LPXLOPER XL_dDeltaFlag,
															  LPXLOPER XL_dSmoothingValue,
															  LPXLOPER XL_sCalcProbaSurvOrNot,
															  LPXLOPER XL_convInVolFwd);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Bootstrapping_VFDK_HW1To3F(LPXLOPER XL_volCurve,
                                                                           LPXLOPER XL_zc,
                                                                           LPXLOPER XL_noticeDates,
                                                                           LPXLOPER XL_swapStartDates,
                                                                           LPXLOPER XL_swapEndDates,
                                                                           LPXLOPER XL_HW3FParameters,
                                                                           LPXLOPER XL_observationDate);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SwaptionPrice_VFDK_HW1To3F(LPXLOPER XL_dSwaptionExpiryInYears,
                                                                           LPXLOPER XL_dSwaptionTenorInYears,
                                                                           LPXLOPER XL_dNoticePeriodInDays,
                                                                           LPXLOPER XL_dStrike,
                                                                           LPXLOPER XL_dCallPut,
                                                                           LPXLOPER XL_zc,
                                                                           LPXLOPER XL_noticeDates,
                                                                           LPXLOPER XL_sigma,
                                                                           LPXLOPER XL_HW3FParameters,
                                                                           LPXLOPER XL_observationDate);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ImpliedFwdCorrelation_VFDK_HW1To3F(LPXLOPER XL_dSwaptionExpiryInYears,
                                                                                   LPXLOPER XL_dSwaptionTenorInYears,
                                                                                   LPXLOPER XL_dSwaptionTenor2InYears,
                                                                                   LPXLOPER XL_dNoticePeriodInDays,
                                                                                   LPXLOPER XL_zc,
                                                                                   LPXLOPER XL_noticeDates,
                                                                                   LPXLOPER XL_sigma,
                                                                                   LPXLOPER XL_HW3FParameters,
                                                                                   LPXLOPER XL_observationDate);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_HW3F_CALIBRATION(   LPXLOPER XL_dSpotDate,
                                                                    LPXLOPER XL_zc,
                                                                    LPXLOPER XL_HW3FParametersIn,
                                                                    LPXLOPER XL_volCurve,
                                                                    LPXLOPER XL_correlationCurve,
                                                                    LPXLOPER XL_volWeightsCurve,
                                                                    LPXLOPER XL_correlationWeightsCurve,
                                                                    LPXLOPER XL_dNoticeDates,
                                                                    LPXLOPER XL_dSwapStartDates,
                                                                    LPXLOPER XL_dSwapEndDates);



#endif
/*----------------------------------------------------------------------------------*/
/*---- End Of File ----*/