#ifndef ARM_XL_SWTION_LOCAL_H
#define ARM_XL_SWTION_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_SWAPTION (LPXLOPER XL_swapId,
													  LPXLOPER XL_receiveOrPay,
													  LPXLOPER XL_strike,
													  LPXLOPER XL_maturity,
													  LPXLOPER XL_exerciseType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SWAPTION (LPXLOPER XL_swapId,
														  LPXLOPER XL_receiveOrPay,
														  LPXLOPER XL_strike,
														  LPXLOPER XL_maturity,
														  LPXLOPER XL_exerciseType);

__declspec(dllexport) LPXLOPER WINAPI Local_SwaptionFromExpiry(LPXLOPER XL_optionExpiry,
													           LPXLOPER XL_swapTerm,
													           LPXLOPER XL_liborType,
													           LPXLOPER XL_strike,
													           LPXLOPER XL_receiveOrPay,
													           LPXLOPER XL_spread,
													           LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SwaptionFromExpiry(LPXLOPER XL_optionExpiry,
													           LPXLOPER XL_swapTerm,
													           LPXLOPER XL_liborType,
													           LPXLOPER XL_strike,
													           LPXLOPER XL_receiveOrPay,
													           LPXLOPER XL_spread,
													           LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_LIBORSWAPTION (LPXLOPER XL_startDate,
														   LPXLOPER XL_endDate,
														   LPXLOPER XL_receiveOrPay,
														   LPXLOPER XL_strike,
														   LPXLOPER XL_maturity,
														   LPXLOPER XL_liborType,
														   LPXLOPER XL_spread,
														   LPXLOPER XL_exerciseType,
														   LPXLOPER XL_resetFreq,
														   LPXLOPER XL_payFreq,
														   LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORSWAPTION (LPXLOPER XL_startDate,
															   LPXLOPER XL_endDate,
															   LPXLOPER XL_receiveOrPay,
															   LPXLOPER XL_strike,
															   LPXLOPER XL_maturity,
															   LPXLOPER XL_liborType,
															   LPXLOPER XL_spread,
															   LPXLOPER XL_exerciseType,
															   LPXLOPER XL_resetFreq,
															   LPXLOPER XL_payFreq,
															   LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_EXOSWAPTION (LPXLOPER XL_swapId,
														 LPXLOPER XL_receiveOrPay,
														 LPXLOPER XL_xStyleId,
														 LPXLOPER XL_kRefValId,
														 LPXLOPER XL_swapYearTerm);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EXOSWAPTION (LPXLOPER XL_swapId,
															 LPXLOPER XL_receiveOrPay,
															 LPXLOPER XL_xStyleId,
															 LPXLOPER XL_kRefValId,
															 LPXLOPER XL_swapYearTerm);

__declspec(dllexport) LPXLOPER WINAPI Local_VARFIXSWAPTION (LPXLOPER XL_startDate,
															LPXLOPER XL_endDate,
															LPXLOPER XL_fixSpreads,
															LPXLOPER XL_XStyle,
															LPXLOPER XL_receiveOrPay,
															LPXLOPER XL_strike,
															LPXLOPER XL_maturity,
															LPXLOPER XL_liborType,
															LPXLOPER XL_spread,
															LPXLOPER XL_resetFreq,
															LPXLOPER XL_payFreq,
															LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VARFIXSWAPTION (LPXLOPER XL_startDate,
																LPXLOPER XL_endDate,
																LPXLOPER XL_fixSpreads,
																LPXLOPER XL_XStyle,
																LPXLOPER XL_receiveOrPay,
																LPXLOPER XL_strike,
																LPXLOPER XL_maturity,
																LPXLOPER XL_liborType,
																LPXLOPER XL_spread,
																LPXLOPER XL_resetFreq,
																LPXLOPER XL_payFreq,
																LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_OPTIONALACCRUALZCBOND(LPXLOPER XL_startDate,
																	  LPXLOPER XL_endDate,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_payfreq,
																	  LPXLOPER XL_nbCurPerforAcc,
																	  LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_OPTIONALACCRUALZCBOND(LPXLOPER XL_startDate,
																		  LPXLOPER XL_endDate,
																		  LPXLOPER XL_strike,
																		  LPXLOPER XL_payfreq,
																		  LPXLOPER XL_nbCurPerforAcc,
																		  LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_EXOCFSWAPTION (LPXLOPER XL_swapId,
														   LPXLOPER XL_receiveOrPay,
														   LPXLOPER XL_capOrFloor,
														   LPXLOPER XL_xStyleId,
														   LPXLOPER XL_kSptionRefValId,
														   LPXLOPER XL_kCFloorRefValId,
														   LPXLOPER XL_floorPosition,
														   LPXLOPER XL_IsBarrierCF);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EXOCFSWAPTION (LPXLOPER XL_swapId,
															   LPXLOPER XL_receiveOrPay,
															   LPXLOPER XL_capOrFloor,
															   LPXLOPER XL_xStyleId,
															   LPXLOPER XL_kSptionRefValId,
															   LPXLOPER XL_kCFloorRefValId,
															   LPXLOPER XL_floorPosition,
															   LPXLOPER XL_IsBarrierCF);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FlexAccretSwaption (LPXLOPER XL_startDate,
																	LPXLOPER XL_endDate,
																	LPXLOPER XL_fixedRate,
																	LPXLOPER XL_nbCurrentPeriodsForAccrued,
																	LPXLOPER XL_receiveOrPay,
																	LPXLOPER XL_freq,
																	LPXLOPER XL_liborType,
																	LPXLOPER XL_spread,
																	LPXLOPER XL_exerciseDates,
																	LPXLOPER XL_currency);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FlexAccretSwaption (LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
																		LPXLOPER XL_fixedRate,
																		LPXLOPER XL_nbCurrentPeriodsForAccrued,
																		LPXLOPER XL_receiveOrPay,
																		LPXLOPER XL_freq,
																		LPXLOPER XL_liborType,
																		LPXLOPER XL_spread,
																		LPXLOPER XL_exerciseDates,
																		LPXLOPER XL_currency);

__declspec(dllexport) LPXLOPER WINAPI Local_SwaptionStickyDelta(LPXLOPER XL_swaptionId,
																LPXLOPER XL_modelId,
																LPXLOPER XL_perturbeDiscountCurvId);

#endif /* ARM_XL_SWTION_LOCAL_H */