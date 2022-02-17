#ifndef ARM_XL_OPTION_LOCAL_H
#define ARM_XL_OPTION_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_EXOPTION (LPXLOPER XL_underlyingId,
											    LPXLOPER XL_optionType,
											    LPXLOPER XL_styleId,
											    LPXLOPER XL_KRefValId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EXOPTION (LPXLOPER XL_underlyingId,
												LPXLOPER XL_optionType,
												LPXLOPER XL_styleId,
												LPXLOPER XL_KRefValId);

__declspec(dllexport) LPXLOPER WINAPI Local_OPTION (LPXLOPER XL_underlyingId,
											  LPXLOPER XL_maturity,
											  LPXLOPER XL_strike,
											  LPXLOPER XL_optionType,
											  LPXLOPER XL_exerciseType,
											  LPXLOPER XL_strikeType,
											  LPXLOPER XL_FstXDate,
											  LPXLOPER XL_PayDate);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_OPTION (LPXLOPER XL_underlyingId,
											  LPXLOPER XL_maturity,
											  LPXLOPER XL_strike,
											  LPXLOPER XL_optionType,
											  LPXLOPER XL_exerciseType,
											  LPXLOPER XL_strikeType,
											  LPXLOPER XL_FstXDate,
											  LPXLOPER XL_PayDate);

__declspec(dllexport) LPXLOPER WINAPI Local_VolImp (LPXLOPER XL_secId,
													LPXLOPER XL_modId,
													LPXLOPER XL_price);

__declspec(dllexport) LPXLOPER WINAPI Local_bsOption (LPXLOPER XL_spot,
													  LPXLOPER XL_strike,
													  LPXLOPER XL_volatility,
													  LPXLOPER XL_dividend,
													  LPXLOPER XL_discountRate,
													  LPXLOPER XL_maturity,
													  LPXLOPER XL_CallPut);

__declspec(dllexport) LPXLOPER WINAPI Local_bsDelta (LPXLOPER XL_spot,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_volatility,
													 LPXLOPER XL_dividend,
													 LPXLOPER XL_discountRate,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_CallPut);

__declspec(dllexport) LPXLOPER WINAPI Local_bsVega (LPXLOPER XL_spot,
													LPXLOPER XL_strike,
													LPXLOPER XL_volatility,
													LPXLOPER XL_dividend,
													LPXLOPER XL_discountRate,
													LPXLOPER XL_maturity,
													LPXLOPER XL_CallPut);

__declspec(dllexport) LPXLOPER WINAPI Local_bsGamma (LPXLOPER XL_spot,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_volatility,
													 LPXLOPER XL_dividend,
													 LPXLOPER XL_discountRate,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_CallPut);

__declspec(dllexport) LPXLOPER WINAPI Local_bsTheta (LPXLOPER XL_spot,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_volatility,
													 LPXLOPER XL_dividend,
													 LPXLOPER XL_discountRate,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_CallPut);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_OPTIONONPORTFOLIO (LPXLOPER XL_portfolioId,
																   LPXLOPER XL_styleId,
																   LPXLOPER XL_strikesId,
																   LPXLOPER XL_optionType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_OPTIONONPORTFOLIO (LPXLOPER XL_portfolioId,
																	   LPXLOPER XL_styleId,
																	   LPXLOPER XL_strikesId,
																	   LPXLOPER XL_optionType);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LPXLOPER XL_CalibratorId,
															LPXLOPER XL_PfTypeId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LPXLOPER XL_CalibratorId,
															LPXLOPER XL_PfTypeId);


_declspec(dllexport) LPXLOPER WINAPI Local_SpreadOptionFormula(
	LPXLOPER XL_fwd1, 
	LPXLOPER XL_fwd2, 
	LPXLOPER XL_vol1, 
	LPXLOPER XL_vol2, 
	LPXLOPER XL_Correl, 
	LPXLOPER XL_strike, 
	LPXLOPER XL_optMat, 
	LPXLOPER XL_optType, 
	LPXLOPER XL_modelType, 
	LPXLOPER XL_spreadVol );

_declspec(dllexport) LPXLOPER WINAPI Local_SumOption(
	LPXLOPER XL_dateStrip,
	LPXLOPER XL_payDate,
	LPXLOPER XL_strike,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_indexDayCount,
	LPXLOPER XL_coeff );

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_SumOption(
	LPXLOPER XL_dateStrip,
	LPXLOPER XL_payDate,
	LPXLOPER XL_strike,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_indexDayCount,
	LPXLOPER XL_coeff );

__declspec(dllexport) LPXLOPER WINAPI Local_STRIPOPTION (LPXLOPER XL_underlyingId,
														 LPXLOPER XL_optionType,
														 LPXLOPER XL_strikesId,
														 LPXLOPER XL_scheduleId,
														 LPXLOPER XL_PorS,
														 LPXLOPER XL_fxFixingsId,
														 LPXLOPER XL_leverageId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STRIPOPTION (LPXLOPER XL_underlyingId,
															 LPXLOPER XL_optionType,
															 LPXLOPER XL_strikesId,
															 LPXLOPER XL_scheduleId,
															 LPXLOPER XL_PorS,
															 LPXLOPER XL_fxFixingsId,
															 LPXLOPER XL_leverageId);

__declspec(dllexport) LPXLOPER WINAPI Local_STRIPDIGITALOPTION ( LPXLOPER XL_underlyingId,
																 LPXLOPER XL_optionType,
																 LPXLOPER XL_strikesId,
																 LPXLOPER XL_scheduleId,
																 LPXLOPER XL_correl,
																 LPXLOPER XL_PorS,
																 LPXLOPER XL_fxFixingsId,
																 LPXLOPER XL_digitalOptionData,
																 LPXLOPER XL_leverageId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STRIPDIGITALOPTION ( LPXLOPER XL_underlyingId,
																	 LPXLOPER XL_optionType,
																	 LPXLOPER XL_strikesId,
																	 LPXLOPER XL_scheduleId,
																	 LPXLOPER XL_correl,
																	 LPXLOPER XL_PorS,
																	 LPXLOPER XL_fxFixingsId,
																	 LPXLOPER XL_digitalOptionData,
																	 LPXLOPER XL_leverageId);

__declspec(dllexport) LPXLOPER WINAPI Local_FxOptionStrip(LPXLOPER XL_underlyingId, 
														  LPXLOPER XL_strikesCurveId, 
														  LPXLOPER XL_optionType, 
														  LPXLOPER XL_startDate, 
														  LPXLOPER XL_endDate, 
														  LPXLOPER XL_notionalId,
														  LPXLOPER XL_resetData,
														  LPXLOPER XL_payData,
														  LPXLOPER XL_dayCount, 
														  LPXLOPER XL_fwdRule, 
														  LPXLOPER XL_intRule, 
														  LPXLOPER XL_stubRule,
														  LPXLOPER XL_PorS,
														  LPXLOPER XL_fxFixingsId,
														  LPXLOPER XL_digitalOptionData,
														  LPXLOPER XL_leverageId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FxOptionStrip(LPXLOPER XL_underlyingId, 
															  LPXLOPER XL_strikesCurveId, 
															  LPXLOPER XL_optionType, 
															  LPXLOPER XL_startDate, 
															  LPXLOPER XL_endDate, 
															  LPXLOPER XL_notionalId,
															  LPXLOPER XL_resetData,
															  LPXLOPER XL_payData,
															  LPXLOPER XL_dayCount, 
															  LPXLOPER XL_fwdRule, 
															  LPXLOPER XL_intRule, 
															  LPXLOPER XL_stubRule,
															  LPXLOPER XL_PorS,
															  LPXLOPER XL_fxFixingsId,
															  LPXLOPER XL_digitalOptionData,
															  LPXLOPER XL_leverageId);

#endif /* ARM_XL_OPTION_LOCAL_H */
