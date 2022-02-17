#ifndef ARM_XL_GLOB_LOCAL_H
#define ARM_XL_GLOB_LOCAL_H


__declspec(dllexport) LPXLOPER WINAPI Local_bsflexible (LPXLOPER XL_F,LPXLOPER XL_V,
														LPXLOPER XL_B,LPXLOPER XL_K,
													   LPXLOPER XL_CallPutFlag);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Price (	LPXLOPER XL_secId,
															LPXLOPER XL_modId);

__declspec(dllexport) LPXLOPER WINAPI Local_FreeObject (LPXLOPER XL_secId);

__declspec(dllexport) LPXLOPER WINAPI Local_FreeAllObjects ();

__declspec(dllexport) LPXLOPER WINAPI Local_NextBusinessDay (LPXLOPER XL_date,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_days);

__declspec(dllexport) LPXLOPER WINAPI Local_IsBusinessDay (LPXLOPER XL_date,
														   LPXLOPER XL_currency);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Accrued (LPXLOPER XL_secId,
														 LPXLOPER XL_fwdDate,
														 LPXLOPER XL_modId);

__declspec(dllexport) LPXLOPER WINAPI Local_IMPLIEDVOL (LPXLOPER XL_instId,
														LPXLOPER XL_modId,
														LPXLOPER XL_price,
														LPXLOPER XL_LnOrNorVol);

__declspec(dllexport) LPXLOPER WINAPI Local_SetNotional (LPXLOPER XL_secId,
														 LPXLOPER XL_refValId,
														 LPXLOPER XL_percentRemainder);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_View (LPXLOPER XL_instId);

__declspec(dllexport) long WINAPI Local_View ();

__declspec(dllexport) long WINAPI Local_View_XML();

__declspec(dllexport) int WINAPI Local_ARM_Help (void);

extern long Local_ViewFile (const CCString& C_instId);

extern long LocalGetArmViewFile (const CCString& sockId);

__declspec(dllexport) int WINAPI Local_DisplayErrorMessage ();

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_today ();

__declspec(dllexport) LPXLOPER WINAPI Local_GetExpiry (LPXLOPER XL_secId);

__declspec(dllexport) LPXLOPER WINAPI Local_Sensitivity (LPXLOPER XL_secId,
														 LPXLOPER XL_modId,
														 LPXLOPER XL_param);

__declspec(dllexport) LPXLOPER WINAPI Local_GetFRMShortRateVols(LPXLOPER XL_modelId);

__declspec(dllexport) LPXLOPER WINAPI Local_XCcyAdjust (LPXLOPER XL_startDate,
														LPXLOPER XL_endDate, 
														LPXLOPER XL_payFreq,
														LPXLOPER XL_domCcy,
														LPXLOPER XL_forIndexType,
														LPXLOPER XL_forCcy,
														LPXLOPER XL_spreadsId,
														LPXLOPER XL_zcDomId,
														LPXLOPER XL_discDomId,
														LPXLOPER XL_zcForId,
														LPXLOPER XL_discForId,
														LPXLOPER XL_FX,
														LPXLOPER XL_CouponId,
														LPXLOPER XL_domDc,
														LPXLOPER XL_forDc);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_XCcyAdjust (LPXLOPER XL_startDate,
															LPXLOPER XL_endDate, 
															LPXLOPER XL_payFreq,
															LPXLOPER XL_domCcy,
															LPXLOPER XL_forIndexType,
															LPXLOPER XL_forCcy,
															LPXLOPER XL_spreadsId,
															LPXLOPER XL_zcDomId,
															LPXLOPER XL_discDomId,
															LPXLOPER XL_zcForId,
															LPXLOPER XL_discForId,
															LPXLOPER XL_FX,
															LPXLOPER XL_CouponId,
															LPXLOPER XL_domDc,
															LPXLOPER XL_forDc);

__declspec(dllexport) LPXLOPER WINAPI Local_ADJUSTTOBUSDATE (LPXLOPER XL_date,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_rule);

__declspec(dllexport) LPXLOPER WINAPI Local_FwdPrice (LPXLOPER XL_secId,
													  LPXLOPER XL_modId,
													  LPXLOPER XL_fwdDate);

__declspec(dllexport) LPXLOPER WINAPI Local_CvSensitivity (LPXLOPER XL_secId,
														   LPXLOPER XL_modId,
														   LPXLOPER XL_param);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BetweenDates (LPXLOPER XL_date1,
															  LPXLOPER XL_date2,
															  LPXLOPER XL_daycount,
															  LPXLOPER XL_isYearFrac);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CountBusinessDays (LPXLOPER XL_date1,
																   LPXLOPER XL_date2,
																   LPXLOPER XL_calendar);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ADDMONTHS (LPXLOPER XL_date,
														   LPXLOPER XL_nb,
														   LPXLOPER XL_rule,
														   LPXLOPER XL_currency);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ADDYEARS (LPXLOPER XL_date,
														  LPXLOPER XL_nb,
														  LPXLOPER XL_rule,
														  LPXLOPER XL_currency);

__declspec(dllexport) LPXLOPER WINAPI Local_FxConvert (LPXLOPER XL_ccy1,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_asOfDate,
													   LPXLOPER XL_amount,
													   LPXLOPER XL_cvname);

__declspec(dllexport) LPXLOPER WINAPI Local_FxConvertFromCalypso (LPXLOPER XL_ccy1,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_asOfDate,
													   LPXLOPER XL_cvname);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCurrency(LPXLOPER XL_Security);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetDefaultCurrency ();

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetDefaultCurrency (LPXLOPER XL_currency);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ADDPERIOD (LPXLOPER XL_date,
														   LPXLOPER XL_freq,
														   LPXLOPER XL_ccy,
														   LPXLOPER XL_nbPeriods,
														   LPXLOPER XL_adjRule,
														   LPXLOPER XL_goToEndOfMonth);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Cover (LPXLOPER XL_secId,
													   LPXLOPER XL_modId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Price_OptUnder (LPXLOPER XL_secId,
																LPXLOPER XL_modId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetPID(void);

__declspec(dllexport) LPXLOPER WINAPI Local_ParallelShift (LPXLOPER XL_curveId,
														   LPXLOPER XL_value);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ParallelShift (LPXLOPER XL_curveId,
															   LPXLOPER XL_value);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ClonedAndSetNotional (LPXLOPER XL_secId,
																	  LPXLOPER XL_refValId,
																	  LPXLOPER XL_percentRemainder);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_ClonedAndSetNotional (LPXLOPER XL_secId,
																		  LPXLOPER XL_refValId,
																		  LPXLOPER XL_percentRemainder);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_Clone (LPXLOPER XL_objectId);

__declspec(dllexport) LPXLOPER WINAPI Local_FIXRATES (LPXLOPER XL_secId,
													  LPXLOPER XL_rate);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FIXRATES (LPXLOPER XL_secId,
														  LPXLOPER XL_rate);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayScheduleValues (LPXLOPER XL_instId,
																	   LPXLOPER XL_typeValues,
																	   LPXLOPER XL_RecOrPay,
																	   LPXLOPER XL_modelId);
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayScheduleDates (LPXLOPER XL_instId,
																	  LPXLOPER XL_typeDates,
																	  LPXLOPER XL_RecOrPay,
																	  LPXLOPER XL_viewInitExch);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayReplicPort(LPXLOPER XL_instId,
                                                                  LPXLOPER XL_WeightOrStrike,
                                                                  LPXLOPER XL_PayoffOrSensi,
                                                                  LPXLOPER XL_RecOrPay,
                                                                  LPXLOPER XL_ModelId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INTERPOL (LPXLOPER XL_vexX,
														  LPXLOPER XL_vecY,
														  LPXLOPER XL_X,
														  LPXLOPER XL_typeInterpol);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIANGULARINTERPOL (LPXLOPER XL_vexX,
																	LPXLOPER XL_vecY,
																	LPXLOPER XL_matZ,
																	LPXLOPER XL_X,
																	LPXLOPER XL_Y);

__declspec(dllexport) LPXLOPER WINAPI Local_KImp (LPXLOPER XL_sec,
												  LPXLOPER XL_model,
												  LPXLOPER XL_price,
												  LPXLOPER XL_param);

__declspec(dllexport) LPXLOPER WINAPI Local_BSSpot (LPXLOPER XL_secId,
													LPXLOPER XL_modId,
													LPXLOPER XL_date);

__declspec(dllexport) void WINAPI Local_ARM_ProdConnect(void);

__declspec(dllexport) void WINAPI Local_ARM_RepliConnect(void);

__declspec(dllexport) void WINAPI Local_ARM_InfocConnect(void);

__declspec(dllexport) void WINAPI Local_ARM_RecConnect(void);

__declspec(dllexport) void WINAPI Local_ARM_ShutDownETK(void);

__declspec(dllexport) void WINAPI Local_ARM_SwitchToETK(void);

__declspec(dllexport) void WINAPI Local_ARM_SwitchToFLATFILE(void);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetMeanRevFromSummit(LPXLOPER XL_ccy,
																	 LPXLOPER XL_index,
																	 LPXLOPER XL_cvname,
																	 LPXLOPER XL_date,
																	 LPXLOPER XL_2or3Factor);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCutOffFromSummit(LPXLOPER XL_ccy,
																	LPXLOPER XL_index,
																	LPXLOPER XL_cvname,
																	LPXLOPER XL_date);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Hedge(LPXLOPER XL_sec,
													  LPXLOPER XL_typeHedge);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETINFOFROMPRCS(LPXLOPER XL_prcs,
																LPXLOPER XL_datatype);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETOBJINFOFROMPRCS(LPXLOPER XL_prcs,
																   LPXLOPER XL_datatype);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFixing(LPXLOPER XL_source,
														  LPXLOPER XL_index,
														  LPXLOPER XL_term,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_date);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFixingFromCalypso(
														  LPXLOPER XL_index,
														  LPXLOPER XL_term,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_source,
														  LPXLOPER XL_cvname,
														  LPXLOPER XL_date);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetSecurityData(
	LPXLOPER XL_secId,
	LPXLOPER XL_data);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SetSecurityData(
	LPXLOPER XL_secId,
	LPXLOPER XL_data);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetSecurityData(
	LPXLOPER XL_secId,
	LPXLOPER XL_data);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetDomCcy(LPXLOPER XL_calcId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetDomCcy(LPXLOPER XL_calcId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetForCcy(LPXLOPER XL_calcId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetForCcy(LPXLOPER XL_calcId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFundingCcy(LPXLOPER XL_calcId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetFundingCcy(LPXLOPER XL_calcId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SecurityFlows(LPXLOPER XL_labels,LPXLOPER XL_values);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ViewCell(LPXLOPER XL_ObjectId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_MatrixVectorViewer(LPXLOPER XL_ObjectId);

__declspec(dllexport) LPXLOPER WINAPI Local_GetFixingFromInstrument(LPXLOPER XL_tradeId);

#endif	/* ARM_XL_GLOB_LOCAL_H */
