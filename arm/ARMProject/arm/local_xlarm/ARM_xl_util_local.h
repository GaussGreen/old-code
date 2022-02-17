#ifndef ARM_XL_UTIL_LOCAL_H
#define ARM_XL_UTIL_LOCAL_H


__declspec(dllexport) LPXLOPER WINAPI Local_GetCorrelInst(LPXLOPER XL_date1,
														  LPXLOPER XL_date2,
														  LPXLOPER XL_ccy1,
														  LPXLOPER XL_index1,
														  LPXLOPER XL_expiry1,
														  LPXLOPER XL_tenor1,
														  LPXLOPER XL_curve1_ccy1,	
														  LPXLOPER XL_curve2_ccy1,	
														  LPXLOPER XL_nbmonths_curve1_ccy1,
														  LPXLOPER XL_ccy2,
														  LPXLOPER XL_index2,
														  LPXLOPER XL_expiry2,
														  LPXLOPER XL_tenor2,
														  LPXLOPER XL_curve1_ccy2,	
														  LPXLOPER XL_curve2_ccy2,	
														  LPXLOPER XL_nbmonths_curve1_ccy2,
														  LPXLOPER XL_type,
														  LPXLOPER XL_lambda,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_precision);

__declspec(dllexport) LPXLOPER WINAPI Local_GetMoyCorrel(LPXLOPER XL_date1,
														 LPXLOPER XL_date2,
														 LPXLOPER XL_ccy1,
														 LPXLOPER XL_index1,
														 LPXLOPER XL_expiry1,
														 LPXLOPER XL_tenor1,
														 LPXLOPER XL_curve1_ccy1,	
														 LPXLOPER XL_curve2_ccy1,	
														 LPXLOPER XL_nbmonths_curve1_ccy1,
														 LPXLOPER XL_ccy2,
														 LPXLOPER XL_index2,
														 LPXLOPER XL_expiry2,
														 LPXLOPER XL_tenor2,
														 LPXLOPER XL_curve1_ccy2,	
														 LPXLOPER XL_curve2_ccy2,	
														 LPXLOPER XL_nbmonths_curve1_ccy2,
														 LPXLOPER XL_type,
														 LPXLOPER XL_lambda,
														 LPXLOPER XL_ccy,
														 LPXLOPER XL_precision);

__declspec(dllexport) LPXLOPER WINAPI Local_GetCorrelQuanto(LPXLOPER XL_date1,
															LPXLOPER XL_date2,
															LPXLOPER XL_ccy,
															LPXLOPER XL_index,
															LPXLOPER XL_expiry,
															LPXLOPER XL_tenor,
															LPXLOPER XL_cvname1,
															LPXLOPER XL_cvname2,
															LPXLOPER XL_switchinmonth,
															LPXLOPER XL_domccy,
															LPXLOPER XL_domindex,
															LPXLOPER XL_forccy,
															LPXLOPER XL_forindex,
															LPXLOPER XL_type,
															LPXLOPER XL_lambda,
															LPXLOPER XL_precision,
															LPXLOPER XL_calccy,
															LPXLOPER XL_fwdOrNot);

__declspec(dllexport) LPXLOPER WINAPI Local_GetMoyCorrelQuanto(LPXLOPER XL_date1,
															   LPXLOPER XL_date2,
															   LPXLOPER XL_ccy,
															   LPXLOPER XL_index,
															   LPXLOPER XL_expiry,
															   LPXLOPER XL_tenor,
															   LPXLOPER XL_cvname1,
															   LPXLOPER XL_cvname2,
															   LPXLOPER XL_switchinmonth,
															   LPXLOPER XL_domccy,
															   LPXLOPER XL_domindex,
															   LPXLOPER XL_forccy,
															   LPXLOPER XL_forindex,
															   LPXLOPER XL_type,
															   LPXLOPER XL_lambda,
															   LPXLOPER XL_precision,
															   LPXLOPER XL_calccy,
															   LPXLOPER XL_fwdOrNot);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetDealsFromSummitFilter (LPXLOPER XL_filter);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetAsOfVolOrRate(LPXLOPER XL_AsOfDate,
																 LPXLOPER XL_ccy,
																 LPXLOPER XL_index,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_expiry,
																 LPXLOPER XL_matu,
																 LPXLOPER XL_yieldOrVol,
																 LPXLOPER XL_calcMod,
																 LPXLOPER XL_volType);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFwdRatesMatrix(LPXLOPER XL_AsOfDate,
																  LPXLOPER XL_ccy,
																  LPXLOPER XL_index,
																  LPXLOPER XL_cvName);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInfo(LPXLOPER XL_secId,
														LPXLOPER XL_type);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetModelFactorFromSummit(LPXLOPER XL_date,
																		 LPXLOPER XL_model,
																		 LPXLOPER XL_type,
																		 LPXLOPER XL_factorName,
																		 LPXLOPER XL_ccy,
																		 LPXLOPER XL_index,
																		 LPXLOPER XL_cvName,
																		 LPXLOPER XL_calcMethod);

#endif