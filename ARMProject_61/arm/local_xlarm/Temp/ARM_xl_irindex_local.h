#ifndef ARM_XL_IRINDEX_LOCAL_H
#define ARM_XL_IRINDEX_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_LIBOR (LPXLOPER XL_liborType,
												   LPXLOPER XL_ccy,
												   LPXLOPER XL_resetFreq,
												   LPXLOPER XL_payFreq,
                                                   LPXLOPER XL_daycount,
												   LPXLOPER XL_intRule);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBOR (LPXLOPER XL_liborType,
													   LPXLOPER XL_ccy,
													   LPXLOPER XL_resetFreq,
													   LPXLOPER XL_payFreq,
                                                       LPXLOPER XL_daycount,
													   LPXLOPER XL_intRule);

__declspec(dllexport) LPXLOPER WINAPI Local_IRINDEX (LPXLOPER XL_dayCount,
													 LPXLOPER XL_payfreq,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_compMethod,
													 LPXLOPER XL_fwdRule,
													 LPXLOPER XL_resetTiming,
													 LPXLOPER XL_resetGap,
													 LPXLOPER XL_payTiming,
													 LPXLOPER XL_payGap,
													 LPXLOPER XL_ccyId,
													 LPXLOPER XL_indexType,
													 LPXLOPER XL_decompFreq,
													 LPXLOPER XL_intRule,
													 LPXLOPER XL_resetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRINDEX (LPXLOPER XL_dayCount,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_maturity,
														 LPXLOPER XL_compMethod,
														 LPXLOPER XL_fwdRule,
														 LPXLOPER XL_resetTiming,
														 LPXLOPER XL_resetGap,
														 LPXLOPER XL_payTiming,
														 LPXLOPER XL_payGap,
														 LPXLOPER XL_ccy,
														 LPXLOPER XL_indexType,
														 LPXLOPER XL_decompFreq,
														 LPXLOPER XL_intRule,
														 LPXLOPER XL_resetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_MultiIrindex (LPXLOPER XL_irIndexVec,
														  LPXLOPER XL_weightVec);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MultiIrindex (LPXLOPER XL_irIndexVec,
														  LPXLOPER XL_weightVec);


__declspec(dllexport) LPXLOPER WINAPI Local_FixedIndex(LPXLOPER XL_dayCount,
													   LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FixedIndex(LPXLOPER XL_dayCount,
														   LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_IRINDEX_MONEY_MARKET (LPXLOPER XL_mmTerm,
																  LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRINDEX_MONEY_MARKET (LPXLOPER XL_mmTerm,
																	  LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_CMS(LPXLOPER XL_CMSType,
												LPXLOPER XL_liborType,
												LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CMS(LPXLOPER XL_CMSType,
													LPXLOPER XL_liborType,
													LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_IRINDEX2 (LPXLOPER XL_payFreq,
													  LPXLOPER XL_maturity,
													  LPXLOPER XL_ccy,
													  LPXLOPER XL_indexType,
													  LPXLOPER XL_resetFreq);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRINDEX2 (LPXLOPER XL_payFreq,
													      LPXLOPER XL_maturity,
													      LPXLOPER XL_ccy,
													      LPXLOPER XL_indexType,
														  LPXLOPER XL_resetFreq);

#endif /* ARM_XL_IRINDEX_LOCAL_H */