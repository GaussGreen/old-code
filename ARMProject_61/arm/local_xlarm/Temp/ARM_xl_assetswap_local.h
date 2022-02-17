#ifndef ARM_XL_ASSETSWAP_LOCAL_H
#define ARM_XL_ASSETSWAP_LOCAL_H


__declspec(dllexport) LPXLOPER WINAPI Local_BasisSwap (LPXLOPER XL_asOfDate,
												 LPXLOPER XL_delivery,
												 LPXLOPER XL_maturity,
												 LPXLOPER XL_spread1,
												 LPXLOPER XL_ccy1,
												 LPXLOPER XL_index1,
												 LPXLOPER XL_forwardCurve1,
												 LPXLOPER XL_discountCurve1,
												 LPXLOPER XL_ccy2,
												 LPXLOPER XL_index2,
												 LPXLOPER XL_forwardCurve2,
												 LPXLOPER XL_discountCurve2,
												 LPXLOPER XL_outMode,
												 LPXLOPER XL_solve,
												 LPXLOPER XL_amortizationId);

__declspec(dllexport) LPXLOPER WINAPI Local_ASWMargin (LPXLOPER XL_bondMaturity,
												 LPXLOPER XL_bondCoupon,
												 LPXLOPER XL_bondFrequency,
												 LPXLOPER XL_bondBase,
												 LPXLOPER XL_bondPrice,
												 LPXLOPER XL_bondRedemptionPrice,
												 LPXLOPER XL_asOfDate,
												 LPXLOPER XL_delivery,
												 LPXLOPER XL_fixDecompFrequency,
												 LPXLOPER XL_ccy1,
												 LPXLOPER XL_index1,
												 LPXLOPER XL_forwardCurve1,
											     LPXLOPER XL_discountCurve1,
												 LPXLOPER XL_ccy2,
											     LPXLOPER XL_index2,
												 LPXLOPER XL_forwardCurve2,
												 LPXLOPER XL_discountCurve2,
												 LPXLOPER XL_solve,
												 LPXLOPER XL_amortizationId);

__declspec(dllexport) LPXLOPER WINAPI Local_ASWPrice (LPXLOPER XL_bondMaturity,
												LPXLOPER XL_bondCoupon,
												LPXLOPER XL_bondFrequency,
												LPXLOPER XL_bondBase,
												LPXLOPER XL_bondMargin,
												LPXLOPER XL_bondRedemptionPrice,
												LPXLOPER XL_asOfDate,
												LPXLOPER XL_delivery,
												LPXLOPER XL_fixDecompFrequency,
												LPXLOPER XL_ccy1,
												LPXLOPER XL_index1,
												LPXLOPER XL_forwardCurve1,
											    LPXLOPER XL_discountCurve1,
												LPXLOPER XL_ccy2,
											    LPXLOPER XL_index2,
												LPXLOPER XL_forwardCurve2,
												LPXLOPER XL_discountCurve2,
												LPXLOPER XL_solve,
												LPXLOPER XL_amortizationId,);


__declspec(dllexport) LPXLOPER WINAPI Local_FRNMargin (LPXLOPER XL_asOfDate,
												 LPXLOPER XL_delivery,
												 LPXLOPER XL_maturity,
												 LPXLOPER XL_ccy1,
												 LPXLOPER XL_index1,
												 LPXLOPER XL_forwardCurve1,
												 LPXLOPER XL_discountCurve1,
												 LPXLOPER XL_facialMargin,
												 LPXLOPER XL_price,
												 LPXLOPER XL_ccy2,
												 LPXLOPER XL_index2,
												 LPXLOPER XL_forwardCurve2,
												 LPXLOPER XL_discountCurve2,
												 LPXLOPER XL_fixing,
												 LPXLOPER XL_spread,
												 LPXLOPER XL_outMode,
												 LPXLOPER XL_solve,
												 LPXLOPER XL_amortizationId,);

__declspec(dllexport) LPXLOPER WINAPI Local_FRNPrice (LPXLOPER XL_asOfDate,
												LPXLOPER XL_delivery,
												LPXLOPER XL_maturity,
												LPXLOPER XL_ccy1,
												LPXLOPER XL_index1,
												LPXLOPER XL_forwardCurve1,
												LPXLOPER XL_discountCurve1,
												LPXLOPER XL_facialMargin,
												LPXLOPER XL_valoMargin,
												LPXLOPER XL_ccy2,
												LPXLOPER XL_index2,
												LPXLOPER XL_forwardCurve2,
												LPXLOPER XL_discountCurve2,
												LPXLOPER XL_fixing,
												LPXLOPER XL_spread,
												LPXLOPER XL_outMode,
												LPXLOPER XL_solve,
												LPXLOPER XL_amortizationId,);

__declspec(dllexport) LPXLOPER WINAPI Local_CptBPV (LPXLOPER XL_asOfDate,
											  LPXLOPER XL_delivery,
											  LPXLOPER XL_maturity,
											  LPXLOPER XL_zc,
											  LPXLOPER XL_frequency,
											  LPXLOPER XL_dayCount,
											  LPXLOPER XL_ccy,
											  LPXLOPER XL_amortizationId);

__declspec(dllexport) LPXLOPER WINAPI Local_BondASWMargin (LPXLOPER XL_bond,
														   LPXLOPER XL_bondPrice,
														   LPXLOPER XL_asOfDate,
														   LPXLOPER XL_delivery,
														   LPXLOPER XL_ccy1,
														   LPXLOPER XL_ccy2,
														   LPXLOPER XL_index2,
														   LPXLOPER XL_forwardCurve1,
														   LPXLOPER XL_discountCurve1,
														   LPXLOPER XL_forwardCurve2,
														   LPXLOPER XL_discountCurve2,
														   LPXLOPER XL_solve,
														   LPXLOPER XL_amortizationId);

__declspec(dllexport) LPXLOPER WINAPI Local_BondASWPrice (LPXLOPER XL_bond,
														  LPXLOPER XL_bondMargin,
														  LPXLOPER XL_asOfDate,
														  LPXLOPER XL_delivery,
														  LPXLOPER XL_ccy1,
														  LPXLOPER XL_ccy2,
														  LPXLOPER XL_index2,
														  LPXLOPER XL_forwardCurve1,
														  LPXLOPER XL_discountCurve1,
														  LPXLOPER XL_forwardCurve2,
														  LPXLOPER XL_discountCurve2,
														  LPXLOPER XL_solve,
														  LPXLOPER XL_amortizationId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_NextCpnDate (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_maturity,
															 LPXLOPER XL_frequency,
															 LPXLOPER XL_rule,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_intrule);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_PrevCpnDate (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_maturity,
															 LPXLOPER XL_frequency,
															 LPXLOPER XL_rule,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_intrule)

#endif /* ARM_XL_ASSETSWAP_LOCAL_H */