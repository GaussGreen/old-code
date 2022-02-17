#ifndef ARM_XL_BOND_LOCAL_H
#define ARM_XL_BOND_LOCAL_H




__declspec(dllexport) LPXLOPER WINAPI Local_BOND (LPXLOPER XL_issueDate,
												  LPXLOPER XL_maturityDate,
												  LPXLOPER XL_firstCouponDate,
												  LPXLOPER XL_couponRate,
												  LPXLOPER XL_redemptionPrice,
												  LPXLOPER XL_periodicity,
												  LPXLOPER XL_dayCount,
												  LPXLOPER XL_settleGap,
												  LPXLOPER XL_couponDateFlag,
												  LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BOND (LPXLOPER XL_issueDate,
													  LPXLOPER XL_maturityDate,
													  LPXLOPER XL_firstCouponDate,
													  LPXLOPER XL_couponRate,
													  LPXLOPER XL_redemptionPrice,
													  LPXLOPER XL_periodicity,
													  LPXLOPER XL_dayCount,
													  LPXLOPER XL_settleGap,
													  LPXLOPER XL_couponDateFlag,
													  LPXLOPER XL_ccy);


__declspec(dllexport) LPXLOPER WINAPI Local_RISKYBOND  (LPXLOPER XL_issueDate,
														LPXLOPER XL_maturityDate,
														LPXLOPER XL_firstCouponDate,
														LPXLOPER XL_couponRate,
														LPXLOPER XL_redemptionPrice,
														LPXLOPER XL_periodicity,
														LPXLOPER XL_dayCount,
														LPXLOPER XL_settleGap,
														LPXLOPER XL_couponDateFlag,
														LPXLOPER XL_ccy,
														LPXLOPER XL_sRepo,
														LPXLOPER XL_ssl,
														LPXLOPER XL_recoveryRate);

__declspec(dllexport) LPXLOPER WINAPI Local_RISKYBONDWITHCF (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_redemptionPrice,
															 LPXLOPER XL_periodicity,
															 LPXLOPER XL_dayCount,
															 LPXLOPER XL_settleGap,
															 LPXLOPER XL_couponDateFlag,
															 LPXLOPER XL_ccyId,
															 LPXLOPER XL_yearTerms,
															 LPXLOPER XL_cashFlows,
															 LPXLOPER XL_sRepo,
															 LPXLOPER XL_ssl,
															 LPXLOPER XL_recoveryRate
															 );

__declspec(dllexport) LPXLOPER WINAPI Local_BONDTEC (LPXLOPER XL_issueDate,
													 LPXLOPER XL_maturityDate,
													 LPXLOPER XL_firstCouponDate,
													 LPXLOPER XL_couponRate,
													 LPXLOPER XL_redemptionPrice,
													 LPXLOPER XL_periodicity,
													 LPXLOPER XL_dayCount,
													 LPXLOPER XL_settleGap,
													 LPXLOPER XL_couponDateFlag,
													 LPXLOPER XL_ccy,
													 LPXLOPER XL_tec,
													 LPXLOPER XL_pfId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BONDTEC (LPXLOPER XL_issueDate,
														 LPXLOPER XL_maturityDate,
														 LPXLOPER XL_firstCouponDate,
														 LPXLOPER XL_couponRate,
														 LPXLOPER XL_redemptionPrice,
														 LPXLOPER XL_periodicity,
														 LPXLOPER XL_dayCount,
														 LPXLOPER XL_settleGap,
														 LPXLOPER XL_couponDateFlag,
														 LPXLOPER XL_ccy,
														 LPXLOPER XL_tec,
														 LPXLOPER XL_pfId );

__declspec(dllexport) LPXLOPER WINAPI Local_YTOPRICE (LPXLOPER XL_bondId,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_yield);

__declspec(dllexport) LPXLOPER WINAPI Local_PTOYIELD (LPXLOPER XL_bondId,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_price);

__declspec(dllexport) LPXLOPER WINAPI Local_ZEROCOUPON (LPXLOPER XL_maturityDate,
														LPXLOPER XL_dayCount,
														LPXLOPER XL_issueDate,
														LPXLOPER XL_settleGap);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZEROCOUPON (LPXLOPER XL_maturityDate,
															LPXLOPER XL_dayCount,
															LPXLOPER XL_issueDate,
															LPXLOPER XL_settleGap);

__declspec(dllexport) LPXLOPER WINAPI Local_BDFAPRICE (LPXLOPER XL_bondId,
													   LPXLOPER XL_settlement,
													   LPXLOPER XL_actuPrice,
													   LPXLOPER XL_forwardDate,
													   LPXLOPER XL_repoRate);

__declspec(dllexport) LPXLOPER WINAPI Local_BDREPORATE (LPXLOPER XL_bondId,
														LPXLOPER XL_settlement, 
														LPXLOPER XL_actuPrice,
														LPXLOPER XL_forwardDate,
														LPXLOPER XL_forwardPrice);

__declspec(dllexport) LPXLOPER WINAPI Local_YTODURATION (LPXLOPER XL_bondId,
														 LPXLOPER XL_settlement,
														 LPXLOPER XL_actuRate,
														 LPXLOPER XL_flagCpn);

__declspec(dllexport) LPXLOPER WINAPI Local_YTOCONVEXITY (LPXLOPER XL_bondId,
														  LPXLOPER XL_settlement,
														  LPXLOPER XL_actuRate);

#endif	/* ARM_XL_CCY_LOCAL_H */
