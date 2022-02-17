#ifndef ARM_XL_ZCCURVE_LOCAL_H
#define ARM_XL_ZCCURVE_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapInt (LPXLOPER XL_date,
															 LPXLOPER XL_matuRate,
															 LPXLOPER XL_mmVsFut,
															 LPXLOPER XL_swapVsFut,
															 LPXLOPER XL_raw,
															 LPXLOPER XL_interp,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_swapFrq,
															 LPXLOPER XL_fixDayCount);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapInt (LPXLOPER XL_date,
																 LPXLOPER XL_matuRate,
																 LPXLOPER XL_mmVsFut,
																 LPXLOPER XL_swapVsFut,
																 LPXLOPER XL_raw,
																 LPXLOPER XL_interp,
																 LPXLOPER XL_ccy,
																 LPXLOPER XL_swapFrq,
															     LPXLOPER XL_fixDayCount);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCLINT (LPXLOPER XL_matu,
													LPXLOPER XL_rate,
													LPXLOPER XL_meth,
													LPXLOPER XL_aDate,
													LPXLOPER XL_ccy,
													LPXLOPER XL_interpMeth);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCLINT (LPXLOPER XL_matu,
														LPXLOPER XL_rate,
														LPXLOPER XL_meth,
														LPXLOPER XL_aDate,
														LPXLOPER XL_ccy,
														LPXLOPER XL_interpMeth);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCFLAT (LPXLOPER XL_zeroFlat,
													LPXLOPER XL_date,
                                                    LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCFLAT (LPXLOPER XL_zeroFlat,
														LPXLOPER XL_date,
                                                        LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCSPREADED (LPXLOPER XL_zcSprId,
														LPXLOPER XL_zcInitId,
														LPXLOPER XL_date,
														LPXLOPER XL_mmFreq,
														LPXLOPER XL_swapFreq,
														LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSPREADED (LPXLOPER XL_zcSprId,
															LPXLOPER XL_zcInitId,
															LPXLOPER XL_date,
															LPXLOPER XL_mmFreq,
															LPXLOPER XL_swapFreq,
															LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSpreadedFromSummit(LPXLOPER XL_index,
																	   LPXLOPER XL_currency,
																	   LPXLOPER XL_cvName,
																	   LPXLOPER XL_aSdate,
																	   LPXLOPER XL_convAdj,
																	   LPXLOPER XL_raw,
																	   LPXLOPER XL_mmFreq,
																	   LPXLOPER XL_swapFrq,
																	   LPXLOPER XL_interp);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSpreadedFromSummit(LPXLOPER XL_index,
																		   LPXLOPER XL_currency,
																		   LPXLOPER XL_cvName,
																		   LPXLOPER XL_aSdate,
																		   LPXLOPER XL_convAdj,
																		   LPXLOPER XL_raw,
																		   LPXLOPER XL_mmFreq,
																		   LPXLOPER XL_swapFrq,
																		   LPXLOPER XL_interp);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapIntSmooth (LPXLOPER XL_date,
																   LPXLOPER XL_matuRate,
																   LPXLOPER XL_mmVsFut,
																   LPXLOPER XL_swapVsFut,
																   LPXLOPER XL_raw,
																   LPXLOPER XL_interp,
																   LPXLOPER XL_ccy,
																   LPXLOPER XL_lambda,
																   LPXLOPER XL_prec);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapIntSmooth (LPXLOPER XL_date,
																	   LPXLOPER XL_matuRate,
																	   LPXLOPER XL_mmVsFut,
																	   LPXLOPER XL_swapVsFut,
																	   LPXLOPER XL_raw,
																	   LPXLOPER XL_interp,
																	   LPXLOPER XL_ccy,
																	   LPXLOPER XL_lambda,
																	   LPXLOPER XL_prec);

__declspec(dllexport) LPXLOPER WINAPI Local_DiscountPrice (LPXLOPER XL_curve,
														   LPXLOPER XL_matu);
__declspec(dllexport) LPXLOPER WINAPI Local_DiscountYield (LPXLOPER XL_curve,
														   LPXLOPER XL_matu,
														   LPXLOPER XL_meth);

__declspec(dllexport) LPXLOPER WINAPI Local_ForwardPrice (LPXLOPER XL_curve,
														  LPXLOPER XL_matu1,
														  LPXLOPER XL_matu2);
__declspec(dllexport) LPXLOPER WINAPI Local_ForwardYield (LPXLOPER XL_curve,
														  LPXLOPER XL_matu1,
														  LPXLOPER XL_matu2,
														  LPXLOPER XL_meth,
														  LPXLOPER XL_adjDaycount,
														  LPXLOPER XL_daycount);

__declspec(dllexport) LPXLOPER WINAPI Local_GenForwardYield (LPXLOPER XL_curve,
														  LPXLOPER XL_matu1,
														  LPXLOPER XL_matu2,
														  LPXLOPER XL_isSwapRate,
														  LPXLOPER XL_decompFreq,
														  LPXLOPER XL_daycount);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapFutInt (LPXLOPER XL_date,
																LPXLOPER XL_matuRate,
																LPXLOPER XL_mmVsFut,
																LPXLOPER XL_swapVsFut,
																LPXLOPER XL_raw,
																LPXLOPER XL_interp,
																LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapFutInt (LPXLOPER XL_date,
																	LPXLOPER XL_matuRate,
																	LPXLOPER XL_mmVsFut,
																	LPXLOPER XL_swapVsFut,
																	LPXLOPER XL_raw,
																	LPXLOPER XL_interp,
																	LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapFutIntSmooth (LPXLOPER XL_date,
																	  LPXLOPER XL_matuRate,
																	  LPXLOPER XL_mmVsFut,
																	  LPXLOPER XL_swapVsFut,
																	  LPXLOPER XL_raw,
																	  LPXLOPER XL_interp,
																	  LPXLOPER XL_ccy,
																	  LPXLOPER XL_lambda,
																	  LPXLOPER XL_prec);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapFutIntSmooth (LPXLOPER XL_date,
																		  LPXLOPER XL_matuRate,
																		  LPXLOPER XL_mmVsFut,
																		  LPXLOPER XL_swapVsFut,
																		  LPXLOPER XL_raw,
																		  LPXLOPER XL_interp,
																		  LPXLOPER XL_ccy,
																		  LPXLOPER XL_lambda,
																		  LPXLOPER XL_prec);

__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedZCSWAPINT (LPXLOPER XL_shiftValue,
															  LPXLOPER XL_nbPlot,
															  LPXLOPER XL_date,
															  LPXLOPER XL_matuRate,
															  LPXLOPER XL_mmVsFut,
															  LPXLOPER XL_swapVsFut,
															  LPXLOPER XL_raw,
															  LPXLOPER XL_interp,
															  LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftedZCSWAPINT (LPXLOPER XL_shiftValue,
																  LPXLOPER XL_nbPlot,
																  LPXLOPER XL_date,
																  LPXLOPER XL_matuRate,
																  LPXLOPER XL_mmVsFut,
																  LPXLOPER XL_swapVsFut,
																  LPXLOPER XL_raw,
																  LPXLOPER XL_interp,
																  LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedZCSWAPINTSmooth (LPXLOPER XL_shiftValue,
																	LPXLOPER XL_nbPlot,
																	LPXLOPER XL_date,
																	LPXLOPER XL_matuRate,
																	LPXLOPER XL_mmVsFut,
																	LPXLOPER XL_swapVsFut,
																	LPXLOPER XL_raw,
																	LPXLOPER XL_interp,
																	LPXLOPER XL_lambda,
																	LPXLOPER XL_prec,
																	LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftedZCSWAPINTSmooth (LPXLOPER XL_shiftValue,
																		LPXLOPER XL_nbPlot,
																		LPXLOPER XL_date,
																		LPXLOPER XL_matuRate,
																		LPXLOPER XL_mmVsFut,
																		LPXLOPER XL_swapVsFut,
																		LPXLOPER XL_raw,
																		LPXLOPER XL_interp,
																		LPXLOPER XL_lambda,
																		LPXLOPER XL_prec,
																		LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_TRANS2SMOOTH (LPXLOPER XL_inCv,
														  LPXLOPER XL_lambda,
														  LPXLOPER XL_prec);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TRANS2SMOOTH (LPXLOPER XL_inCv,
															  LPXLOPER XL_lambda,
															  LPXLOPER XL_prec);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCSWAPCUBDIFF (LPXLOPER XL_date,
														   LPXLOPER XL_matuRate,
														   LPXLOPER XL_mmVsFut,
														   LPXLOPER XL_swapVsFut,
														   LPXLOPER XL_raw,
														   LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSWAPCUBDIFF (LPXLOPER XL_date,
															   LPXLOPER XL_matuRate,
															   LPXLOPER XL_mmVsFut,
															   LPXLOPER XL_swapVsFut,
															   LPXLOPER XL_raw,
															   LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCINTSMOOTH (LPXLOPER XL_matu,
														 LPXLOPER XL_rate,
														 LPXLOPER XL_aDate,
														 LPXLOPER XL_meth,
														 LPXLOPER XL_lambda,
														 LPXLOPER XL_prec);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCINTSMOOTH (LPXLOPER XL_matu,
															 LPXLOPER XL_rate,
															 LPXLOPER XL_aDate,
															 LPXLOPER XL_meth,
															 LPXLOPER XL_lambda,
															 LPXLOPER XL_prec);

__declspec(dllexport) LPXLOPER WINAPI Local_GetZCFromSummit (LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_cvName,
															 LPXLOPER XL_aSdate);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetZCFromSummit (LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_aSdate);

__declspec(dllexport) LPXLOPER WINAPI Local_GetMaturitiesFromZC(LPXLOPER XL_cvName);

__declspec(dllexport) LPXLOPER WINAPI Local_GetInitialCurveFromSummit (LPXLOPER XL_index,
																	   LPXLOPER XL_currency,
																	   LPXLOPER XL_cvName,
																	   LPXLOPER XL_aSdate,
																	   LPXLOPER XL_adjOrNot);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CreateZCFromSummit (LPXLOPER XL_index,
																	LPXLOPER XL_currency,
																	LPXLOPER XL_cvName,
																	LPXLOPER XL_aSdate,
																	LPXLOPER XL_convAdj,
																	LPXLOPER XL_fwd,
																	LPXLOPER XL_swapFrq);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CreateZCFromSummit (LPXLOPER XL_index,
																		LPXLOPER XL_currency,
																		LPXLOPER XL_cvName,
																		LPXLOPER XL_aSdate,
																		LPXLOPER XL_convAdj,
																		LPXLOPER XL_fwd,
																		LPXLOPER XL_swapFrq);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCVSK (LPXLOPER XL_param,
												   LPXLOPER XL_date);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCVSK (LPXLOPER XL_param,
													   LPXLOPER XL_date);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCSPLICUB (LPXLOPER XL_matu,
													   LPXLOPER XL_rate,
													   LPXLOPER XL_meth,
													   LPXLOPER XL_date,
													   LPXLOPER XL_lastBucket);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSPLICUB (LPXLOPER XL_matu,
														   LPXLOPER XL_rate,
														   LPXLOPER XL_meth,
														   LPXLOPER XL_date,
														   LPXLOPER XL_lastBucket);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCSPL (LPXLOPER XL_param,
												   LPXLOPER XL_date);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSPL (LPXLOPER XL_param,
													   LPXLOPER XL_date);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpCurve (LPXLOPER XL_inCv,
														   LPXLOPER XL_epsilon);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpCurve (LPXLOPER XL_inCv,
															   LPXLOPER XL_epsilon);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpSpreadedCurve(	LPXLOPER XL_inCv,
																	LPXLOPER XL_epsilon,
																	LPXLOPER XL_curveToBump,
																	LPXLOPER XL_meth);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpSpreadedCurve(	LPXLOPER XL_inCv,
																		LPXLOPER XL_epsilon,
																		LPXLOPER XL_curveToBump,
																		LPXLOPER XL_meth);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateZeroCurveLin (LPXLOPER XL_matuRate,
																LPXLOPER XL_meth,
																LPXLOPER XL_date,
																LPXLOPER XL_ccy,
                                                                LPXLOPER XL_interpMeth);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZeroCurveLin (LPXLOPER XL_matuRate,
																	LPXLOPER XL_meth,
																	LPXLOPER XL_date,
																	LPXLOPER XL_ccy,
                                                                    LPXLOPER XL_interpMeth);

__declspec(dllexport) LPXLOPER WINAPI Local_ShiftZeroCurveLin (LPXLOPER XL_value,
                                                               LPXLOPER XL_nbPlot,
                                                               LPXLOPER XL_matuRate,
																LPXLOPER XL_meth,
																LPXLOPER XL_date,
																LPXLOPER XL_ccy,
																LPXLOPER XL_interpMeth)

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftZeroCurveLin (LPXLOPER XL_value,
                                                                    LPXLOPER XL_nbPlot,
                                                                    LPXLOPER XL_matuRate,
																	LPXLOPER XL_meth,
																	LPXLOPER XL_date,
																	LPXLOPER XL_ccy,
																	LPXLOPER XL_interpMeth);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCCashInt (LPXLOPER XL_date,
															 LPXLOPER XL_matuRate,
															 LPXLOPER XL_bonds,
															 LPXLOPER XL_mmVsFut,
															 LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCCashInt (LPXLOPER XL_date,
																 LPXLOPER XL_matuRate,
																 LPXLOPER XL_bonds,
																 LPXLOPER XL_mmVsFut,
																 LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_GenerateBasisAdjCurve(LPXLOPER XL_DomCrv1,
																  LPXLOPER XL_DomBSCrv1,
																  LPXLOPER XL_ForCrv2,
																  LPXLOPER XL_ForBSCrv2,
																  LPXLOPER XL_BSasSprds,
																  LPXLOPER XL_RetSprds,
																  LPXLOPER XL_matuVec);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenerateBasisAdjCurve(LPXLOPER XL_Crv1,
																	  LPXLOPER XL_BSCrv1,
																	  LPXLOPER XL_Crv2,
																	  LPXLOPER XL_BSCrv2,
																	  LPXLOPER XL_RetSprds,
																	  LPXLOPER XL_matuVec,
																	  LPXLOPER XL_BSasSprds);

__declspec(dllexport) LPXLOPER WINAPI Local_ZCSWAPSPLSUM (LPXLOPER XL_date,
														   LPXLOPER XL_matuRate,
														   LPXLOPER XL_mmVsFut,
														   LPXLOPER XL_swapVsFut,
														   LPXLOPER XL_raw,
														   LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSWAPSPLSUM (LPXLOPER XL_date,
															   LPXLOPER XL_matuRate,
															   LPXLOPER XL_mmVsFut,
															   LPXLOPER XL_swapVsFut,
															   LPXLOPER XL_raw,
															   LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FixingSched(LPXLOPER XL_asOfDate,
															LPXLOPER XL_LiborFixing,
															LPXLOPER XL_FXFixing);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFixingSchedFromSummit(LPXLOPER XL_ListOfKeys,
																		 LPXLOPER XL_AsOf,
																		 LPXLOPER XL_Source,
																		 LPXLOPER XL_dateStrip);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FixingSched(LPXLOPER XL_asOfDate,
																LPXLOPER XL_LiborFixing,
																LPXLOPER XL_FXFixing);

#endif	/* ARM_XL_ZCCURVE_LOCAL_H */
