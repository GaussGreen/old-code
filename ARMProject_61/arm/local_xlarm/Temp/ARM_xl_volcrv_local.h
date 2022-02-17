#ifndef ARM_XL_VOLCRV_LOCAL_H
#define ARM_XL_VOLCRV_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_VolCurv(LPXLOPER XL_matu,
													    LPXLOPER XL_strikes,
													    LPXLOPER XL_vols,
													    LPXLOPER XL_date,
													    LPXLOPER XL_volType,
													    LPXLOPER XL_strikeType,
                                                        LPXLOPER XL_InterpolType,
                                                        LPXLOPER XL_ccy,
														LPXLOPER XL_indexName);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FxVolCurv(LPXLOPER XL_date,
                                                          LPXLOPER XL_matus,
                                                          LPXLOPER XL_fxFwds,
													      LPXLOPER XL_pivotVols,
                                                          LPXLOPER XL_pivotTypes,
													      LPXLOPER XL_deltaPuts,
													      LPXLOPER XL_deltaCalls,
                                                          LPXLOPER XL_volsPuts,
                                                          LPXLOPER XL_volsCalls,
                                                          LPXLOPER XL_InterpolTypes,
													      LPXLOPER XL_WhatIsInterpolated,
													      LPXLOPER XL_correctSplineWithLinear,
                                                          LPXLOPER XL_fxSpot,
                                                          LPXLOPER XL_domZcCurve,
                                                          LPXLOPER XL_forZcCurve,
                                                          LPXLOPER XL_inRRSTR,
                                                          LPXLOPER XL_isATM);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FxVolCurv(LPXLOPER XL_date,
	                                                          LPXLOPER XL_matus,
	                                                          LPXLOPER XL_fxFwds,
														      LPXLOPER XL_pivotVols,
	                                                          LPXLOPER XL_pivotTypes,
														      LPXLOPER XL_deltaPuts,
														      LPXLOPER XL_deltaCalls,
	                                                          LPXLOPER XL_volsPuts,
	                                                          LPXLOPER XL_volsCalls,
	                                                          LPXLOPER XL_InterpolTypes,
														      LPXLOPER XL_WhatIsInterpolated,
														      LPXLOPER XL_correctSplineWithLinear,
	                                                          LPXLOPER XL_fxSpot,
	                                                          LPXLOPER XL_domZcCurve,
	                                                          LPXLOPER XL_forZcCurve,
	                                                          LPXLOPER XL_inRRSTR,
	                                                          LPXLOPER XL_isATM);

__declspec(dllexport) LPXLOPER WINAPI Local_NewComputeFxVolatility(LPXLOPER XL_curve,
															       LPXLOPER XL_matu,
															       LPXLOPER XL_moneyness);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolCurv(LPXLOPER XL_matu,
														LPXLOPER XL_strikes,
														LPXLOPER XL_vols,
														LPXLOPER XL_date,
														LPXLOPER XL_volType,
														LPXLOPER XL_strikeType,
                                                        LPXLOPER XL_InterpolType,
                                                        LPXLOPER XL_ccy,
														LPXLOPER XL_indexName);

__declspec(dllexport) LPXLOPER WINAPI Local_VolFlat(LPXLOPER XL_vol,
													LPXLOPER XL_date,
                                                    LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolFlat(LPXLOPER XL_vol,
														LPXLOPER XL_date,
                                                        LPXLOPER XL_ccy);



__declspec(dllexport) LPXLOPER WINAPI Local_ComputeVolatility (LPXLOPER XL_curve,
															   LPXLOPER XL_matu,
															   LPXLOPER XL_strike,
															   LPXLOPER XL_tenor);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeModelVolatility (LPXLOPER XL_model,
																	LPXLOPER XL_matu,
																	LPXLOPER XL_tenor,
																	LPXLOPER XL_fwd,
																	LPXLOPER XL_strike);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrel (LPXLOPER XL_Correl,
															   LPXLOPER XL_X,
															   LPXLOPER XL_Y );

__declspec(dllexport) LPXLOPER WINAPI Local_VolCube(LPXLOPER XL_ATMVolId,
													LPXLOPER XL_volCurves,
													LPXLOPER XL_tenors,
                                                    LPXLOPER XL_volType,
													LPXLOPER XL_CheckCcy);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolCube(LPXLOPER XL_ATMVolId,
														LPXLOPER XL_volCurves,
														LPXLOPER XL_tenors,
                                                        LPXLOPER XL_volType,
														LPXLOPER XL_CheckCcy);

__declspec(dllexport) LPXLOPER WINAPI Local_GetVolFromCalypso(LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_cvName,
															 LPXLOPER XL_date,
															 LPXLOPER XL_vtype,
															 LPXLOPER XL_matIndex,
															 LPXLOPER XL_impOrHist,
															 LPXLOPER XL_indexid);

__declspec(dllexport) LPXLOPER WINAPI Local_GetVolFromSummit(LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_cvName,
															 LPXLOPER XL_date,
															 LPXLOPER XL_vtype,
															 LPXLOPER XL_matIndex,
															 LPXLOPER XL_impOrHist,
															 LPXLOPER XL_indexid);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetVolFromSummit(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
																 LPXLOPER XL_matIndex,
																 LPXLOPER XL_impOrHist,
																 LPXLOPER XL_indexid);

__declspec(dllexport) LPXLOPER WINAPI Local_GetVolCubeFromCalypso(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
                                                                 LPXLOPER XL_suffix,
                                                                 LPXLOPER XL_type,
																 LPXLOPER XL_tenors,
																 LPXLOPER XL_SmileOrNot
																 LPXLOPER XL_indexid);
__declspec(dllexport) LPXLOPER WINAPI Local_GetVolCubeFromSummit(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
																 LPXLOPER XL_tenors,
																 LPXLOPER XL_SmileOrNot
																 LPXLOPER XL_indexid,
																 LPXLOPER XL_SmileType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetVolCubeFromSummit(LPXLOPER XL_index,
																	LPXLOPER XL_currency,
																	LPXLOPER XL_cvName,
																	LPXLOPER XL_date,
																	LPXLOPER XL_vtype,
																	LPXLOPER XL_tenors,
																	LPXLOPER XL_SmileOrNot,
																	LPXLOPER XL_indexid,
																	LPXLOPER XL_SmileType);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpVolatility (LPXLOPER XL_vol,
																LPXLOPER XL_value,
																LPXLOPER XL_nthLine,
																LPXLOPER XL_nthCol,
																LPXLOPER XL_isCumul,
																LPXLOPER XL_isAbsolute);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpVolatility (LPXLOPER XL_vol,
																	LPXLOPER XL_value,
																	LPXLOPER XL_nthLine,
																	LPXLOPER XL_nthCol,
																	LPXLOPER XL_isCumul,
																	LPXLOPER XL_isAbsolute);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpSmile (LPXLOPER XL_vol,
														   LPXLOPER XL_value,
														   LPXLOPER XL_tenor,
														   LPXLOPER XL_nthLine,
														   LPXLOPER XL_nthCol,
														   LPXLOPER XL_isCumul,
														   LPXLOPER XL_isAbsolute);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpSmile (LPXLOPER XL_vol,
															   LPXLOPER XL_value,
															   LPXLOPER XL_tenor,
															   LPXLOPER XL_nthLine,
															   LPXLOPER XL_nthCol,
															   LPXLOPER XL_isCumul,
															   LPXLOPER XL_isAbsolute);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpHyperCubeSmile(LPXLOPER XL_vol,
																  LPXLOPER XL_value,
																  LPXLOPER XL_cubeTenor,
																  LPXLOPER XL_smileTenor,
																  LPXLOPER XL_nthLine,
																  LPXLOPER XL_nthCol,
																  LPXLOPER XL_isCumul,
																  LPXLOPER XL_isAbsolute);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpHyperCubeSmile(LPXLOPER XL_vol,
																	  LPXLOPER XL_value,
																	  LPXLOPER XL_cubeTenor,
																	  LPXLOPER XL_smileTenor,
																	  LPXLOPER XL_nthLine,
																	  LPXLOPER XL_nthCol,
																	  LPXLOPER XL_isCumul,
																	  LPXLOPER XL_isAbsolute);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FXBumpRRorSTR ( LPXLOPER XL_vol,
																LPXLOPER XL_value,
																LPXLOPER XL_nthLine,
																LPXLOPER XL_nthCol,
																LPXLOPER XL_spotFX,
																LPXLOPER XL_isCumul,
																LPXLOPER XL_isAbsolute,
																LPXLOPER XL_isRR);

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FXBumpRRorSTR ( LPXLOPER XL_vol,
																LPXLOPER XL_value,
																LPXLOPER XL_nthLine,
																LPXLOPER XL_nthCol,
																LPXLOPER XL_spotFX,
																LPXLOPER XL_isCumul,
																LPXLOPER XL_isAbsolute,
																LPXLOPER XL_isRR);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFXVolFromSummit(LPXLOPER XL_ccy1,
																   LPXLOPER XL_ccy2,
																   LPXLOPER XL_date,
																   LPXLOPER XL_cvName,
																   LPXLOPER XL_impOrHist,
																   LPXLOPER XL_volType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetFXVolFromSummit(LPXLOPER XL_ccy1,
																	   LPXLOPER XL_ccy2,
																	   LPXLOPER XL_date,
																	   LPXLOPER XL_cvName,
																	   LPXLOPER XL_impOrHist,
																	   LPXLOPER XL_volType);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetNewFXVolFromSummit(LPXLOPER XL_ccy1,
																	  LPXLOPER XL_ccy2,
																	  LPXLOPER XL_date,
																	  LPXLOPER XL_cvName,
																	  LPXLOPER XL_domZc,
																	  LPXLOPER XL_forZc,
																	  LPXLOPER XL_fxSpot,
																	  LPXLOPER XL_forwards,
																	  LPXLOPER XL_WhatIsInterpolated,
																	  LPXLOPER XL_correctSplineWithLinear,
																	  LPXLOPER XL_isATM);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInitialFXVolFromSummit(LPXLOPER XL_ccy1,
																		  LPXLOPER XL_ccy2,
																		  LPXLOPER XL_date,
																		  LPXLOPER XL_cvName,
																		  LPXLOPER XL_impOrHist,
																		  LPXLOPER XL_volType);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInitialVolFromSummit(LPXLOPER XL_index,
																		LPXLOPER XL_currency,
																		LPXLOPER XL_cvName,
																		LPXLOPER XL_date,
																		LPXLOPER XL_vtype,
																		LPXLOPER XL_matIndex);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFXCorrelFromSummit(LPXLOPER XL_ccy1,
																	  LPXLOPER XL_index,
																	  LPXLOPER XL_ccy2,
																	  LPXLOPER XL_asof,
																	  LPXLOPER XL_cvName,
																	  LPXLOPER XL_tenors);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetFXCorrelFromSummit(LPXLOPER XL_ccy1,
																		  LPXLOPER XL_index,
																		  LPXLOPER XL_ccy2,
																		  LPXLOPER XL_asof,
																		  LPXLOPER XL_cvName,
																		  LPXLOPER XL_tenors);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCorrelFromSummit(LPXLOPER XL_ccy1,
																	LPXLOPER XL_index1,
																	LPXLOPER XL_ccy2,
																	LPXLOPER XL_index2,
																	LPXLOPER XL_asof,
																	LPXLOPER XL_cvName);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetCorrelFromSummit(LPXLOPER XL_ccy1,
																		LPXLOPER XL_index1,
																		LPXLOPER XL_ccy2,
																		LPXLOPER XL_index2,
																		LPXLOPER XL_asof,
																		LPXLOPER XL_cvName);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetNthMaturity (LPXLOPER XL_curve,
																LPXLOPER XL_nLine);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONVERTFROMBSTONORMALVOL (LPXLOPER XL_vol,
																		  LPXLOPER XL_zc,
																		  LPXLOPER XL_isSwoptVol,
																		  LPXLOPER XL_inPct);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONVERTFROMBSTONORMALVOL (LPXLOPER XL_vol,
																			  LPXLOPER XL_zc,
																			  LPXLOPER XL_isSwoptVol,
																			  LPXLOPER XL_inPct);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONVERTFROMNORMALTOBSVOL (LPXLOPER XL_vol,
																		  LPXLOPER XL_zc,
																		  LPXLOPER XL_isSwoptVol,
																		  LPXLOPER XL_inPct,
																		  LPXLOPER XL_post100);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONVERTFROMNORMALTOBSVOL (LPXLOPER XL_vol,
																			  LPXLOPER XL_zc,
																			  LPXLOPER XL_isSwoptVol,
																			  LPXLOPER XL_inPct,
																			  LPXLOPER XL_post100);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS(LPXLOPER XL_tree3f,
																	           LPXLOPER XL_PrcsObject,
                                                                               LPXLOPER XL_ForwardVolDates);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS( LPXLOPER XL_tree3f,
																	                LPXLOPER XL_PrcsObject,
                                                                                    LPXLOPER XL_ForwardVolDates);


__declspec(dllexport) LPXLOPER WINAPI Local_InterpolInStrikeFwdTime (LPXLOPER XL_curve,
																	 LPXLOPER XL_forward,
																	 LPXLOPER XL_strike,
																	 LPXLOPER XL_matu,
																	 LPXLOPER XL_precision,
																	 LPXLOPER XL_sigmaATM,
																	 LPXLOPER XL_y2NULL);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeFxVol (LPXLOPER XL_fxvol,
														  LPXLOPER XL_asof,
														  LPXLOPER XL_matu,
														  LPXLOPER XL_calcmatu,
														  LPXLOPER XL_fxspot,
														  LPXLOPER XL_strike,
														  LPXLOPER XL_discCrv,
														  LPXLOPER XL_divCrv);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMSPOTTOFWDVOL( LPXLOPER XL_asOf,
                                                                        LPXLOPER XL_dZcId,
                                                                        LPXLOPER XL_fZcId,
                                                                        LPXLOPER XL_dBSZcId,
                                                                        LPXLOPER XL_fBSZcId,
                                                                        LPXLOPER XL_volSwopBaseId,
                                                                        LPXLOPER XL_volSwopForeignId,
                                                                        LPXLOPER XL_fxVolId,
                                                                        LPXLOPER XL_dMeanReversionBase,
                                                                        LPXLOPER XL_dMeanReversionForeign,
                                                                        LPXLOPER XL_dFxRdCorrId,
                                                                        LPXLOPER XL_dFxRfCorrId,
                                                                        LPXLOPER XL_dRdRfCorrId,
                                                                        LPXLOPER XL_dCutOff,
                                                                        LPXLOPER XL_dVolLongTerm,
                                                                        LPXLOPER XL_calibBasisIncluded,
                                                                        LPXLOPER XL_ForwardVolDates);  

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL( LPXLOPER XL_asOf,
                                                                            LPXLOPER XL_dZcId,
                                                                            LPXLOPER XL_fZcId,
                                                                            LPXLOPER XL_dBSZcId,
                                                                            LPXLOPER XL_fBSZcId,
                                                                            LPXLOPER XL_volSwopBaseId,
                                                                            LPXLOPER XL_volSwopForeignId,
                                                                            LPXLOPER XL_fxVolId,
                                                                            LPXLOPER XL_dMeanReversionBase,
                                                                            LPXLOPER XL_dMeanReversionForeign,
                                                                            LPXLOPER XL_dFxRdCorrId,
                                                                            LPXLOPER XL_dFxRfCorrId,
                                                                            LPXLOPER XL_dRdRfCorrId,
                                                                            LPXLOPER XL_dCutOff,
                                                                            LPXLOPER XL_dVolLongTerm,
                                                                            LPXLOPER XL_calibBasisIncluded,
                                                                            LPXLOPER XL_ForwardVolDates);  


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS( LPXLOPER XL_tree3f,
																	            LPXLOPER XL_PrcsObject);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS (LPXLOPER XL_tree3f,
																	                LPXLOPER XL_PrcsObject);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMFWDTOSPOTVOL( LPXLOPER XL_asOf,
                                                                        LPXLOPER XL_dZcId,
                                                                        LPXLOPER XL_fZcId,
                                                                        LPXLOPER XL_dBSZcId,
                                                                        LPXLOPER XL_fBSZcId,
                                                                        LPXLOPER XL_volSwopBaseId,
                                                                        LPXLOPER XL_volSwopForeignId,
                                                                        LPXLOPER XL_fxVolId,
                                                                        LPXLOPER XL_dMeanReversionBase,
                                                                        LPXLOPER XL_dMeanReversionForeign,
                                                                        LPXLOPER XL_dFxRdCorrId,
                                                                        LPXLOPER XL_dFxRfCorrId,
                                                                        LPXLOPER XL_dRdRfCorrId,
                                                                        LPXLOPER XL_dCutOff,
                                                                        LPXLOPER XL_dVolLongTerm,
                                                                        LPXLOPER XL_calibBasisIncluded);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL( LPXLOPER XL_asOf,
                                                                            LPXLOPER XL_dZcId,
                                                                            LPXLOPER XL_fZcId,
                                                                            LPXLOPER XL_dBSZcId,
                                                                            LPXLOPER XL_fBSZcId,
                                                                            LPXLOPER XL_volSwopBaseId,
                                                                            LPXLOPER XL_volSwopForeignId,
                                                                            LPXLOPER XL_fxVolId,
                                                                            LPXLOPER XL_dMeanReversionBase,
                                                                            LPXLOPER XL_dMeanReversionForeign,
                                                                            LPXLOPER XL_dFxRdCorrId,
                                                                            LPXLOPER XL_dFxRfCorrId,
                                                                            LPXLOPER XL_dRdRfCorrId,
                                                                            LPXLOPER XL_dCutOff,
                                                                            LPXLOPER XL_dVolLongTerm,
                                                                            LPXLOPER XL_calibBasisIncluded);

__declspec(dllexport) LPXLOPER WINAPI Local_SetVolCurveName(LPXLOPER XL_curve,
															LPXLOPER XL_name);

// Hypercube
__declspec(dllexport) LPXLOPER WINAPI Local_HyperCube(LPXLOPER XL_volCubes,
													  LPXLOPER XL_keys);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HyperCube(LPXLOPER XL_volCubes,
														  LPXLOPER XL_keys);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateCorrelCubeByExpiry(LPXLOPER XL_hyperCube,
																	 LPXLOPER XL_TenorList,
																	 LPXLOPER XL_ExpiryList,
																	 LPXLOPER XL_IntersurfaceInterpol);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateCorrelCubeByExpiry(LPXLOPER XL_hyperCube,
																	     LPXLOPER XL_TenorList,
																	     LPXLOPER XL_ExpiryList,
																	     LPXLOPER XL_IntersurfaceInterpol);


__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrelByExpiry(LPXLOPER XL_correlCube,
																  LPXLOPER XL_Expiry,
																  LPXLOPER XL_Tenor1,
																  LPXLOPER XL_Tenor2);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeHyperCorrel(LPXLOPER XL_correlCube,
														       LPXLOPER XL_Tenor1,
															   LPXLOPER XL_Tenor2,
															   LPXLOPER XL_Expiry,
															   LPXLOPER XL_Moneyness);


// IndexIndexCorrelCube
__declspec(dllexport) LPXLOPER WINAPI Local_IndexIndexCorrelCube(LPXLOPER XL_correlCurves,
																 LPXLOPER XL_Tenors1List,
																 LPXLOPER XL_Tenors2List);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IndexIndexCorrelCube(LPXLOPER XL_correlCurves,
																	 LPXLOPER XL_Tenors1List,
																	 LPXLOPER XL_Tenors2List);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeIndexIndexCorrel(LPXLOPER XL_correlCube,
																	LPXLOPER XL_Tenor1,
																	LPXLOPER XL_Tenor2,
																	LPXLOPER XL_Expiry1,
																	LPXLOPER XL_Expiry2);

__declspec(dllexport) LPXLOPER WINAPI Local_CreateCorrelator( LPXLOPER XL_mktTags,
															  LPXLOPER XL_HyperDiagCurveIds,
															  LPXLOPER XL_IndexIndexCurveIds,
															  LPXLOPER XL_CorrelCurveIds,
															  LPXLOPER XL_IndexCurveIds,
															  LPXLOPER XL_IRVolHyperCubeIds,
															  LPXLOPER XL_VolVolHyperCubeIds,
															  LPXLOPER XL_FXVolHyperCubeIds);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateCorrelator( LPXLOPER XL_mktTags,
															      LPXLOPER XL_HyperDiagCurveIds,
															      LPXLOPER XL_IndexIndexCurveIds,
															      LPXLOPER XL_CorrelCurveIds,
															      LPXLOPER XL_IndexCurveIds,
																  LPXLOPER XL_IRVolHyperCubeIds,
																  LPXLOPER XL_VolVolHyperCubeIds,
																  LPXLOPER XL_FXVolHyperCubeIds);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeIdIdCorrelFromManager(LPXLOPER XL_correlManagerId,
																		 LPXLOPER XL_tenor1,
																		 LPXLOPER XL_tenor2,
																		 LPXLOPER XL_expiry1,
																		 LPXLOPER XL_expiry2,
																		 LPXLOPER XL_currecy);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeHyperCorrelFromManager(LPXLOPER XL_correlManagerId,
																		  LPXLOPER XL_tenor1,
																		  LPXLOPER XL_tenor2,
																		  LPXLOPER XL_expiry,
																		  LPXLOPER XL_currecy,
																		  LPXLOPER XL_byExpiry);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrSimpleFromManager(LPXLOPER XL_correlManagerId,
															      LPXLOPER XL_tenor,
															      LPXLOPER XL_expiry,
																  LPXLOPER XL_currency,
																  LPXLOPER XL_Corr);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeFXCorrelFromManager(LPXLOPER XL_correlManagerId,
																	   LPXLOPER XL_currency1,
																	   LPXLOPER XL_currency2,
																	   LPXLOPER XL_tenor,
																	   LPXLOPER XL_expiry);

__declspec(dllexport) LPXLOPER WINAPI Local_BumpVolatilityCorrelManager(LPXLOPER XL_CorrelManager,
																		LPXLOPER XL_Value,
																		LPXLOPER XL_nthLine,
																		LPXLOPER XL_nthCol,
																		LPXLOPER XL_isCumul,
																		LPXLOPER XL_isAbsolute,
																		LPXLOPER XL_TypeCorr,
																		LPXLOPER XL_Currency/*,
																		LPXLOPER XL_isToClone*/);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BumpVolatilityCorrelManager(LPXLOPER XL_CorrelManager,
																		    LPXLOPER XL_Value,
																		    LPXLOPER XL_nthLine,
																		    LPXLOPER XL_nthCol,
																		    LPXLOPER XL_isCumul,
																		    LPXLOPER XL_isAbsolute,
																		    LPXLOPER XL_TypeCorr,
																		    LPXLOPER XL_Currency/*,
																		    LPXLOPER XL_isToClone*/);


__declspec(dllexport) LPXLOPER WINAPI Local_GetCorrelDiagFromIndex(LPXLOPER XL_IndexIndexCubeId,
																   LPXLOPER XL_Tenor,
														      	   LPXLOPER XL_ListTenor);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetCorrelDiagFromIndex(LPXLOPER XL_IndexIndexCubeId,
																       LPXLOPER XL_Tenor,
														      	       LPXLOPER XL_ListTenor);

__declspec(dllexport) LPXLOPER WINAPI Local_GetMixtureParamsFromSummit ( LPXLOPER XL_index,
																		 LPXLOPER XL_currency,
																		 LPXLOPER XL_cvName,
																		 LPXLOPER XL_date,
																		 LPXLOPER XL_interpolMeth );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetMixtureParamsFromSummit ( LPXLOPER XL_index,
																			 LPXLOPER XL_currency,
																			 LPXLOPER XL_cvName,
																			 LPXLOPER XL_date,
																			 LPXLOPER XL_interpolMeth );

//----------------------------------------------------------------------------------//
// ARM_SABRVol Functions															//
//----------------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_SABRVol(LPXLOPER XL_SigmaOrAlphaId,
													LPXLOPER XL_RhoId,
													LPXLOPER XL_NuId,
													LPXLOPER XL_BetaId,
													LPXLOPER XL_SOrAFlag,
													LPXLOPER XL_ModelType,
													LPXLOPER XL_Weight);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABRVol(LPXLOPER XL_SigmaOrAlphaId,
														LPXLOPER XL_RhoId,
														LPXLOPER XL_NuId,
														LPXLOPER XL_BetaId,
														LPXLOPER XL_SOrAFlag,
														LPXLOPER XL_ModelType,
														LPXLOPER XL_Weight);


__declspec(dllexport) LPXLOPER WINAPI Local_GetSABRVolFromSummit(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
															     LPXLOPER XL_matIndex,
															     LPXLOPER XL_impOrHist,
															     LPXLOPER XL_indexid,
																 LPXLOPER XL_sigmaOrAlpha,
																 LPXLOPER XL_modelType,
																 LPXLOPER XL_weight);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetSABRVolFromSummit(LPXLOPER XL_index,
																	 LPXLOPER XL_currency,
																	 LPXLOPER XL_cvName,
																	 LPXLOPER XL_date,
																	 LPXLOPER XL_vtype,
																	 LPXLOPER XL_matIndex,
																	 LPXLOPER XL_impOrHist,
																	 LPXLOPER XL_indexid,
																	 LPXLOPER XL_sigmaOrAlpha,
																	 LPXLOPER XL_modelType,
																	 LPXLOPER XL_weight);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ConvIndexInYearTerm(LPXLOPER XL_index,
															           LPXLOPER XL_asof,
															           LPXLOPER XL_ccy);

#endif /* ARM_XL_VOLCRV_LOCAL_H */
