#ifndef ARM_XL_PF_LOCAL_H
#define ARM_XL_PF_LOCAL_H



#include "XL_local_api.h"



__declspec(dllexport) LPXLOPER WINAPI Local_PF (LPXLOPER XL_insts,
												LPXLOPER XL_coeffs,
												LPXLOPER XL_marketPrices,
												LPXLOPER XL_precisions);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PF (LPXLOPER XL_insts,
													LPXLOPER XL_coeffs,
													LPXLOPER XL_marketPrices,
													LPXLOPER XL_precisions);

__declspec(dllexport) LPXLOPER WINAPI Local_PF_FILL (LPXLOPER XL_Assets,
												LPXLOPER XL_Weights,
												LPXLOPER XL_MktPrices,
												LPXLOPER XL_VegaPrecisions,
                                                LPXLOPER XL_PortfolioId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PF_FILL (LPXLOPER XL_Assets,
												LPXLOPER XL_Weights,
												LPXLOPER XL_MktPrices,
												LPXLOPER XL_VegaPrecisions,
                                                LPXLOPER XL_PortfolioId);

__declspec(dllexport) LPXLOPER WINAPI Local_PF_Merge (LPXLOPER XL_PortfolioIds);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PF_Merge (LPXLOPER XL_PortfolioIds);

__declspec(dllexport) LPXLOPER WINAPI Local_PFGYCSIGVARFIT (LPXLOPER XL_pf, 
                                                      LPXLOPER XL_curve,
                                                      LPXLOPER XL_matCurve,
                                                      LPXLOPER XL_in_min_meanRev,
                                                      LPXLOPER XL_in_max_meanRev,
                                                      LPXLOPER XL_in_min_vol,
                                                      LPXLOPER XL_in_max_vol,
                                                      LPXLOPER XL_in_precision_meanRev,
                                                      LPXLOPER XL_in_precision_vol);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFGYCSIGVARFIT (LPXLOPER XL_pf, 
                                                      LPXLOPER XL_curve,
                                                      LPXLOPER XL_matCurve,
                                                      LPXLOPER XL_in_min_meanRev,
                                                      LPXLOPER XL_in_max_meanRev,
                                                      LPXLOPER XL_in_min_vol,
                                                      LPXLOPER XL_in_max_vol,
                                                      LPXLOPER XL_in_precision_meanRev,
                                                      LPXLOPER XL_in_precision_vol);

__declspec(dllexport) LPXLOPER WINAPI Local_PFLOGDECVOLFIT(LPXLOPER XL_pfId, 
                                                     LPXLOPER XL_curveId,
                                                     LPXLOPER XL_proba,
                                                     LPXLOPER XL_accuracy,
                                                     LPXLOPER XL_shapeType,
                                                     LPXLOPER XL_decay,
                                                     LPXLOPER XL_slope,
                                                     LPXLOPER XL_asymptote,                                                    
                                                     LPXLOPER XL_matCurve,
                                                     LPXLOPER XL_volinit_vect,
                                                     LPXLOPER XL_coeff_vect);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFLOGDECVOLFIT(LPXLOPER XL_pfId, 
                                                     LPXLOPER XL_curveId,
                                                     LPXLOPER XL_proba,
                                                     LPXLOPER XL_accuracy,
                                                     LPXLOPER XL_shapeType,
                                                     LPXLOPER XL_decay,
                                                     LPXLOPER XL_slope,
                                                     LPXLOPER XL_asymptote,
                                                     LPXLOPER XL_matCurve,
                                                     LPXLOPER XL_volinit_vect,
                                                     LPXLOPER XL_coeff_vect);

__declspec(dllexport) LPXLOPER WINAPI Local_PFGYCSIGVARPENALFIT(LPXLOPER XL_pfId, 
                                                          LPXLOPER XL_curveId,
                                                          LPXLOPER XL_matCurve,
                                                          LPXLOPER XL_start_meanRev,
                                                          LPXLOPER XL_start_vol,
                                                          LPXLOPER XL_penal_vect,
                                                          LPXLOPER XL_coeff_vect,
                                                          LPXLOPER XL_accuracy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFGYCSIGVARPENALFIT(LPXLOPER XL_pfId, 
                                                          LPXLOPER XL_curveId,
                                                          LPXLOPER XL_matCurve,
                                                          LPXLOPER XL_start_meanRev,
                                                          LPXLOPER XL_start_vol,
                                                          LPXLOPER XL_penal_vect,
                                                          LPXLOPER XL_coeff_vect,
                                                          LPXLOPER XL_accuracy);

__declspec(dllexport) LPXLOPER WINAPI Local_PFINSTLOGDECVOLFIT(LPXLOPER XL_pfId,
															   LPXLOPER XL_secId,
															   LPXLOPER XL_curveId,
															   LPXLOPER XL_proba,
															   LPXLOPER XL_UsePFResetDates,
															   LPXLOPER XL_accuracy,
															   LPXLOPER XL_shapeType,
															   LPXLOPER XL_decay,
															   LPXLOPER XL_slope,
															   LPXLOPER XL_asymptote,
															   LPXLOPER XL_VolBSId,
															   LPXLOPER XL_matCurve,
															   LPXLOPER XL_volinit_vect,
															   LPXLOPER XL_coeff_vect);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFINSTLOGDECVOLFIT(LPXLOPER XL_pfId,
																   LPXLOPER XL_secId,
																   LPXLOPER XL_curveId,
																   LPXLOPER XL_proba,
																   LPXLOPER XL_UsePFResetDates,
																   LPXLOPER XL_accuracy,
																   LPXLOPER XL_shapeType,
																   LPXLOPER XL_decay,
																   LPXLOPER XL_slope,
																   LPXLOPER XL_asymptote,
																   LPXLOPER XL_VolBSId,
																   LPXLOPER XL_matCurve,
																   LPXLOPER XL_volinit_vect,
																   LPXLOPER XL_coeff_vect);


__declspec(dllexport) LPXLOPER WINAPI Local_PFMODFIT (LPXLOPER XL_modName,
                                                LPXLOPER XL_pf,
                                                LPXLOPER XL_settlement,
                                                LPXLOPER XL_curve,
                                                LPXLOPER XL_vList,
                                                LPXLOPER XL_fList,
                                                LPXLOPER XL_nagAlgo,
                                                LPXLOPER XL_step,
                                                LPXLOPER XL_horizon);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFMODFIT (LPXLOPER XL_modName,
                                                    LPXLOPER XL_pf,
                                                    LPXLOPER XL_settlement,
                                                    LPXLOPER XL_curve,
                                                    LPXLOPER XL_vList,
                                                    LPXLOPER XL_fList,
                                                    LPXLOPER XL_nagAlgo,
                                                    LPXLOPER XL_step,
                                                    LPXLOPER XL_horizon);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CalibrationPF (LPXLOPER XL_curve,
														   LPXLOPER XL_vol,
														   LPXLOPER XL_sec,
														   LPXLOPER XL_modeType,
														   LPXLOPER XL_potfolioType,
														   LPXLOPER XL_shift);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CalibrationPF (LPXLOPER XL_curve,
															   LPXLOPER XL_vol,
															   LPXLOPER XL_sec,
															   LPXLOPER XL_modeType,
															   LPXLOPER XL_potfolioType,
															   LPXLOPER XL_shift);




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_MRSCalibrationPF (LPXLOPER XL_curve,
														   LPXLOPER XL_vol,
														   LPXLOPER XL_sec,
														   LPXLOPER XL_portfolio,
														   LPXLOPER XL_freq,
														   LPXLOPER XL_atmFlag);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_MRSCalibrationPF (LPXLOPER XL_curve,
															   LPXLOPER XL_vol,
															   LPXLOPER XL_sec,
															   LPXLOPER XL_portfolio,
															   LPXLOPER XL_freq,
															   LPXLOPER XL_atmFlag);

__declspec(dllexport) LPXLOPER WINAPI Local_GetAssetFromPF(
	                                                LPXLOPER XL_PortfolioId,
                                                    LPXLOPER XL_AssetIdIndex);

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeAll(	LPXLOPER XL_pfId,
													    LPXLOPER XL_modelId);


#endif /* ARM_XL_PF_LOCAL_H */