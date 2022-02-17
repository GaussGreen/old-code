#ifndef ARM_XL_XXXPROJECT_LOCAL_H
#define ARM_XL_XXXPROJECT_LOCAL_H

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Price			(		LPXLOPER XL_secId,
																			LPXLOPER XL_mktId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Hedge_Create	(		LPXLOPER XL_secId,
																			LPXLOPER XL_scenId,
																			LPXLOPER XL_mktDt);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Hedge_GetData	(		LPXLOPER XL_HedgeId,
																			LPXLOPER XL_Key);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Scenario		(		LPXLOPER XL_Shift,
																			LPXLOPER XL_Currency,
																			LPXLOPER XL_Type_Scenario,
																			LPXLOPER XL_SubType_Scenario,
																			LPXLOPER XL_Stress_Order,
																			LPXLOPER XL_Relatif,																	
																			LPXLOPER XL_CumulInv,
																			LPXLOPER XL_Perturbative);	
												
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Scenari_Compose	(		LPXLOPER XL_scenId1,
																			LPXLOPER XL_scenId2);

__declspec(dllexport) LPXLOPER WINAPI Local_MktDatas_Create(				LPXLOPER XL_asOfDate,
																			LPXLOPER XL_mktDatasKeys,
																			LPXLOPER XL_mktDatasId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDatas_Create(			LPXLOPER XL_asOfDate,
																			LPXLOPER XL_mktDatasKeys,
																			LPXLOPER XL_mktDatasId);

#endif	/* ARM_XL_XXXPROJECT_LOCAL_H */
