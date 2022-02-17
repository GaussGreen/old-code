#ifndef ARM_XL_REFVAL_LOCAL_H
#define ARM_XL_REFVAL_LOCAL_H


__declspec(dllexport) LPXLOPER WINAPI Local_REFVALUE (LPXLOPER XL_dates,
													  LPXLOPER XL_values,
													  LPXLOPER XL_valueType,
													  LPXLOPER XL_conversion,
													  LPXLOPER XL_calcMethod,
													  LPXLOPER XL_values2);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REFVALUE (LPXLOPER XL_dates,
														  LPXLOPER XL_values,
														  LPXLOPER XL_valueType,
														  LPXLOPER XL_conversion,
														  LPXLOPER XL_calcMethod,
														  LPXLOPER XL_values2);


__declspec(dllexport) LPXLOPER WINAPI Local_CONSTREFVALUE (LPXLOPER XL_value);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CONSTREFVALUE (LPXLOPER XL_value)

__declspec(dllexport) LPXLOPER WINAPI Local_IATHREELEVREFVAL (LPXLOPER XL_level1,
															  LPXLOPER XL_amort1,
															  LPXLOPER XL_level2,
															  LPXLOPER XL_amort2,
															  LPXLOPER XL_level3,
															  LPXLOPER XL_amort3,
															  LPXLOPER XL_notionnal);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IATHREELEVREFVAL (LPXLOPER XL_level1,
																  LPXLOPER XL_amort1,
																  LPXLOPER XL_level2,
																  LPXLOPER XL_amort2,
																  LPXLOPER XL_level3,
																  LPXLOPER XL_amort3,
																  LPXLOPER XL_notionnal);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CptRefValue (LPXLOPER XL_refval,
															 LPXLOPER XL_date);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SumRefValue (LPXLOPER XL_refval1,
															 LPXLOPER XL_refval2,
															 LPXLOPER XL_coef);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SumRefValue (LPXLOPER XL_refval1,
																 LPXLOPER XL_refval2,
																 LPXLOPER XL_coef);

#endif /* ARM_XL_REFVAL_LOCAL_H */
