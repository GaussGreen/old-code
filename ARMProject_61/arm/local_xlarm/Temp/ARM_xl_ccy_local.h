#ifndef ARM_XL_CCY_LOCAL_H
#define ARM_XL_CCY_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_ISOCCY (LPXLOPER XL_name);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ISOCCY (LPXLOPER XL_name);

__declspec(dllexport) LPXLOPER WINAPI Local_CCY (LPXLOPER XL_name, LPXLOPER XL_idCurve,
										   LPXLOPER XL_crossValue, LPXLOPER XL_dayCount);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CCY (LPXLOPER XL_name, LPXLOPER XL_idCurve, 
											   LPXLOPER XL_crossValue, LPXLOPER XL_dayCount);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInfoFromCcy (LPXLOPER XL_ccyId,
																LPXLOPER XL_type);

#endif	/* ARM_XL_CCY_LOCAL_H */
