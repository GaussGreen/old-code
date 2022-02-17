#ifndef ARM_XL_IRFUT_LOCAL_H
#define ARM_XL_IRFUT_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_IRFUT (LPXLOPER XL_delivery, 
													LPXLOPER XL_underlying);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRFUT (LPXLOPER XL_delivery,
													LPXLOPER XL_underlying);

__declspec(dllexport) LPXLOPER WINAPI Local_THREE_MONTH_FUT (LPXLOPER XL_delivery,
															 LPXLOPER XL_market,
															 LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_THREE_MONTH_FUT (LPXLOPER XL_delivery,
																 LPXLOPER XL_market,
																 LPXLOPER XL_ccy);


#endif /* ARM_XL_IRFUT_LOCAL_H */