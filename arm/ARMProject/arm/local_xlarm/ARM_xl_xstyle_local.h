#ifndef ARM_XL_XSTYLE_LOCAL_H
#define ARM_XL_XSTYLE_LOCAL_H


__declspec(dllexport) LPXLOPER WINAPI Local_BERMUDANXSTYLE (LPXLOPER XL_xDates, LPXLOPER XL_xExpiryDates);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BERMUDANXSTYLE (LPXLOPER XL_xDates, LPXLOPER XL_xExpiryDates);

__declspec(dllexport) LPXLOPER WINAPI Local_EUROPEANXSTYLE (LPXLOPER XL_xdate);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EUROPEANXSTYLE (LPXLOPER XL_xdate);

__declspec(dllexport) LPXLOPER WINAPI Local_AMERICANXSTYLE (LPXLOPER XL_xStartDate,
															LPXLOPER XL_xEndDate);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AMERICANXSTYLE (LPXLOPER XL_xStartDate,
																LPXLOPER XL_xEndDate);

#endif /*ARM_XL_XSTYLE_LOCAL_H*/