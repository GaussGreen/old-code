#ifndef ARM_XL_FOREX_LOCAL_H
#define ARM_XL_FOREX_LOCAL_H

__declspec(dllexport) LPXLOPER WINAPI Local_FOREX (LPXLOPER XL_LeftCcy,
												   LPXLOPER XL_RightCcy,
                                                   LPXLOPER XL_SpotValue);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FOREX (LPXLOPER XL_LeftCcy,
													   LPXLOPER XL_RightCcy,
                                                       LPXLOPER XL_SpotValue);

#endif	/* ARM_XL_FOREX_LOCAL_H */

/* EOF %M% */ 


