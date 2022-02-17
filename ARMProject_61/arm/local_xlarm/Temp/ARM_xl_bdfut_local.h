#ifndef ARM_XL_BDFUT_H
#define ARM_XL_BDFUT_H

#include "XL_local_api.h"



__declspec(dllexport) LPXLOPER WINAPI Local_BDFUT (LPXLOPER XL_delivery,
												   LPXLOPER XL_underlyingId,
												   LPXLOPER XL_coupon,
												   LPXLOPER XL_cfVal);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BDFUT (LPXLOPER XL_delivery,
													   LPXLOPER XL_underlyingId,
													   LPXLOPER XL_coupon,
													   LPXLOPER XL_cfVal);

__declspec(dllexport) LPXLOPER WINAPI Local_GetConversionFactor (LPXLOPER XL_bdFutId,
																 LPXLOPER XL_factId);

__declspec(dllexport) LPXLOPER WINAPI Local_GetCheapest (LPXLOPER XL_bdFutId);

__declspec(dllexport) LPXLOPER WINAPI Local_GILT_NOTIONNAL_BUND (LPXLOPER XL_delivery,
																 LPXLOPER XL_underlyingId,
																 LPXLOPER XL_bdFutType,
																 LPXLOPER XL_market);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GILT_NOTIONNAL_BUND (LPXLOPER XL_delivery,
																	 LPXLOPER XL_underlyingId,
																	 LPXLOPER XL_bdFutType,
																	 LPXLOPER XL_market);

__declspec(dllexport) LPXLOPER WINAPI Local_GILT (LPXLOPER XL_delivery,
												  LPXLOPER XL_underlyingId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GILT (LPXLOPER XL_delivery,
													  LPXLOPER XL_underlyingId);

__declspec(dllexport) LPXLOPER WINAPI Local_NOTIONNAL (LPXLOPER XL_delivery,
													   LPXLOPER XL_underlyingId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NOTIONNAL (LPXLOPER XL_delivery,
														   LPXLOPER XL_underlyingId);

__declspec(dllexport) LPXLOPER WINAPI Local_BUND_LIFFE (LPXLOPER XL_delivery,
														LPXLOPER XL_underlyingId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BUND_LIFFE (LPXLOPER XL_delivery,
															LPXLOPER XL_underlyingId);

__declspec(dllexport) LPXLOPER WINAPI Local_BUND_DTB (LPXLOPER XL_delivery,
													  LPXLOPER XL_underlyingId);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BUND_DTB (LPXLOPER XL_delivery,
														  LPXLOPER XL_underlyingId);

__declspec(dllexport) LPXLOPER WINAPI Local_BDFUTBASKET (LPXLOPER XL_delivery,
														 LPXLOPER XL_underlyingId,
														 LPXLOPER XL_coupon,
														 LPXLOPER XL_futurePrice);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BDFUTBASKET (LPXLOPER XL_delivery,
															 LPXLOPER XL_underlyingId,
															 LPXLOPER XL_coupon,
															 LPXLOPER XL_futurePrice);



#endif	/* ARM_XL_BDFUT_H */

/* EOF %M% */ 
