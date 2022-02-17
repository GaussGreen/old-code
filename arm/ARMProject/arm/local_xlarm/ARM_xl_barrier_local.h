
#ifndef ARM_XL_BARRIER_LOCAL_H
#define ARM_XL_BARRIER_LOCAL_H

#include "XL_local_api.h"

__declspec(dllexport) LPXLOPER WINAPI Local_CONSTBARRIER (LPXLOPER XL_underlying,
														  LPXLOPER XL_tAsset,
														  LPXLOPER XL_maturity,
														  LPXLOPER XL_barrier,
														  LPXLOPER XL_upDown,
														  LPXLOPER XL_inOut,
														  LPXLOPER XL_triggerVar,
														  LPXLOPER XL_rebate,
														  LPXLOPER XL_firstX);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CONSTBARRIER (LPXLOPER XL_underlying,
															  LPXLOPER XL_tAsset,
															  LPXLOPER XL_maturity,
															  LPXLOPER XL_barrier,
															  LPXLOPER XL_upDown,
															  LPXLOPER XL_inOut,
															  LPXLOPER XL_triggerVar,
															  LPXLOPER XL_rebate,
															  LPXLOPER XL_firstX);

__declspec(dllexport) LPXLOPER WINAPI Local_BARRIER (LPXLOPER XL_underlying,
													 LPXLOPER XL_tAsset,
													 LPXLOPER XL_xStyle,
													 LPXLOPER XL_refVal,
													 LPXLOPER XL_upDown,
													 LPXLOPER XL_inOut,
													 LPXLOPER XL_triggerVar,
													 LPXLOPER XL_rebate);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BARRIER (LPXLOPER XL_underlying,
														 LPXLOPER XL_tAsset,
														 LPXLOPER XL_xStyle,
														 LPXLOPER XL_refVal,
														 LPXLOPER XL_upDown,
														 LPXLOPER XL_inOut,
														 LPXLOPER XL_triggerVar,
														 LPXLOPER XL_rebate);


#endif