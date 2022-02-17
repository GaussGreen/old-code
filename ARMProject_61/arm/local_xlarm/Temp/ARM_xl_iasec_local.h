#ifndef ARM_XL_IASEC_LOCAL_H
#define ARM_XL_IASEC_LOCAL_H



__declspec(dllexport) LPXLOPER WINAPI Local_IASEC (LPXLOPER XL_underlying,
												   LPXLOPER XL_iaCtrl,
												   LPXLOPER XL_refVal,
												   LPXLOPER XL_iaCtrlType);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IASEC (LPXLOPER XL_underlying,
													   LPXLOPER XL_iaCtrl,
													   LPXLOPER XL_refVal,
													   LPXLOPER XL_iaCtrlType);

#endif
