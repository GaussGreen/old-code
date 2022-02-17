#ifndef XL_LOCAL_XLARM_INTERFACE_H
#define XL_LOCAL_XLARM_INTERFACE_H

#include <string.h>

#include "XL_local_api.h"
#include "XL_local_xlarm_common.h"

CCEXTERN_FUNCTION (int XLLOCALARM_xlAutoOpen, (void));
CCEXTERN_FUNCTION (int XLLOCALARM_xlAutoClose, (void));

//--------------------------------------------

#define g_rgWorksheetFuncsCols         	29
#define g_rgCommandFuncsCols            7
#define g_rgMenuCols                    5


extern LPSTR g_rgWorksheetFuncs[][g_rgWorksheetFuncsCols];
extern int g_rgWorksheetFuncsRows;

extern LPSTR g_rgCommandFuncs[][g_rgCommandFuncsCols];
extern int g_rgCommandFuncsRows;

extern LPSTR g_rgMenu[][g_rgMenuCols];
extern int g_rgMenuRows;


__declspec(dllexport) int WINAPI XLLOCALARM_xlCalculateNow (void);
__declspec(dllexport) int WINAPI XLLOCALARM_xlFirstCalculate (void);

extern BOOL GetHwnd (HWND* pHwnd);

#endif	/* XL_LOCAL_XLARM_INTERFACE_H */

/* EOF %M% */ 
