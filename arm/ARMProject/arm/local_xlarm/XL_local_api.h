/* -------------------------------------------------------------------------

   File: %M% 
   Path: %P% 
   Description: include du module XL d'acces aux fonctions Arm 
   Created: 99/12/02
   Author: Charles-Emmanuel MUSY 
   Modified: %E% %U% 
   Last maintained by: Charles-Emmanuel MUSY 
   Revision: %I% 

   -------------------------------------------------------------------------

   Note: Source: generic.c de l'Excel Developer's Kit

   Platform:    Microsoft Windows

   ------------------------------------------------------------------------- */

#ifndef XL_LOCAL_API_H
#define XL_LOCAL_API_H

#include <libCCcommon\CCcommon.h>

#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <xlcall.h>
#include <framewrk.h>
#ifdef __cplusplus
}
#endif

CCEXTERN_FUNCTION (void FAR PASCAL HookExcelWindow, (HANDLE hWndExcel));
CCEXTERN_FUNCTION (void FAR PASCAL UnhookExcelWindow, (HANDLE hWndExcel));
/*BOOL C_GetHwnd (HWND* pHwnd);*/

__declspec(dllexport) int WINAPI XLLOCALARM_Exit (void);

#ifdef __cplusplus
extern "C"
{
#endif	// __cplusplus
extern HWND g_hWndMain;
#ifdef __cplusplus
}
#endif	// __cplusplus

#ifdef __cplusplus
extern "C"
{
#endif	// __cplusplus
extern HANDLE g_hInst;
#ifdef __cplusplus
}
#endif	// __cplusplus

#ifdef XL_LOCAL_API_C

#include <winuser.h>

#ifdef XLLOCALARM
#include "XL_LOCAL_xlarm_interface.h"
#define TITLE_STR	XLLOCALARM_TITLE_STR
#endif 	/* XLLOCALARM */

#ifdef MERCURE_ARM_XLL
#include "XL_LOCAL_xlarm_interface.h"
#define TITLE_STR	MERCURE_ARM_XLL_TITLE_STR
#endif

#include <CCmessage.h>

extern char* XLLOCALARM_revision;

#endif 	/* XL_LOCAL_API_C */

extern int isAuthorized();

#endif 	/* XL_LOCAL_API_H */

/* EOF %M% */ 
