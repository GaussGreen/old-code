/* -------------------------------------------------------------------------

   File: %M% 
   Path: %P% 
   Description: API XL de creation d'XLL
   Created: 99/12/02
   Author: EDK
   Modified: %E% %U% 
   Last maintained by: Charles-Emmanuel MUSY 
   Revision: %I% 

   -------------------------------------------------------------------------

   Note: Source: generic.c de l'Excel Developer's Kit

   Platform: Microsoft Windows

   ------------------------------------------------------------------------- */

#include <libCCcommon\CCcommon.h>
#include <libCCmessage\CCMessage.h>

#define XL_LOCAL_API_C
#include "XL_local_api.h"
//#include "ARM_local_license.h"

HWND g_hWndMain = NULL;
HANDLE g_hInst = NULL;

LPSTR lpTemp = " ";

BOOL WINAPI DllMain (HANDLE hDLL, DWORD dwReason, LPVOID lpReserved)
{
	int i, j;

// FIXMEFRED: mig.vc8 (25/05/2007 10:46:05): stderr is not a constant (see libCCmessage.c)
	MSG_set_stdlog(stderr);

	switch(dwReason)
   	{
    		case DLL_PROCESS_ATTACH:
      		{
				g_hInst = hDLL;
				for(i = 0; i < g_rgWorksheetFuncsRows; i++) 
	  			{	

	    			for(j = 0; j < g_rgWorksheetFuncsCols; j++) 
	      			{
						if (g_rgWorksheetFuncs[i][j] == NULL)
							g_rgWorksheetFuncs[i][j] = lpTemp;
						g_rgWorksheetFuncs[i][j][0] = (BYTE)lstrlen (g_rgWorksheetFuncs[i][j] + 1);
	      			}
	  			}

				break;
      		}
    		case DLL_PROCESS_DETACH:
    		case DLL_THREAD_ATTACH:
    		case DLL_THREAD_DETACH:
    		default:
      			break;
    	}

  	return TRUE;
}

__declspec(dllexport) int WINAPI xlAutoOpen (void)
{
  	static XLOPER xDLL,		// name of this DLL //
  	xMenu,   				// xltypeMulti containing the menu //
  	xTool,   				// xltypeMulti containing the toolbar //
  	xTest;   				// used for menu test //
  	LPXLOPER pxMenu;       	// Points to first menu item //
  	LPXLOPER px;           	// Points to the current item //
  	int i, j, k;              	// Loop indices //
  	HANDLE hMenu;        	// global memory holding menu //
#ifndef NO_RG_TOOL
  	LPXLOPER pxTool;       	// Points to first toolbar item //
  	HANDLE hTool;         	// global memory holding toolbar //
#endif 	/* NO_RG_TOOL */
  	HMODULE hmodDLL1,hmodDLL2;
  	char szFile1[256], szFile2[256];
  	DWORD dwRet; 

#ifdef XLLOCALARM
  	hmodDLL1 = GetModuleHandle (XLLOCALARM_XLL_NAME);
#elif defined(MERCURE_ARM_XLL)
	hmodDLL1 = GetModuleHandle (MERCURE_ARM_XLL_NAME);
#else
  	hmodDLL1 = NULL;
#endif /* XLLOCALARM */

  	hmodDLL2 = GetModuleHandle (NULL);

  	dwRet = GetModuleFileName (hmodDLL1, szFile1, 256);
  	dwRet = GetModuleFileName (hmodDLL2, szFile2, 256);

  	Excel (xlfGetBar, &xTest, 3, TempInt (10), TempStr (TITLE_STR), TempInt (0));

  	if(xTest.xltype != xltypeErr) 
	{
		static char szBuf[255];

		wsprintf ((LPSTR)szBuf, " You are already using Local ARM \n Close Excel and Reopen your XLL");

		Excel (xlcAlert, 0, 2, TempStr (szBuf), TempInt (2));

  		Excel (xlFree, 0, 2, (LPXLOPER)&xTest, (LPXLOPER)&xDLL);

  		return 1;
	}

	k = 1;

	if (k==1)
	{
//	#ifdef XLLOCALARM 
  		XLLOCALARM_xlAutoOpen ();
//	#endif /* XLLOCALARM */

  		Excel (xlGetName, &xDLL, 0);

  		for(i = 0; i < g_rgWorksheetFuncsRows; i++) 
		{
			Excel(xlfRegister, 0, 1 + g_rgWorksheetFuncsCols,
	   		(LPXLOPER)&xDLL,
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][0]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][1]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][2]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][3]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][4]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][5]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][6]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][7]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][8]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][9]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][10]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][11]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][12]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][13]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][14]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][15]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][16]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][17]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][18]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][19]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][20]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][21]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][22]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][23]),
			(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][24]),
			(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][25]),
			(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][26]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][27]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][28]));
/*	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][29]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][30]),
	   		(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][31]),
			(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][32]),
			(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][33]));*/
		}
  		for(i = 0; i < g_rgCommandFuncsRows; i++) 
		{
    		Excel (xlfRegister, 0, 1 + g_rgCommandFuncsCols,
	   		(LPXLOPER)&xDLL,
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][0]),
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][1]),
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][2]),
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][3]),
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][4]),
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][5]),
	   		(LPXLOPER)TempStr (g_rgCommandFuncs[i][6]));
   		}

  		Excel (xlfGetBar, &xTest, 3, TempInt (10), TempStr (TITLE_STR), TempInt (0));

  		if(xTest.xltype == xltypeErr) 
   		{
   			hMenu = GlobalAlloc (GMEM_MOVEABLE, sizeof (XLOPER) * 5 * g_rgMenuRows);
   			px = pxMenu = (LPXLOPER)GlobalLock (hMenu);

   			for(i = 0; i < g_rgMenuRows; i++) 
			{
	  			for(j = 0; j < g_rgMenuCols; j++) 
    			{
      				px->xltype = xltypeStr;
      				px->val.str = g_rgMenu[i][j];
      				px++;
    			}
			}

   			xMenu.xltype = xltypeMulti;
   			xMenu.val.array.lparray = pxMenu;
   			xMenu.val.array.rows = g_rgMenuRows;
   			xMenu.val.array.columns = 5;

   			Excel (xlfAddMenu, 0, 3, TempNum (10), (LPXLOPER)&xMenu, TempStr (" Help"));

   			GlobalUnlock (hMenu);
   			GlobalFree (hMenu);
   		}

	#ifndef NO_RG_TOOL
  		Excel (xlfGetToolbar, &xTest, 2, TempInt (1), TempStr (TITLE_STR));

  		if(xTest.xltype == xltypeErr) 
   		{
   			hTool = GlobalAlloc (GMEM_MOVEABLE, sizeof (XLOPER) * 8 * g_rgToolRows);
   			px = pxTool = (LPXLOPER)GlobalLock (hTool);

   			for(i = 0; i < g_rgToolRows; i++) 
			{
	  			for(j = 0; j < 8; j++) 
    			{
      				px->xltype = xltypeStr;
      				px->val.str = g_rgTool[i][j];
      				px++;
    			}
			}

   			xTool.xltype = xltypeMulti;
   			xTool.val.array.lparray = pxTool;
   			xTool.val.array.rows = g_rgToolRows;
   			xTool.val.array.columns = 8;

   			Excel (xlfAddToolbar, 0, 2, TempStr (TITLE_STR), (LPXLOPER)&xTool);

      		Excel (xlcShowToolbar, 0, 6, TempStr (TITLE_STR), TempBool (1), TempInt (5), TempMissing (), TempMissing (), TempInt (999));

      		GlobalUnlock (hTool);
      		GlobalFree (hTool);
  		}
	#endif 	/* NO_RG_TOOL */

	}
	else
	{
		static char szBuf[255];

		wsprintf ((LPSTR)szBuf, " You are not authorized to use %s \n ", 
					XLLOCALARM_XLL_DESCRIPTION);

		Excel (xlcAlert, 0, 2, TempStr (szBuf), TempInt (2));
	}

  	Excel (xlFree, 0, 2, (LPXLOPER)&xTest, (LPXLOPER)&xDLL);

  	return 1;
}

__declspec(dllexport) int WINAPI xlAutoClose (void)
{
  	int i;
  	XLOPER xRes;

  	for(i = 0; i < g_rgWorksheetFuncsRows; i++)
	{
    	Excel (xlfSetName, 0, 1, TempStr (g_rgWorksheetFuncs[i][2]));
	}

  	for(i = 0; i < g_rgCommandFuncsRows; i++)
	{
		Excel (xlfSetName, 0, 1, TempStr (g_rgCommandFuncs[i][2]));
	}

  	Excel (xlfGetBar, &xRes, 3, TempInt (10), TempStr (TITLE_STR), TempInt (0));

  	if(xRes.xltype != xltypeErr) 
    {
    	Excel (xlfDeleteMenu, 0, 2, TempNum (10), TempStr (TITLE_STR));
    	Excel (xlFree, 0, 1, (LPXLOPER)&xRes);
    }

#ifndef NO_RG_TOOL
  	Excel (xlfGetToolbar, &xRes, 2, TempInt (7), TempStr (TITLE_STR));

  	if(xRes.xltype != xltypeErr) 
    {
    	Excel (xlfDeleteToolbar, 0, 1, TempStr (TITLE_STR));
    	Excel (xlFree, 0, 1, (LPXLOPER)&xRes);
    }
#endif 	/* NO_RG_TOOL */

#ifdef XLLOCALARM 
 	XLLOCALARM_xlAutoClose ();
#endif /* XLLOCALARM */

  	return 1;
}


int lpstricmp (LPSTR s, LPSTR t)
{
  	int i;

  	if(*s != *t)
	{
    		return 1;
	}

  	for(i = 1; i <= s[0]; i++) 
    {
    	if(tolower(s[i]) != tolower(t[i]))
		{
			return 1;
		}
    }

  	return 0;
}

__declspec(dllexport) LPXLOPER WINAPI xlAutoRegister (LPXLOPER pxName)
{
  	static XLOPER xDLL, xRegId;
  	int i;

  	xRegId.xltype = xltypeErr;
  	xRegId.val.err = xlerrValue;

  	for(i = 0; i < g_rgWorksheetFuncsRows; i++) 
    {
    	if(!lpstricmp (g_rgWorksheetFuncs[i][0], pxName->val.str))	
		{
	  		Excel (xlfRegister, 0, 1 + g_rgWorksheetFuncsCols,
				(LPXLOPER)&xDLL,
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][0]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][1]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][2]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][3]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][4]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][5]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][6]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][7]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][8]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][9]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][10]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][11]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][12]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][13]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][14]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][15]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][16]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][17]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][18]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][19]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][20]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][21]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][22]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][23]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][24]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][25]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][26]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][27]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][28]));
/*				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][29]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][30]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][31]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][32]),
				(LPXLOPER)TempStr (g_rgWorksheetFuncs[i][33]));*/
				Excel (xlFree, 0, 1, (LPXLOPER)&xDLL);

	  		return (LPXLOPER)&xRegId;
		}
  	}

  	for(i = 0; i < g_rgCommandFuncsRows; i++) 
    {
    	if(!lpstricmp (g_rgCommandFuncs[i][0], pxName->val.str))	
		{
	  		Excel (xlfRegister, 0, 1 + g_rgCommandFuncsCols,
				(LPXLOPER)&xDLL,
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][0]),
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][1]),
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][2]),
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][3]),
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][4]),
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][5]),
				(LPXLOPER)TempStr (g_rgCommandFuncs[i][6]));
	    		Excel (xlFree, 0, 1, (LPXLOPER)&xDLL);

	  		return (LPXLOPER) &xRegId;
		}
	}

	return (LPXLOPER)&xRegId;
}

__declspec(dllexport) int WINAPI xlAutoAdd (void)
{        
  	static char szBuf[255];

#ifdef XLLOCALARM 
  	wsprintf ((LPSTR)szBuf, " Thank you for using%s revision %s\n ", 
	        XLLOCALARM_XLL_DESCRIPTION, XLLOCALARM_revision);
#endif 	/* XLLOCALARM */

  	Excel (xlcAlert, 0, 2, TempStr (szBuf), TempInt (2));

  	return 1;
}

__declspec(dllexport) int WINAPI xlAutoRemove (void)
{
  	static char szBuf[255];

#ifdef XLLOCALARM 
  	wsprintf ((LPSTR)szBuf, " Thank you for using%s revision %s\n ", 
	        XLLOCALARM_XLL_DESCRIPTION, XLLOCALARM_revision);
#endif 	/* XLLOCALARM */

  	Excel (xlcAlert, 0, 2, TempStr (szBuf), TempInt (2));

  	return 1;
}

__declspec(dllexport) LPXLOPER WINAPI xlAddInManagerInfo (LPXLOPER xAction)
{
  	static XLOPER xInfo, xIntAction;

  	Excel (xlCoerce, &xIntAction, 2, xAction, TempInt (xltypeInt));

  	if(xIntAction.val.w == 1) 
    {
    	xInfo.xltype = xltypeStr;
#ifdef XLLOCALARM 
	  	xInfo.val.str = XLLOCALARM_XLL_TOPIC;
#else
	  	xInfo.val.str = "\026Generic Standalone DLL";
#endif 	/* XLLOCALARM */
    }
  	else 
    {
    	xInfo.xltype = xltypeErr;
    	xInfo.val.err = xlerrValue;
    }

  	return (LPXLOPER)&xInfo;
}

static FARPROC g_lpfnExcelWndProc = NULL;

LRESULT CALLBACK ExcelCursorProc (HWND hwnd, UINT wMsg, WPARAM wParam, LPARAM lParam)
{
  	if (wMsg == WM_SETCURSOR) 
    {
    	SetCursor (LoadCursor (NULL, IDC_ARROW));
    	return 0L;
    }
  	else 
    {
    	return CallWindowProc ((WNDPROC)g_lpfnExcelWndProc, hwnd, wMsg, wParam, lParam);
    }
}

extern void FAR PASCAL HookExcelWindow (HANDLE hWndExcel)
{
  	g_lpfnExcelWndProc = (FARPROC)GetWindowLong (hWndExcel, GWL_WNDPROC);
  	SetWindowLong (hWndExcel, GWL_WNDPROC, (LONG)(FARPROC)ExcelCursorProc);
}

extern void FAR PASCAL UnhookExcelWindow (HANDLE hWndExcel)
{
  	SetWindowLong (hWndExcel, GWL_WNDPROC, (LONG)g_lpfnExcelWndProc);
  	g_lpfnExcelWndProc = NULL;
}

__declspec(dllexport) int WINAPI XLLOCALARM_Exit (void)
{
  	XLOPER xDLL, xFunc, xRegId;
  	int i;

	xFunc.xltype = xltypeStr;

  	Excel (xlGetName, &xDLL, 0);

  	for(i = 0; i < g_rgWorksheetFuncsRows; i++) 
    {
    	xFunc.val.str = (LPSTR)(g_rgWorksheetFuncs[i][0]);
    	Excel (xlfRegisterId, &xRegId, 2, (LPXLOPER)&xDLL, (LPXLOPER)&xFunc);
    	Excel (xlfUnregister, 0, 1, (LPXLOPER)&xRegId);
    }

  	Excel (xlFree, 0, 1, (LPXLOPER)&xDLL);

  	return xlAutoClose ();
}

/*
__declspec(dllexport) void WINAPI xlAutoFree (LPXLOPER px)
{
	// only for strings
	if(px)
	{
		if(px->val.str)
		{
			free (px->val.str);
			px->val.str = NULL;
		}

		if (px->xltype==xltypeMulti) 
		{
			unsigned int rows,cols; 
			unsigned int i; 
			//	Deleting an array.
			XLOPER* rgx=(XLOPER*)px->val.array.lparray ;
			if (!rgx) return ;			
			if (px->val.array.rows*px->val.array.columns==0) return; 
			rows = px->val.array.rows; 
			cols = px->val.array.columns;
			//	checking element type and deleteling if string
			for (i=0;i<rows*cols;i++) 
			{
				if (rgx[i].xltype==xltypeStr) 
				{
					char* ptr=rgx[i].val.str ;
					free( ptr) ;
				}
			}
			delete [] rgx ; 
		}
	}

	return;
}
*/

/* EOF %M% */
