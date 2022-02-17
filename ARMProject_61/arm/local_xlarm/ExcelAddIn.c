
//
// Includes
//
#include <windows.h>
#include "ExcelAddIn.h"
#define CATEGORY L"New Add-In"
//
// Later, the instance handle is required to create dialog boxes.
// g_hInst holds the instance handle passed in by DllMain so that it is
// available for later use. hWndMain is used in several routines to
// store Microsoft Excel's hWnd. This is used to attach dialog boxes as
// children of Microsoft Excel's main window. A buffer is used to store
// the free space that DIALOGMsgProc will put into the dialog box.
//
//
// Global Variables
//
HWND g_hWndMain = NULL;
HANDLE g_hInst = NULL;
XCHAR g_szBuffer[20] = L"";

//
// Syntax of the Register Command:
//      REGISTER(module_text, procedure, type_text, function_text, 
//               argument_text, macro_type, category, shortcut_text,
//               help_topic, function_help, argument_help1, argument_help2,...)
//
//
// g_rgWorksheetFuncs will use only the first 11 arguments of 
// the Register function.
//
// This is a table of all the worksheet functions exported by this module.
// These functions are all registered (in xlAutoOpen) when you
// open the XLL. Before every string, leave a space for the
// byte count. The format of this table is the same as 
// arguments two through eleven of the REGISTER function.
// g_rgWorksheetFuncsRows define the number of rows in the table. The
// g_rgWorksheetFuncsCols represents the number of columns in the table.
//
#define g_rgWorksheetFuncsRows 1
#define g_rgWorksheetFuncsCols 10

static LPWSTR g_rgWorksheetFuncs
[g_rgWorksheetFuncsRows][g_rgWorksheetFuncsCols] =
{
	{ L"FuncFib",
		L"UU",	
		L"FuncFib",
		L"Compute to...",
		L"1",
		CATEGORY,
		L"",
		L"",
		L"Number to compute to"
		L"Computes the nth fibonacci number",
	},
};

//
// g_rgCommandFuncs
//
// This is a table of all the command functions exported by this module.
// These functions are all registered (in xlAutoOpen) when you
// open the XLL. Before every string, leave a space for the
// byte count. The format of this table is the same as 
// arguments two through eight of the REGISTER function.
// g_rgFuncsRows define the number of rows in the table. The
// g_rgCommandFuncsCols represents the number of columns in the table.
//
#define g_rgCommandFuncsRows 4
#define g_rgCommandFuncsCols 7

static LPWSTR g_rgCommandFuncs[g_rgCommandFuncsRows][g_rgCommandFuncsCols] =
{
	{ L"fDialog",                   // Procedure
		L"A",                   // type_text
		L"fDialog",             // function_text
		L"",                    // argument_text
		L"2",                   // macro_type
		CATEGORY,      // category
		L"l"                    // shortcut_text
	},
	{ L"fDance",
		L"A",
		L"fDance",
		L"",
		L"2",
		CATEGORY,
		L"m"
	},
	{ L"fShowDialog",
		L"A",
		L"fShowDialog",
		L"",
		L"2",
		CATEGORY,
		L"n"},
	{ L"fExit",
		L"A",
		L"fExit",
		L"",
		L"2",
		CATEGORY,
		L"o"
	},
};

//
// g_rgMenu
//
// This is a table describing the Generic drop-down menu. It is in
// the same format as the Microsoft Excel macro language menu tables.
// The first column contains the name of the menu or command, the
// second column contains the function to be executed, the third
// column contains the (Macintosh only) shortcut key, the fourth
// column contains the help text for the status bar, and
// the fifth column contains the help text index. Leave a space
// before every string so the byte count can be inserted. g_rgMenuRows
// defines the number of menu items. 5 represents the number of
// columns in the table.
//

#define g_rgMenuRows 5
#define g_rgMenuCols 5

static LPWSTR g_rgMenu[g_rgMenuRows][g_rgMenuCols] =
{
	{L"&Generic",          L"",            L"",
		L"The Generic XLL Add-In",         	L""},
	{L"&Dialog...",        L"fDialog",     L"",
		L"Run a sample generic dialog",    	L""},
	{L"D&ance",            L"fDance",      L"",
		L"Make the selection dance around",	L""},
	{L"&Native Dialog...", L"fShowDialog", L"",
		L"Run a sample native dialog",     	L""},
	{L"E&xit",             L"fExit",       L"",
		L"Close the Generic XLL",          	L""},
};

//
// g_rgTool
//
// This is a table describing the toolbar. It is in the same format
// as the Microsoft Excel macro language toolbar tables. The first column
// contains the ID of the tool, the second column contains the function
// to be executed, the third column contains a logical value specifying
// the default image of the tool, the fourth column contains a logical
// value specifying whether the tool can be used, the fifth column contains
// a face for the tool, the sixth column contains the help_text that
// is displayed in the status bar, the seventh column contains the Balloon
// text (Macintosh Only), and the eighth column contains the help topics
// as quoted text. Leave a space before every string so the byte count
// can be inserted. g_rgToolRows defines the number of tools on the toolbar.
// 8 represents the number of columns in the table.
//

#define g_rgToolRows 3
#define g_rgToolCols 8

static LPWSTR g_rgTool[g_rgToolRows][g_rgToolCols] =
{
	{L"211", L"fDance", L"FALSE", L"TRUE", L"", L"Dance the Selection", L"", L""},
	{L"0",   L"",       L"",      L"",     L"", L"",                    L"", L""},
	{L"212", L"fExit",  L"FALSE", L"TRUE", L"", L"Exit this example",   L"", L""},
};



///***************************************************************************
// DllMain()
//
// Purpose:
//
//      Windows calls DllMain, for both initialization and termination.
//		It also makes calls on both a per-process and per-thread basis,
//		so several initialization calls can be made if a process is multithreaded.
//
//      This function is called when the DLL is first loaded, with a dwReason
//      of DLL_PROCESS_ATTACH.
//
// Parameters:
//
//      HANDLE hDLL         Module handle.
//      DWORD dwReason,     Reason for call
//      LPVOID lpReserved   Reserved
//
// Returns: 
//      The function returns TRUE (1) to indicate success. If, during
//      per-process initialization, the function returns zero, 
//      the system cancels the process.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

BOOL APIENTRY DllMain( HANDLE hDLL, 
					   DWORD dwReason, 
					   LPVOID lpReserved )
{
	switch (dwReason)
	{
	case DLL_PROCESS_ATTACH:

		// The instance handle passed into DllMain is saved
		// in the global variable g_hInst for later use.

		g_hInst = hDLL;
		break;
	case DLL_PROCESS_DETACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	default:
		break;
	}
	return TRUE;
}


///***************************************************************************
// xlAutoOpen()
//
// Purpose: 
//      Microsoft Excel call this function when the DLL is loaded.
//
//      Microsoft Excel uses xlAutoOpen to load XLL files.
//      When you open an XLL file, the only action
//      Microsoft Excel takes is to call the xlAutoOpen function.
//
//      More specifically, xlAutoOpen is called:
//
//       - when you open this XLL file from the File menu,
//       - when this XLL is in the XLSTART directory, and is
//         automatically opened when Microsoft Excel starts,
//       - when Microsoft Excel opens this XLL for any other reason, or
//       - when a macro calls REGISTER(), with only one argument, which is the
//         name of this XLL.
//
//      xlAutoOpen is also called by the Add-in Manager when you add this XLL 
//      as an add-in. The Add-in Manager first calls xlAutoAdd, then calls
//      REGISTER("EXAMPLE.XLL"), which in turn calls xlAutoOpen.
//
//      xlAutoOpen should:
//
//       - register all the functions you want to make available while this
//         XLL is open,
//
//       - add any menus or menu items that this XLL supports,
//
//       - perform any other initialization you need, and
//
//       - return 1 if successful, or return 0 if your XLL cannot be opened.
//
// Parameters:
//
// Returns: 
//
//      int         1 on success, 0 on failure
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoOpen(void)
{

	static XLOPER12 xDLL,	   // name of this DLL //
	xMenu,	 // xltypeMulti containing the menu //
	xTool,	 // xltypeMulti containing the toolbar //
	xTest;	 // used for menu test //
	LPXLOPER12 pxMenu;	   // Points to first menu item //
	LPXLOPER12 px;		   // Points to the current item //
	LPXLOPER12 pxTool;	   // Points to first toolbar item //
	int i, j;			   // Loop indices //
	HANDLE   hMenu;		   // global memory holding menu //
	HANDLE  hTool;		   // global memory holding toolbar //

	//
	// In the following block of code the name of the XLL is obtained by
	// calling xlGetName. This name is used as the first argument to the
	// REGISTER function to specify the name of the XLL. Next, the XLL loops
	// through the g_rgWorksheetFuncs[] table, and the g_rgCommandFuncs[]
	// tableregistering each function in the table using xlfRegister. 
	// Functions must be registered before you can add a menu item.
	//
	
	Excel12f(xlGetName, &xDLL, 0);

	for (i=0; i<g_rgWorksheetFuncsRows; i++)
	{
		Excel12f(xlfRegister, 0, 1+ g_rgWorksheetFuncsCols,
			  (LPXLOPER12) &xDLL,
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][0]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][1]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][2]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][3]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][4]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][5]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][6]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][7]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][8]),
			  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][9]));
	}

	/*for (i=0; i<g_rgCommandFuncsRows; i++)
	{
		Excel12f(xlfRegister, 0, 1+ g_rgCommandFuncsCols,
			  (LPXLOPER12) &xDLL,
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][0]),
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][1]),
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][2]),
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][3]),
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][4]),
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][5]),
			  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][6]));
	}*/

	//
	// In the following block of code, the Generic drop-down menu is created.
	// Before creation, a check is made to determine if Generic already
	// exists. If not, it is added. If the menu needs to be added, memory is
	// allocated to hold the array of menu items. The g_rgMenu[] table is then
	// transferred into the newly created array. The array is passed as an
	// argument to xlfAddMenu to actually add the drop-down menu before the
	// help menu. As a last step the memory allocated for the array is
	// released.
	//
	// This block uses TempStr12() and TempNum12(). Both create a temporary
	// XLOPER12. The XLOPER12 created by TempStr12() contains the string passed to
	// it. The XLOPER12 created by TempNum12() contains the number passed to it.
	// The Excel12f() function frees the allocated temporary memory. Both
	// functions are part of the framework library.
	//

	Excel12f(xlfGetBar, &xTest, 3, TempInt12(10), TempStr12(L"AddinXLL"), TempInt12(0));

	if (xTest.xltype == xltypeErr)
	{
		hMenu = GlobalAlloc(GMEM_MOVEABLE,sizeof(XLOPER12) * g_rgMenuCols * g_rgMenuRows);
		px = pxMenu = (LPXLOPER12) GlobalLock(hMenu);

		for (i=0; i < g_rgMenuRows; i++)
		{
			for (j=0; j < g_rgMenuCols; j++)
			{
				px->xltype = xltypeStr;
				px->val.str = TempStr12(g_rgMenu[i][j])->val.str;
				px++;
			}
		}

		xMenu.xltype = xltypeMulti;
		xMenu.val.array.lparray = pxMenu;
		xMenu.val.array.rows = g_rgMenuRows;
		xMenu.val.array.columns = g_rgMenuCols;

		Excel12f(xlfAddMenu,0,3,TempNum12(10),(LPXLOPER12)&xMenu,TempStr12(L"Help"));

		GlobalUnlock(hMenu);
		GlobalFree(hMenu);
	}

	//
	// In this block of code, the Test toolbar is created. Before
	// creation, a check is made to ensure that Test doesn't already
	// exist. If it does not, it is created. Memory is allocated to hold
	// the array containing the toolbar. The information from the g_rgTool[]
	// table is then transferred into this array. The toolbar is added with
	// xlfAddToolbar and subsequently displayed with xlcShowToolbar. Finally,
	// the memory allocated for the toolbar and the XLL filename is released.
	//
	// This block uses TempInt12(), TempBool12(), and TempMissing12(). All three
	// create a temporary XLOPER12. The XLOPER12 created by TempInt() contains
	// the integer passed to it. TempBool12() creates an XLOPER12 containing the
	// boolean value passed to it. TempMissing12() creates an XLOPER12 that
	// simulates a missing argument. The Excel12f() function frees the temporary
	// memory associated with these functions. All three are part of the
	// framework library.
	//

	Excel12f(xlfGetToolbar, &xTest, 2, TempInt12(1), TempStr12(L"Test"));

	if (xTest.xltype == xltypeErr)
	{
		hTool = GlobalAlloc(GMEM_MOVEABLE, sizeof(XLOPER12) * g_rgToolCols * g_rgToolRows);
		px = pxTool = (LPXLOPER12) GlobalLock(hTool);

		for (i = 0; i < g_rgToolRows; i++)
		{
			for (j = 0; j < g_rgToolCols; j++)
			{
				px->xltype = xltypeStr;
				px->val.str = g_rgTool[i][j];
				px++;
			}
		}

		xTool.xltype = xltypeMulti;
		xTool.val.array.lparray = pxTool;
		xTool.val.array.rows = g_rgToolRows;
		xTool.val.array.columns = g_rgToolCols;

		Excel12f(xlfAddToolbar,0,2,TempStr12(L"Test"),(LPXLOPER)&xTool);

		Excel12f(xlcShowToolbar, 0, 6, TempStr12(L"Test"), TempBool12(1),
			  TempInt12(5), TempMissing12(), TempMissing12(), TempInt12(999));

		GlobalUnlock(hTool);
		GlobalFree(hTool);
	}

	// Free the XLL filename //
	Excel12f(xlFree, 0, 2, (LPXLOPER12) &xTest, (LPXLOPER12) &xDLL);

	return 1;
}


///***************************************************************************
// xlAutoClose()
//
// Purpose: Microsoft Excel call this function when the DLL is unloaded.
//
//      xlAutoClose is called by Microsoft Excel:
//
//       - when you quit Microsoft Excel, or 
//       - when a macro sheet calls UNREGISTER(), giving a string argument
//         which is the name of this XLL.
//
//      xlAutoClose is called by the Add-in Manager when you remove this XLL from
//      the list of loaded add-ins. The Add-in Manager first calls xlAutoRemove,
//      then calls UNREGISTER("GENERIC.XLL"), which in turn calls xlAutoClose.
// 
//      xlAutoClose is called by GENERIC.XLL by the function fExit. This function
//      is called when you exit Generic.
// 
//      xlAutoClose should:
// 
//       - Remove any menus or menu items that were added in xlAutoOpen,
// 
//       - do any necessary global cleanup, and
// 
//       - delete any names that were added (names of exported functions, and 
//         so on). Remember that registering functions may cause names to 
//         be created.
// 
//      xlAutoClose does NOT have to unregister the functions that were registered
//      in xlAutoOpen. This is done automatically by Microsoft Excel after
//      xlAutoClose returns.
// 
//      xlAutoClose should return 1.
//
// Parameters:
//
// Returns: 
//
//      int         1
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoClose(void)
{
	int i;
	XLOPER12 xRes;

	//
	// This block first deletes all names added by xlAutoOpen or
	// xlAutoRegister12. Next, it checks if the drop-down menu Generic still
	// exists. If it does, it is deleted using xlfDeleteMenu. It then checks
	// if the Test toolbar still exists. If it is, xlfDeleteToolbar is
	// used to delete it.
	//

	//
	// Due to a bug in Excel the following code to delete the defined names
	// does not work.  There is no way to delete these
	// names once they are Registered
	// The code is left in, in hopes that it will be
	// fixed in a future version.
	//

	for (i = 0; i < g_rgWorksheetFuncsRows; i++)
		Excel12f(xlfSetName, 0, 1, TempStr12(g_rgWorksheetFuncs[i][2]));

	for (i = 0; i < g_rgCommandFuncsRows; i++)
		Excel12f(xlfSetName, 0, 1, TempStr12(g_rgCommandFuncs[i][2]));
	//
	// Everything else works as documented
	//
	Excel12f(xlfGetBar, &xRes, 3, TempInt12(10), TempStr12(L"Generic"), TempInt12(0));

	if (xRes.xltype != xltypeErr)
	{
		Excel12f(xlfDeleteMenu, 0, 2, TempNum12(10), TempStr12(L"Generic"));

		// Free the XLOPER12 returned by xlfGetBar //
		Excel12f(xlFree, 0, 1, (LPXLOPER12) &xRes);
	}

	Excel12f(xlfGetToolbar, &xRes, 2, TempInt12(7), TempStr12(L"Test"));

	if (xRes.xltype != xltypeErr)
	{
		Excel12f(xlfDeleteToolbar, 0, 1, TempStr12(L"Test"));

		// Free the XLOPER12 returned by xlfGetToolbar //
		Excel12f(xlFree, 0, 1, (LPXLOPER12) &xRes);
	}

	return 1;
}


///***************************************************************************
// lpwstricmp()
//
// Purpose: 
//
//      Compares a pascal string and a null-terminated C-string to see
//      if they are equal.  Method is case insensitive
//
// Parameters:
//
//      LPWSTR s    First string (null-terminated)
//      LPWSTR t    Second string (byte counted)
//
// Returns: 
//
//      int         0 if they are equal
//                  Nonzero otherwise
//
// Comments:
//
//      Unlike the usual string functions, lpwstricmp
//      doesn't care about collating sequence.
//
// History:  Date       Author        Reason
///***************************************************************************

int lpwstricmp(LPWSTR s, LPWSTR t)
{
	int i;

	if (wcslen(s) != *t)
		return 1;

	for (i = 1; i <= s[0]; i++)
	{
		if (towlower(s[i-1]) != towlower(t[i]))
			return 1;
	}										  
	return 0;
}


///***************************************************************************
// xlAutoRegister12()
//
// Purpose:
//
//      This function is called by Microsoft Excel if a macro sheet tries to
//      register a function without specifying the type_text argument. If that
//      happens, Microsoft Excel calls xlAutoRegister12, passing the name of the
//      function that the user tried to register. xlAutoRegister12 should use the
//      normal REGISTER function to register the function, only this time it must
//      specify the type_text argument. If xlAutoRegister12 does not recognize the
//      function name, it should return a #VALUE! error. Otherwise, it should
//      return whatever REGISTER returned.
//
// Parameters:
//
//      LPXLOPER12 pxName   xltypeStr containing the
//                          name of the function
//                          to be registered. This is not
//                          case sensitive.
//
// Returns: 
//
//      LPXLOPER12          xltypeNum containing the result
//                          of registering the function,
//                          or xltypeErr containing #VALUE!
//                          if the function could not be
//                          registered.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER12 WINAPI xlAutoRegister12(LPXLOPER12 pxName)
{
	static XLOPER12 xDLL, xRegId;
	int i;

	//
	// This block initializes xRegId to a #VALUE! error first. This is done in
	// case a function is not found to register. Next, the code loops through 
	// the functions in g_rgFuncs[] and uses lpwstricmp to determine if the 
	// current row in g_rgFuncs[] represents the function that needs to be 
	// registered. When it finds the proper row, the function is registered 
	// and the register ID is returned to Microsoft Excel. If no matching 
	// function is found, an xRegId is returned containing a #VALUE! error.
	//

	xRegId.xltype = xltypeErr;
	xRegId.val.err = xlerrValue;


	for (i=0; i<g_rgWorksheetFuncsRows; i++)
	{
		if (!lpwstricmp(g_rgWorksheetFuncs[i][0], pxName->val.str))
		{
			Excel12f(xlfRegister, 0, 1+ g_rgWorksheetFuncsCols,
				  (LPXLOPER12) &xDLL,
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][0]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][1]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][2]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][3]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][4]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][5]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][6]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][7]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][8]),
				  (LPXLOPER12) TempStr12(g_rgWorksheetFuncs[i][9]));
			/// Free oper returned by xl //
			Excel12f(xlFree, 0, 1, (LPXLOPER12) &xDLL);

			return(LPXLOPER12) &xRegId;
		}
	}

	for (i=0; i<g_rgCommandFuncsRows; i++)
	{
		if (!lpwstricmp(g_rgCommandFuncs[i][0], pxName->val.str))
		{
			Excel12f(xlfRegister, 0, 1+ g_rgCommandFuncsCols,
				  (LPXLOPER12) &xDLL,
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][0]),
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][1]),
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][2]),
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][3]),
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][4]),
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][5]),
				  (LPXLOPER12) TempStr12(g_rgCommandFuncs[i][6]));
			/// Free oper returned by xl //
			Excel12f(xlFree, 0, 1, (LPXLOPER12) &xDLL);

			return(LPXLOPER12) &xRegId;
		}
	}     
	return(LPXLOPER12) &xRegId;
}

///***************************************************************************
// xlAutoAdd()
//
// Purpose:
//
//      This function is called by the Add-in Manager only. When you add a
//      DLL to the list of active add-ins, the Add-in Manager calls xlAutoAdd()
//      and then opens the XLL, which in turn calls xlAutoOpen.
//
// Parameters:
//
// Returns: 
//
//      int         1
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoAdd(void)
{
	XCHAR szBuf[255];

	wsprintfW((LPWSTR)szBuf, L"Thank you for adding GENERIC.XLL\n "
			 L"built on %hs at %hs", __DATE__, __TIME__);

	// Display a dialog box indicating that the XLL was successfully added //
	Excel12f(xlcAlert, 0, 2, TempStr12(szBuf), TempInt12(2));
	return 1;
}

///***************************************************************************
// xlAutoRemove()
//
// Purpose:
//
//      This function is called by the Add-in Manager only. When you remove
//      an XLL from the list of active add-ins, the Add-in Manager calls
//      xlAutoRemove() and then UNREGISTER("GENERIC.XLL").
//   
//      You can use this function to perform any special tasks that need to be
//      performed when you remove the XLL from the Add-in Manager's list
//      of active add-ins. For example, you may want to delete an
//      initialization file when the XLL is removed from the list.
//
// Parameters:
//
// Returns: 
//
//      int         1
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoRemove(void)
{
	// Show a dialog box indicating that the XLL was successfully removed //
	Excel12f(xlcAlert, 0, 2, TempStr12(L"Thank you for removing GENERIC.XLL!"),
		  TempInt12(2));
	return 1;
}

///***************************************************************************
// xlAddInManagerInfo12()
//
// Purpose:
//
//      This function is called by the Add-in Manager to find the long name
//      of the add-in. If xAction = 1, this function should return a string
//      containing the long name of this XLL, which the Add-in Manager will use
//      to describe this XLL. If xAction = 2 or 3, this function should return
//      #VALUE!.
//
// Parameters:
//
//      LPXLOPER12 xAction  What information you want. One of:
//                            1 = the long name of the
//                                add-in
//                            2 = reserved
//                            3 = reserved
//
// Returns: 
//
//      LPXLOPER12          The long name or #VALUE!.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER12 WINAPI xlAddInManagerInfo12(LPXLOPER12 xAction)
{
	static XLOPER12 xInfo, xIntAction;

	//
	// This code coerces the passed-in value to an integer. This is how the
	// code determines what is being requested. If it receives a 1, 
	// it returns a string representing the long name. If it receives 
	// anything else, it returns a #VALUE! error.
	//

	Excel12f(xlCoerce, &xIntAction, 2, xAction, TempInt12(xltypeInt));

	if (xIntAction.val.w == 1)
	{
		xInfo.xltype = xltypeStr;
		xInfo.val.str = L"\022Generic Standalone DLL";
	}
	else
	{
		xInfo.xltype = xltypeErr;
		xInfo.val.err = xlerrValue;
	}

	//Word of caution - returning static XLOPERs/XLOPER12s is not thread safe
	//for UDFs declared as thread safe, use alternate memory allocation mechanisms
	return(LPXLOPER12) &xInfo;
}

///***************************************************************************
// DIALOGMsgProc()
//
// Purpose:
//
//     This procedure is associated with the native Windows dialog box that
//     fShowDialog displays. It provides the service routines for the events
//     (messages) that occur when the user operates one of the dialog
//     box's buttons, entry fields, or controls.
//
// Parameters:
//
//      HWND hWndDlg        Contains the HWND of the dialog box
//      UINT message        The message to respond to
//      WPARAM wParam       Arguments passed by Windows
//      LPARAM lParam
//
// Returns: 
//
//      INT_PTR                TRUE if message processed, FALSE if not.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

INT_PTR CALLBACK DIALOGMsgProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM lParam)
{

	//
	// This block is a very simple message loop for a dialog box. It checks for
	// only three messages. When it receives WM_INITDIALOG, it uses a buffer to
	// set a static text item to the amount of free space on the stack. When it
	// receives WM_CLOSE it posts a CANCEL message to the dialog box. When it
	// receives WM_COMMAND it checks if the OK button was pressed. If it was,
	// the dialog box is closed and control returned to fShowDialog.
	//

	switch (message)
	{
	
	case WM_INITDIALOG:
		SetDlgItemText(hWndDlg, FREE_SPACE, (LPSTR)g_szBuffer);
		break;

	case WM_CLOSE:
		PostMessage(hWndDlg, WM_COMMAND, IDCANCEL, 0L);
		break;

	case WM_COMMAND:
		switch (wParam)
		{
		case IDOK:
			EndDialog(hWndDlg, FALSE);
			break;
		}
		break;

	default:
		return FALSE;
	}

	return TRUE;
}

///***************************************************************************
// ExcelCursorProc()
//
// Purpose:
//
//      When a modal dialog box is displayed over Microsoft Excel's window, the
//      cursor is a busy cursor over Microsoft Excel's window. This WndProc traps
//      WM_SETCURSORs and changes the cursor back to a normal arrow.
//
// Parameters:
//
//      HWND hWndDlg        Contains the HWND Window
//      UINT message        The message to respond to
//      WPARAM wParam       Arguments passed by Windows
//      LPARAM lParam
//
// Returns: 
//
//      LRESULT             0 if message handled, otherwise the result of the
//                          default WndProc
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

// Create a place to store Microsoft Excel's WndProc address //
static WNDPROC g_lpfnExcelWndProc = NULL;

LRESULT CALLBACK ExcelCursorProc(HWND hwnd, 
                                 UINT wMsg, 
                                 WPARAM wParam, 
                                 LPARAM lParam)
{
	//
	// This block checks to see if the message that was passed in is a
	// WM_SETCURSOR message. If so, the cursor is set to an arrow; if not,
	// the default WndProc is called.
	//

	if (wMsg == WM_SETCURSOR)
	{
		SetCursor(LoadCursor(NULL, IDC_ARROW));
		return 0L;
	}
	else
	{
		return CallWindowProc(g_lpfnExcelWndProc, hwnd, wMsg, wParam, lParam);
	}
}

///***************************************************************************
// HookExcelWindow()
//
// Purpose:
//
//     This is the function that installs ExcelCursorProc so that it is
//     called before Microsoft Excel's main WndProc.
//
// Parameters:
//
//      HANDLE hWndExcel    This is a handle to Microsoft Excel's hWnd
//
// Returns: 
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

extern void FAR PASCAL HookExcelWindow(HANDLE hWndExcel)
{
	//
	// This block obtains the address of Microsoft Excel's WndProc through the
	// use of GetWindowLongPtr(). It stores this value in a global that can be
	// used to call the default WndProc and also to restore it. Finally, it
	// replaces this address with the address of ExcelCursorProc using
	// SetWindowLongPtr().
	//

	g_lpfnExcelWndProc = (WNDPROC) GetWindowLongPtr(hWndExcel, GWLP_WNDPROC);
	SetWindowLongPtr(hWndExcel, GWLP_WNDPROC, (LONG_PTR)(FARPROC) ExcelCursorProc);
}

///***************************************************************************
// UnhookExcelWindow()
//
// Purpose:
//
//      This is the function that removes the ExcelCursorProc that was
//      called before Microsoft Excel's main WndProc.
//
// Parameters:
//
//      HANDLE hWndExcel    This is a handle to Microsoft Excel's hWnd
//
// Returns: 
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

extern void FAR PASCAL UnhookExcelWindow(HANDLE hWndExcel)
{
	//
	// This function restores Microsoft Excel's default WndProc using
	// SetWindowLongPtr to restore the address that was saved into
	// g_lpfnExcelWndProc by HookExcelWindow(). It then sets g_lpfnExcelWndProc
	// to NULL.
	//

	SetWindowLongPtr(hWndExcel, GWLP_WNDPROC, (LONG_PTR) g_lpfnExcelWndProc);
	g_lpfnExcelWndProc = NULL;
}


///***************************************************************************
// GetHwnd()
//
// Purpose:
//
//      This function returns the hWnd of Excel's main window. 
//
// Parameters:
//
//      HWND * phwnd    Will contain Excel's hWnd
//
// Returns: 
//
//      BOOL            FALSE  Could not find Excel's hWnd
//                      TRUE   Found Excel's hWnd
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************
BOOL GetHwnd(HWND * pHwnd)
{
	XLOPER12 x;

	if (Excel12f(xlGetHwnd, &x, 0) == xlretSuccess)
	{
		*pHwnd = (HWND)x.val.w;
		return TRUE;
	}
	return FALSE;
}

///***************************************************************************
// FuncFib()
//
// Purpose:
//
//      A sample function that computes the nth Fibonacci number.
//      Features a call to several wrapper functions.
//
// Parameters:
//
//      LPXLOPER12 n    int to compute to
//
// Returns: 
//
//      LPXLOPER12      nth Fibonacci number
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER12 WINAPI FuncFib (LPXLOPER12 n)
{
	static XLOPER12 xResult;
	XLOPER12 xlt;
	int val, max, error = -1;
	int fib[2] = {1,1};
	switch (n->xltype)
	{
	case xltypeNum:
		max = (int)n->val.num;
		if (max < 0)
			error = xlerrValue;
		for (val = 3; val <= max; val++)
		{
			fib[val%2] += fib[(val+1)%2];
		}
		xResult.xltype = xltypeNum;
		xResult.val.num = fib[(val+1)%2];
		break;
	case xltypeSRef:
		error = Excel12f(xlCoerce, &xlt, 2, n, TempInt12(xltypeInt));
		if (!error)
		{
			error = -1;
			max = xlt.val.w;
			if (max < 0)
				error = xlerrValue;
			for (val = 3; val <= max; val++)
			{
				fib[val%2] += fib[(val+1)%2];
			}
			xResult.xltype = xltypeNum;
			xResult.val.num = fib[(val+1)%2];
		}
		Excel12f(xlFree, 0, 1, &xlt);
		break;
	default:
		error = xlerrValue;
		break;
	}

	if ( error != - 1 )
	{
		xResult.xltype = xltypeErr;
		xResult.val.err = error;
	}

	//Word of caution - returning static XLOPERs/XLOPER12s is not thread safe
	//for UDFs declared as thread safe, use alternate memory allocation mechanisms
    return(LPXLOPER12) &xResult;
}