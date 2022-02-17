///***************************************************************************
// File:		GENERIC.C
//
// Purpose:		Template for creating 32-bit XLLs for Microsoft Excel 97.
//
//              This file contains sample code you can use as 
//              a template for writing your own Microsoft Excel XLLs. 
//              An XLL is a DLL that stands alone, that is, you
//              can open it by choosing the Open command from the
//              File menu. This code demonstrates many of the features 
//              of the Microsoft Excel C API.
//              
//              When you open GENERIC.XLL, it
//              creates a new Generic menu with four
//              commands:
//              
//                  Dialog...          displays a Microsoft Excel dialog box
//                  Dance              moves the selection around
//                                     until you press ESC
//                  Native Dialog...   displays a Windows dialog box
//                  Exit               Closes GENERIC.XLL and
//                                     removes the menu
//              
//              GENERIC.XLL also provides two functions,
//              Func1 and FuncSum, which can be used whenever
//              GENERIC.XLL is open. They can also be
//              registered without loading all of GENERIC.XLL
//              simply by calling REGISTER("GENERIC.XLL","FUNCx").
//              
//              GENERIC.XLL can also be added with the
//              Add-in Manager.
//              
//              This file uses the framework library
//              (FRMWRK32.LIB).
// 
// Platform:    Microsoft Windows
//
// Functions:		
//              DllMain
//              xlAutoOpen
//              xlAutoClose
//              lpstricmp
//              xlAutoRegister
//              xlAutoRegister
//              xlAutoAdd
//              xlAutoRemove
//              xlAddInManagerInfo
//              DIALOGMsgProc
//              ExcelCursorProc
//              HookExcelWindow
//              UnhookExcelWindow
//              fShowDialog
//              EnumProc
//              GetHwnd
//              Func1
//              FuncSum
//              fDance
//              fDialog
//              fExit
//
///***************************************************************************

//
// Includes
//
#include <windows.h>
#include "..\..\include\xlcall.h"
#include "framewrk.h"
#include "generic.h"

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
char g_szBuffer[20] = "";
//
// Syntax of the Register Command:
// REGISTER(module_text, procedure, type_text, function_text, 
//          argument_text, macro_type, category, shortcut_text,
//          help_topic, function_help, argument_help1, argument_help2,...)
//
//
// g_rgWorksheetFuncs will use only the first 11 arguments of 
//                  the Register function.
//
// This is a table of all the worksheet functions exported by this module.
// These functions are all registered (in xlAutoOpen) when you
// open the XLL. Before every string, leave a space for the
// byte count. The format of this table is the same as 
// arguments two through eleven of the REGISTER function.
// g_rgWorksheetFuncsRows define the number of rows in the table. The
// g_rgWorksheetFuncsCols represents the number of columns in the table.
//
#define g_rgWorksheetFuncsRows 2
#define g_rgWorksheetFuncsCols 10

static LPSTR g_rgWorksheetFuncs
                [g_rgWorksheetFuncsRows][g_rgWorksheetFuncsCols] =
{
    { " Func1",                               // Procedure
      " RR",                                  // type_text
      " Func1",                               // function_text
      " Arg",                                 // argument_text
      " 1",                                   // macro_type
      " Generic Add-In",                      // category
      " ",                                    // shortcut_text
      " ",                                    // help_topic
      " Always returns the string 'Func1'",   // function_help
      " Argument ignored"                     // argument_help1
    },
    { " FuncSum",
      " RRRRRRRRRRRRRRRRRRRRRRRRRRRRRR", // up to 29 args//
      " FuncSum",
      " number1,number2,...",
      " 1",
      " Generic Add-In",
      " ",                                    
      " ",                                  
      " Adds the arguments",   
      " Number1,number2,... are 1 to 29 arguments for which you want to sum."                    
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

static LPSTR g_rgCommandFuncs[g_rgCommandFuncsRows][g_rgCommandFuncsCols] =
{
  { " fDialog",               // Procedure
    " A",                     // type_text
    " fDialog",               // function_text
    " ",                      // argument_text
    " 2",                     // macro_type
    " Generic Add-In",        // category
    " l"                      // shortcut_text
  },
  { " fDance",
    " A",
    " fDance",
    " ",
    " 2",
    " Generic Add-In",
    " m"
  },
  { " fShowDialog",
    " A",
    " fShowDialog",
    " ",
    " 2",
    " Generic Add-In",
    " n"},
  { " fExit",
    " A",
    " fExit",
    " ",
    " 2",
    " Generic Add-In",
    " o"
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

static LPSTR g_rgMenu[g_rgMenuRows][g_rgMenuCols] =
{
    {" &Generic",          " ",           " ",
    " The Generic XLL Add-In",         " "},
    {" &Dialog...",        " fDialog",    " ",
    " Run a sample generic dialog",    " "},
    {" D&ance",            " fDance",     " ",
    " Make the selection dance around"," "},
    {" &Native Dialog...", " fShowDialog"," ",
    " Run a sample native dialog",     " "},
    {" E&xit",             " fExit",      " ",
    " Close the Generic XLL",          " "},
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

static LPSTR g_rgTool[g_rgToolRows][g_rgToolCols] =
{
    {" 211", " fDance", " FALSE", " TRUE", " ", " Dance the Selection", " ", " "},
    {" 0",   " ",       " ", " ",     " ", " ",                    " ", " "},
    {" 212", " fExit",  " FALSE", " TRUE", " ", " Exit this example",   " ", " "},
};


//
// g_rgDialog
//
// This is a table describing the sample dialog box used in the fDialog()
// function. Admittedly, it would be more efficient to use ints for
// the first 5 columns, but that makes the code much more complicated.
// Storing the text in string tables is another method that could be used.
// Leave a space before every string so the byte count can be inserted.
// g_rgDialogRows determines the number of rows in the dialog box.
// 7 represents the number of columns in the table.
//

#define g_rgDialogRows 16
#define g_rgDialogCols 7

static LPSTR g_rgDialog[g_rgDialogRows][g_rgDialogCols] =
{
    {" ",   " ",    " ",    " 494", " 210", " Generic Sample Dialog", " "},
    {" 1",  " 330", " 174", " 88",  " ",    " OK",                    " "},
    {" 2",  " 225", " 174", " 88",  " ",    " Cancel",                " "},
    {" 5",  " 19",  " 11",  " ",    " ",    " &Name:",                " "},
    {" 6",  " 19",  " 29",  " 251", " ",    " ",                      " "},
    {" 14", " 305", " 15",  " 154", " 73",  " &College",              " "},
    {" 11", " ",    " ",    " ",    " ",    " ",                      " 1"},
    {" 12", " ",    " ",    " ",    " ",    " &Harvard",              " 1"},
    {" 12", " ",    " ",    " ",    " ",    " &Other",                " "},
    {" 5",  " 19",  " 50",  " ",    " ",    " &Reference:",           " "},
    {" 10", " 19",  " 67",  " 253", " ",    " ",                      " "},
    {" 14", " 209", " 93",  " 250", " 63",  " &Qualifications",       " "},
    {" 13", " ",    " ",    " ",    " ",    " &BA / BS",              " 1"},
    {" 13", " ",    " ",    " ",    " ",    " &MA / MS",              " 1"},
    {" 13", " ",    " ",    " ",    " ",    " &PhD / Other Grad",     " 0"},
    {" 15", " 19",  " 99",  " 160", " 96",  " GENERIC_List1",         " 1"},
};

///***************************************************************************
// DllMain()
//
// Purpose:
//
//   Windows 95 and Windows NT call one function, DllMain, for both 
//   initialization and termination. It also makes calls on both a 
//   per-process and per-thread basis, so several initialization calls
//   can be made if a process is multithreaded.
//   This function is called when the DLL is first loaded, with a dwReason
//   of DLL_PROCESS_ATTACH. In this example, we byte-count all the strings
//   in the preceding tables.
//
// Parameters:
//
//        HANDLE hDLL       - Module handle.
//        DWORD dwReason,   - Reason for call
//        LPVOID lpReserved - Reserved
//
// Returns: The function returns TRUE (1) to indicate success. If, during
//          per-process initialization, the function returns zero, 
//          the system cancels the process.
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
    {
      //
      // In the following five loops, the strings in g_rgWorksheetFuncs[], 
      // g_rgCommandFuncs[], g_rgTool[], g_rgMenu[] and g_rgDialog[] are 
      // byte-counted so that they won't need to be byte-counted later.
      // In addition, the instance handle passed into DllMain is saved
      // in the global variable g_hInst for later use.
      //
      int i,j;

   	  g_hInst = hDLL;
    
      for (i=0; i<g_rgWorksheetFuncsRows; i++) 
	  {
		  for (j=0; j<g_rgWorksheetFuncsCols; j++) 
		  {
	    	  g_rgWorksheetFuncs[i][j][0] = 
                (BYTE) lstrlen (g_rgWorksheetFuncs[i][j]+1);
		  }
      }

      for (i=0; i<g_rgCommandFuncsRows; i++) 
	  {
		  for (j=0; j<g_rgCommandFuncsCols; j++) 
		  {
	    	  g_rgCommandFuncs[i][j][0] =
                (BYTE) lstrlen (g_rgCommandFuncs[i][j]+1);
		  }
      }

      for (i=0; i<g_rgMenuRows; i++) 
	  {
		  for (j=0; j<g_rgMenuCols; j++) 
		  {
	    	  g_rgMenu[i][j][0] = (BYTE) lstrlen (g_rgMenu[i][j]+1);
		  }
      }

      for (i=0; i<g_rgDialogRows; i++) 
	  {
		  for (j=0; j<g_rgDialogCols; j++)
		  {
	    	  g_rgDialog[i][j][0] = (BYTE) lstrlen (g_rgDialog[i][j]+1);
		  }
      }

      for (i=0; i<g_rgToolRows; i++) 
	  {
		  for (j=0; j<g_rgToolCols; j++)	
		  {   
	    	  g_rgTool[i][j][0] = (BYTE) lstrlen (g_rgTool[i][j]+1);
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


///***************************************************************************
// xlAutoOpen()
//
// Purpose: Microsoft Excel call this function when the DLL is loaded.
//
//  Microsoft Excel uses xlAutoOpen to load XLL files.
//  When you open an XLL file, the only action
//  Microsoft Excel takes is to call the xlAutoOpen function.
//
//  More specifically, xlAutoOpen is called:
//
//   - when you open this XLL file from the File menu,
//   - when this XLL is in the XLSTART directory, and is
//     automatically opened when Microsoft Excel starts,
//   - when Microsoft Excel opens this XLL for any other reason, or
//   - when a macro calls REGISTER(), with only one argument, which is the
//     name of this XLL.
//
//  xlAutoOpen is also called by the Add-in Manager when you add this XLL 
//  as an add-in. The Add-in Manager first calls xlAutoAdd, then calls
//  REGISTER("EXAMPLE.XLL"), which in turn calls xlAutoOpen.
//
//  xlAutoOpen should:
//
//   - register all the functions you want to make available while this
//     XLL is open,
//
//   - add any menus or menu items that this XLL supports,
//
//   - perform any other initialization you need, and
//
//   - return 1 if successful, or return 0 if your XLL cannot be opened.
//
// Parameters: None
//
// Returns: 1 if successful, or 0 if the XLL cannot be opened.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoOpen(void)
{

	static XLOPER xDLL,    // name of this DLL //
				  xMenu,   // xltypeMulti containing the menu //
				  xTool,   // xltypeMulti containing the toolbar //
				  xTest;   // used for menu test //
	LPXLOPER pxMenu;       // Points to first menu item //
	LPXLOPER px;           // Points to the current item //
	LPXLOPER pxTool;       // Points to first toolbar item //
	int i, j;              // Loop indices //
	HANDLE   hMenu;        // global memory holding menu //
	HANDLE  hTool;         // global memory holding toolbar //

  //
  // In the following block of code the name of the XLL is obtained by
  // calling xlGetName. This name is used as the first argument to the
  // REGISTER function to specify the name of the XLL. Next, the XLL loops
  // through the g_rgWorksheetFuncs[] table, and the g_rgCommandFuncs[]
  // tableregistering each function in the table using xlfRegister. 
  // Functions must be registered before you can add a menu item.
  //
	HMODULE hmodDLL1,hmodDLL2;
	char szFile1[256], szFile2[256];
	DWORD dwRet; 

	hmodDLL1 = GetModuleHandle("generic.dll");
	hmodDLL2 = GetModuleHandle(NULL);

	dwRet = GetModuleFileName(hmodDLL1, szFile1, 256);
	dwRet = GetModuleFileName(hmodDLL2, szFile2, 256);
	
	Excel(xlGetName, &xDLL, 0);

	for (i=0; i<g_rgWorksheetFuncsRows; i++) 
	{
		Excel(xlfRegister, 0, 1+ g_rgWorksheetFuncsCols,
			(LPXLOPER) &xDLL,
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][0]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][1]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][2]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][3]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][4]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][5]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][6]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][7]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][8]),
			(LPXLOPER) TempStr(g_rgWorksheetFuncs[i][9]));
	}

	for (i=0; i<g_rgCommandFuncsRows; i++) 
	{
		Excel(xlfRegister, 0, 1+ g_rgCommandFuncsCols,
			(LPXLOPER) &xDLL,
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][0]),
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][1]),
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][2]),
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][3]),
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][4]),
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][5]),
			(LPXLOPER) TempStr(g_rgCommandFuncs[i][6]));
	}

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
  // This block uses TempStr() and TempNum(). Both create a temporary
  // XLOPER. The XLOPER created by TempStr() contains the string passed to
  // it. The XLOPER created by TempNum() contains the number passed to it.
  // The Excel() function frees the allocated temporary memory. Both
  // functions are part of the framework library.
  //

	Excel(xlfGetBar, &xTest, 3, TempInt(10), TempStr(" Generic"), TempInt(0));

	if (xTest.xltype == xltypeErr) 
	{
		hMenu = GlobalAlloc(GMEM_MOVEABLE,sizeof(XLOPER)*5*g_rgMenuRows);
			px = pxMenu = (LPXLOPER) GlobalLock(hMenu);

		for (i=0; i<g_rgMenuRows; i++) 
		{
			for (j=0; j<5; j++) 
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

		Excel(xlfAddMenu,0,3,TempNum(10),(LPXLOPER)&xMenu,TempStr(" Help"));

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
  // This block uses TempInt(), TempBool(), and TempMissing(). All three
  // create a temporary XLOPER. The XLOPER created by TempInt() contains
  // the integer passed to it. TempBool() creates an XLOPER containing the
  // boolean value passed to it. TempMissing() creates an XLOPER that
  // simulates a missing argument. The Excel() function frees the temporary
  // memory associated with these functions. All three are part of the
  // framework library.
  //

	Excel(xlfGetToolbar, &xTest, 2, TempInt(1), TempStr(" Test"));

	if (xTest.xltype == xltypeErr) 
	{
		hTool = GlobalAlloc(GMEM_MOVEABLE, sizeof(XLOPER)*8*g_rgToolRows);
		px = pxTool = (LPXLOPER) GlobalLock(hTool);

		for (i = 0; i < g_rgToolRows; i++) 
		{
			for (j = 0; j < 8; j++) 
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

		Excel(xlfAddToolbar,0,2,TempStr(" Test"),(LPXLOPER)&xTool);

		Excel(xlcShowToolbar, 0, 6, TempStr(" Test"), TempBool(1),
	      TempInt(5), TempMissing(), TempMissing(), TempInt(999));

		GlobalUnlock(hTool);
		GlobalFree(hTool);
	}

  // Free the XLL filename //
	Excel(xlFree, 0, 2, (LPXLOPER) &xTest, (LPXLOPER) &xDLL);

	return 1;
}


///***************************************************************************
// xlAutoClose()
//
// Purpose: Microsoft Excel call this function when the DLL is unloaded.
//
//  xlAutoClose is called by Microsoft Excel:
//
//   - when you quit Microsoft Excel, or
//   - when a macro sheet calls UNREGISTER(), giving a string argument
//     which is the name of this XLL.
//
//  xlAutoClose is called by the Add-in Manager when you remove this XLL from
//  the list of loaded add-ins. The Add-in Manager first calls xlAutoRemove,
//  then calls UNREGISTER("EXAMPLE.XLL"), which in turn calls xlAutoClose.
// 
//  xlAutoClose is called by EXAMPLE.XLL by the function fExit. This function
//  is called when you exit Example.
// 
//  xlAutoClose should:
// 
//   - Remove any menus or menu items that were added in xlAutoOpen,
// 
//   - do any necessary global cleanup, and
// 
//   - delete any names that were added (names of exported functions, and 
//     so on). Remember that registering functions may cause names to 
//     be created.
// 
//  xlAutoClose does NOT have to unregister the functions that were registered
//  in xlAutoOpen. This is done automatically by Microsoft Excel after
//  xlAutoClose returns.
// 
//  xlAutoClose should return 1.
//
// Parameters: None
//
// Returns: 1
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoClose(void)
{
    int i;
    XLOPER xRes;

    //
    // This block first deletes all names added by xlAutoOpen or
    // xlAutoRegister. Next, it checks if the drop-down menu Generic still
    // exists. If it does, it is deleted using xlfDeleteMenu. It then checks
    // if the Test toolbar still exists. If it is, xlfDeleteToolbar is
    // used to delete it.
    //

    //
    // Due to a bug in Excel 5 for NT and Excel 7
    // the following code to delete the defined names
    // does not work.  There is no way to delete these
    // names once they are Registered
    // The code is left in, in hopes that it will be
    // fixed in a future version.
    //
    for (i = 0; i < g_rgWorksheetFuncsRows; i++)
    	Excel(xlfSetName, 0, 1, TempStr(g_rgWorksheetFuncs[i][2]));
    
    for (i = 0; i < g_rgCommandFuncsRows; i++)
    	Excel(xlfSetName, 0, 1, TempStr(g_rgCommandFuncs[i][2]));
    //
    // Everything else works as documented
    //
    Excel(xlfGetBar, &xRes, 3, TempInt(10), TempStr(" Generic"), TempInt(0));

    if (xRes.xltype != xltypeErr) 
	{
      Excel(xlfDeleteMenu, 0, 2, TempNum(10), TempStr(" Generic"));

      // Free the XLOPER returned by xlfGetBar //
      Excel(xlFree, 0, 1, (LPXLOPER) &xRes);
    }

    Excel(xlfGetToolbar, &xRes, 2, TempInt(7), TempStr(" Test"));

    if (xRes.xltype != xltypeErr) 
	{
		  Excel(xlfDeleteToolbar, 0, 1, TempStr(" Test"));
    
		  // Free the XLOPER returned by xlfGetToolbar //
		  Excel(xlFree, 0, 1, (LPXLOPER) &xRes);
    }

    return 1;
}


///***************************************************************************
// lpstricmp()
//
// Purpose: Compares two Pascal strings and checks if they are equal,
//          without case sensitivity. This routine is useful for
//          xlAutoRegister         
//
// Parameters:
//
//      LPSTR s     First string
//      LPSTR t     Second string
//
// Returns: 
//
//      int         0 if they are equal
//                  Nonzero otherwise
//
// Comments:
//
//      Unlike the usual string functions, lpstricmp
//      doesn't care about collating sequence.
//
// History:  Date       Author        Reason
///***************************************************************************

int lpstricmp(LPSTR s, LPSTR t)
{
  int i;

  if (*s != *t)
    return 1;

  for (i = 1; i <= s[0]; i++) 
  {
	  if (tolower(s[i]) != tolower(t[i]))
		    return 1;
  }
  return 0;
}


///***************************************************************************
// xlAutoRegister()
//
// Purpose:
//
//   This function is called by Microsoft Excel if a macro sheet tries to
//   register a function without specifying the type_text argument. If that
//   happens, Microsoft Excel calls xlAutoRegister, passing the name of the
//   function that the user tried to register. xlAutoRegister should use the
//   normal REGISTER function to register the function, only this time it must
//   specify the type_text argument. If xlAutoRegister does not recognize the
//   function name, it should return a #VALUE! error. Otherwise, it should
//   return whatever REGISTER returned.
//
// Parameters:
//
//      LPXLOPER pxName     xltypeStr containing the
//                          name of the function
//                          to be registered. This is not
//                          case sensitive.
//
// Returns: 
//
//      LPXLOPER            xltypeNum containing the result
//                          of registering the function,
//                          or xltypeErr containing #VALUE!
//                          if the function could not be
//                          registered.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER WINAPI xlAutoRegister(LPXLOPER pxName)
{
  static XLOPER xDLL, xRegId;
  int i;

  //
  // This block initializes xRegId to a #VALUE! error first. This is done in
  // case a function is not found to register. Next, the code loops through 
  // the functions in g_rgFuncs[] and uses lpstricmp to determine if the 
  // current row in g_rgFuncs[] represents the function that needs to be 
  // registered. When it finds the proper row, the function is registered 
  // and the register ID is returned to Microsoft Excel. If no matching 
  // function is found, an xRegId is returned containing a #VALUE! error.
  //

  xRegId.xltype = xltypeErr;
  xRegId.val.err = xlerrValue;


  for (i=0; i<g_rgWorksheetFuncsRows; i++) 
  {
    if (!lpstricmp(g_rgWorksheetFuncs[i][0], pxName->val.str))	
	{
		Excel(xlfRegister, 0, 1+ g_rgWorksheetFuncsCols,
	          (LPXLOPER) &xDLL,
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][0]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][1]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][2]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][3]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][4]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][5]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][6]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][7]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][8]),
	          (LPXLOPER) TempStr(g_rgWorksheetFuncs[i][9]));
 	    /// Free oper returned by xl //
	    Excel(xlFree, 0, 1, (LPXLOPER) &xDLL);

	    return (LPXLOPER) &xRegId;
	}
  }

  for (i=0; i<g_rgCommandFuncsRows; i++) 
  {
	  if (!lpstricmp(g_rgCommandFuncs[i][0], pxName->val.str))	
	  {
		Excel(xlfRegister, 0, 1+ g_rgCommandFuncsCols,
	          (LPXLOPER) &xDLL,
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][0]),
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][1]),
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][2]),
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][3]),
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][4]),
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][5]),
	          (LPXLOPER) TempStr(g_rgCommandFuncs[i][6]));
 	    /// Free oper returned by xl //
	    Excel(xlFree, 0, 1, (LPXLOPER) &xDLL);

	    return (LPXLOPER) &xRegId;
	  }
  }     
  return (LPXLOPER) &xRegId;
}

///***************************************************************************
// xlAutoAdd()
//
// Purpose:
//
//    This function is called by the Add-in Manager only. When you add a
//    DLL to the list of active add-ins, the Add-in Manager calls xlAutoAdd()
//    and then opens the XLL, which in turn calls xlAutoOpen.
//
// Parameters:
//
// Returns: 
//
//      int       1
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoAdd(void)
{                                 
	char szBuf[255];
	
	wsprintf((LPSTR)szBuf, " Thank you for adding GENERIC.XLL\n "
    "built on %s at %s", __DATE__, __TIME__);
	
  // Display a dialog box indicating that the XLL was successfully added //
  Excel(xlcAlert, 0, 2, TempStr(szBuf), TempInt(2));
  return 1;
}

///***************************************************************************
// xlAutoRemove()
//
// Purpose:
//
//    This function is called by the Add-in Manager only. When you remove
//    an XLL from the list of active add-ins, the Add-in Manager calls
//    xlAutoRemove() and then UNREGISTER("EXAMPLE.XLL").
//   
//    You can use this function to perform any special tasks that need to be
//    performed when you remove the XLL from the Add-in Manager's list
//    of active add-ins. For example, you may want to delete an
//    initialization file when the XLL is removed from the list.
//
// Parameters:
//
// Returns: 
//
//      int       1
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI xlAutoRemove(void)
{
    // Show a dialog box indicating that the XLL was successfully removed //
    Excel(xlcAlert, 0, 2, TempStr(" Thank you for removing GENERIC.XLL!"),
    	TempInt(2));
    return 1;
}

///***************************************************************************
// xlAddInManagerInfo()
//
// Purpose:
//
//    This function is called by the Add-in Manager to find the long name
//    of the add-in. If xAction = 1, this function should return a string
//    containing the long name of this XLL, which the Add-in Manager will use
//    to describe this XLL. If xAction = 2 or 3, this function should return
//    #VALUE!.
//
// Parameters:
//
//      LPXLOPER xAction    What information you want. One of:
//                          1 = the long name of the
//                              add-in
//                          2 = reserved
//                          3 = reserved
//
// Returns: 
//
//      LPXLOPER            The long name or #VALUE!.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER WINAPI xlAddInManagerInfo(LPXLOPER xAction)
{
    static XLOPER xInfo, xIntAction;

    //
    // This code coerces the passed-in value to an integer. This is how the
    // code determines what is being requested. If it receives a 1, 
    // it returns a string representing the long name. If it receives 
    // anything else, it returns a #VALUE! error.
    //

    Excel(xlCoerce, &xIntAction, 2, xAction, TempInt(xltypeInt));

    if(xIntAction.val.w == 1) 
	{
		xInfo.xltype = xltypeStr;
		xInfo.val.str = "\026Generic Standalone DLL";
    }
    else 
	{
		xInfo.xltype = xltypeErr;
		xInfo.val.err = xlerrValue;
    }

    return (LPXLOPER) &xInfo;
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
//      BOOL                TRUE if message processed, FALSE if not.
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

BOOL CALLBACK DIALOGMsgProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM lParam)
{

 //
 // This block is a very simple message loop for a dialog box. It checks for
 // only three messages. When it receives WM_INITDIALOG, it uses a buffer to
 // set a static text item to the amount of free space on the stack. When it
 // receives WM_CLOSE it posts a CANCEL message to the dialog box. When it
 // receives WM_COMMAND it checks if the OK button was pressed. If it was,
 // the dialog box is closed and control returned to fShowDialog.
 //

  switch(message) 
  {

	case WM_INITDIALOG:
	  SetDlgItemText(hWndDlg, FREE_SPACE, (LPSTR)g_szBuffer);
	  break;
        
	case WM_CLOSE:
	  PostMessage(hWndDlg, WM_COMMAND, IDCANCEL, 0L);
	  break;

	case WM_COMMAND:
	  switch(wParam) 
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
//  When a modal dialog box is displayed over Microsoft Excel's window, the
//  cursor is a busy cursor over Microsoft Excel's window. This WndProc traps
//  WM_SETCURSORs and changes the cursor back to a normal arrow.
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
//      LRESULT            0 if message handled, otherwise the result of the
//                      default WndProc
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

// Create a place to store Microsoft Excel's WndProc address //
static FARPROC g_lpfnExcelWndProc = NULL;

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
//      HANDLE hWndExcel        This is a handle to Microsoft Excel's hWnd
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
  // use of GetWindowLong(). It stores this value in a global that can be
  // used to call the default WndProc and also to restore it. Finally, it
  // replaces this address with the address of ExcelCursorProc using
  // SetWindowLong().
  //

  g_lpfnExcelWndProc = (FARPROC) GetWindowLong(hWndExcel, GWL_WNDPROC);
  SetWindowLong(hWndExcel, GWL_WNDPROC, (LONG)(FARPROC) ExcelCursorProc);
}

///***************************************************************************
// UnhookExcelWindow()
//
// Purpose:
//
//   This is the function that removes the ExcelCursorProc that was
//   called before Microsoft Excel's main WndProc.
//
// Parameters:
//
//      HANDLE hWndExcel        This is a handle to Microsoft Excel's hWnd
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
  // SetWindowLong to restore the address that was saved into
  // g_lpfnExcelWndProc by HookExcelWindow(). It then sets g_lpfnExcelWndProc
  // to NULL.
  //

  SetWindowLong(hWndExcel, GWL_WNDPROC, (LONG) g_lpfnExcelWndProc);
  g_lpfnExcelWndProc = NULL;
}

///***************************************************************************
// fShowDialog()
//
// Purpose:
//          This function loads and shows the native Windows dialog box.
//
// Parameters:
//
// Returns: 
//
//      int 0                   To indicate successful completion
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) int WINAPI fShowDialog(void)
{
  static  XLOPER xRes, xStr;
  int     i, nRc = 0;              // return code //

  //
  // This block first gets the amount of space left on the stack, coerces it
  // to a string, and then copies it into a buffer. It obtains this
  // information so that it can be displayed in the dialog box. The coerced
  // string is then released with xlFree. The buffer is then converted into
  // a null-terminated string. The hWnd of Microsoft Excel is then obtained
  // using GetHwnd. Messages are enabled using xlEnableXLMsgs in
  // preparation for showing the dialog box. hWndMain is passed to
  // HookExcelWindow so that an arrow cursor is displayed over Microsoft
  // Excel's window. The dialog box is then displayed. The hWnd of Microsoft
  // Excel is then obtained. It is obtained a second time in case there
  // are multiple instances of Microsoft Excel using this XLL. hWndMain may
  // have changed while the dialog box was displayed. If the wrong hWnd is
  // passed to UnhookExcelWindow(), Microsoft Excel will most likely crash.
  // The obtained hWnd is then passed to UnhookExcelWindow(). Messages are
  // then disabled using xlDisableXLMsgs.
  //

  Excel(xlStack, (LPXLOPER) &xRes, 0);

  Excel(xlCoerce, (LPXLOPER) &xStr, 2, (LPXLOPER) &xRes, (LPXLOPER) TempInt(xltypeStr));

  lstrcpy(g_szBuffer, xStr.val.str);

  Excel(xlFree, 0, 1, (LPXLOPER) &xStr);

  nRc = g_szBuffer[0];
  for (i = 0; i <= nRc; i++)
    g_szBuffer[i] = g_szBuffer[i+1];
  g_szBuffer[nRc] = '\0';

  GetHwnd(&g_hWndMain);

  Excel(xlEnableXLMsgs, 0, 0);
  HookExcelWindow(g_hWndMain);

  nRc = DialogBox(g_hInst, "TESTSTACK", g_hWndMain, DIALOGMsgProc);

  UnhookExcelWindow(g_hWndMain);
  Excel(xlDisableXLMsgs, 0, 0);

  return(0);
}


///***************************************************************************
// EnumProc()
//
// Purpose:
//
//      This function compares the class name of Excel's main window to
//      the window currently being enumerated.  If it matches it then checks
//      to see if the Low Word of the hWnd matches the value Excel returned
//      using xlGetHwnd.  If it matches it returns the Excel hWnd, otherwise
//      it contiues the enumeration. 
//
// Parameters:
//
//      HWND hwnd         handle to parent window
//      LPARAM lParam     pointer to an EnumStruct
//
// Returns: 
//
//      BOOL    FALSE  Ends enumeration
//              TRUE   Continues enumeration
//
// Comments:
//
// History:  Date       Author        Reason
//           5/31/96     AS            xlGetHwnd only returns 16-bit value
///***************************************************************************

//
// Structure to pass the the EnumWindowsProc
//
typedef struct _EnumStruct {
	HWND		hwnd;   // Return value for Excel's main hWnd
	unsigned short	wLoword; //Contains LowWord of Excel's main hWnd
} EnumStruct;

#define CLASS_NAME_BUFFER	50

BOOL CALLBACK EnumProc(HWND hwnd, EnumStruct * pEnum)
{
	// first check the class of the window.  Must be "XLMAIN"
	char rgsz[CLASS_NAME_BUFFER];
	
  GetClassName(hwnd, rgsz, CLASS_NAME_BUFFER);
	
  if (!lstrcmpi(rgsz, "XLMAIN")) 
  {
		// if that hits, check the loword of the window handle
		if (LOWORD((DWORD) hwnd) == pEnum->wLoword) 
		{
			// We have a match, return Excel's main hWnd
			pEnum->hwnd = hwnd;
			return FALSE;
		}
	}

	// no luck - continue the enumeration
	return TRUE;
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
//      HWND * phwnd      Will contain Excel's hWnd
//
// Returns: 
//
//      BOOL    FALSE  Could not find Excel's hWnd
//              TRUE   Found Excel's hWnd
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************
BOOL GetHwnd(HWND * pHwnd)
{
	XLOPER x;
  
  //
  // xlGetHwnd only returns the LoWord of Excel's hWnd
  // so all the windows have to be enumerated to see
  // which match the LoWord retuned by xlGetHwnd
  //
  if (Excel4(xlGetHwnd, &x, 0) == xlretSuccess) 
  {
		EnumStruct enm;
		
		enm.hwnd = NULL;
		enm.wLoword = x.val.w;

		EnumWindows((WNDENUMPROC) EnumProc, (LPARAM) &enm);
		
		if (enm.hwnd != NULL) 
		{
			*pHwnd = enm.hwnd;
			return TRUE;
		}
  }
  return FALSE;
}


///***************************************************************************
// Func1()
//
// Purpose:
//
//      This is a typical user-defined function provided by an XLL.
//
// Parameters:
//
//      LPXLOPER x      (Ignored)
//
// Returns: 
//
//      LPXLOPER        Always the string "Func1"
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER WINAPI Func1 (LPXLOPER x)
{
  static XLOPER xResult;

  //
  // This function demonstrates returning a string value. The return
  // type is set to a string and filled with the name of the function.
  //

  xResult.xltype = xltypeStr;
  xResult.val.str = "\005Func1";

  return (LPXLOPER) &xResult;
}

///***************************************************************************
// FuncSum()
//
// Purpose:
//
//     This is a typical user-defined function provided by an XLL. This
//     function takes 1-29 arguments and computes their sum. Each argument
//     can be a single number, a range, or an array.
//
// Parameters:
//
//      LPXLOPER ...    1 to 29 arguments
//                      (can be references or values)
//
// Returns: 
//
//      LPXLOPER        The sum of the arguments
//                      or #VALUE! if there are
//                      nonnumerics in the
//
// Comments:
//
// History:  Date       Author        Reason
///***************************************************************************

__declspec(dllexport) LPXLOPER WINAPI FuncSum(
          LPXLOPER px1,LPXLOPER px2,LPXLOPER px3,LPXLOPER px4,
          LPXLOPER px5,LPXLOPER px6,LPXLOPER px7,LPXLOPER px8,
          LPXLOPER px9,LPXLOPER px10,LPXLOPER px11,LPXLOPER px12,
          LPXLOPER px13,LPXLOPER px14,LPXLOPER px15,LPXLOPER px16,
          LPXLOPER px17,LPXLOPER px18,LPXLOPER px19,LPXLOPER px20,
          LPXLOPER px21,LPXLOPER px22,LPXLOPER px23,LPXLOPER px24,
          LPXLOPER px25,LPXLOPER px26,LPXLOPER px27,LPXLOPER px28,
          LPXLOPER px29)
{
  static XLOPER xResult;  // Return value 
  double d=0;             // Accumulate result 
  int iArg;               // The argument being processed 
  LPXLOPER FAR *ppxArg;   // Pointer to the argument being processed 
  XLOPER xMulti;          // Argument coerced to xltypeMulti 
  WORD i;                 // Row and column counters for arrays 
  LPXLOPER px;            // Pointer into array 
  int error = -1;         // -1 if no error; error code otherwise 

  //
  // This block accumulates the arguments passed in. Because FuncSum is a
  // Pascal function, the arguments will be evaluated right to left. For
  // each argument, this code checks the type of the argument and then does
  // things necessary to accumulate that argument type. Unless the caller
  // actually specified all 29 arguments, there will be some missing
  // arguments. The first case handles this. The second case handles
  // arguments that are numbers. This case just adds the number to the
  // accumulator. The third case handles references or arrays. It coerces
  // references to an array. It then loops through the array, switching on
  // XLOPER types and adding each number to the accumulator. The fourth
  // case handles being passed an error. If this happens, the error is
  // stored in error. The fifth and default case handles being passed
  // something odd, in which case an error is set. Finally, a check is made
  // to see if error was ever set. If so, an error of the same type is
  // returned; if not, the accumulator value is returned.
  //

  for (iArg = 0; iArg < 29; iArg++) 
  {
    ppxArg = &px1 + iArg;

	  switch ((*ppxArg)->xltype) 
	  {
	  
		case xltypeMissing:
			break;

		case xltypeNum:
			d += (*ppxArg)->val.num;
			break;

		case xltypeRef:
		case xltypeSRef:
		case xltypeMulti:
			if (xlretUncalced == Excel(xlCoerce, &xMulti, 2,
	          (LPXLOPER) *ppxArg, TempInt(xltypeMulti))) 
			{
			//
			// That coerce might have failed due to an 
			// uncalced cell, in which case, we need to 
			// return immediately. Microsoft Excel will
			// call us again in a moment after that cell
			// has been calced.
			//
			   return 0;
			}

			for (i = 0;
				i < (xMulti.val.array.rows * xMulti.val.array.columns);
				i++) 
			{

				// obtain a pointer to the current item //
				px = xMulti.val.array.lparray + i;

				// switch on XLOPER type //
				switch (px->xltype) 
				{

					// if a num accumulate it //
					case xltypeNum:
						d += px->val.num;
						break;

					// if an error store in error //
					case xltypeErr:
						error = px->val.err;
						break;

					// if missing do nothing //
					case xltypeNil:
						break;

					// if anything else set error //
					default:
						error = xlerrValue;
						break;
				}
			}

			// free the returned array //
			Excel(xlFree, 0, 1, (LPXLOPER) &xMulti);
			break;

		case xltypeErr:
			error = (*ppxArg)->val.err;
	        break;

		default:
			error = xlerrValue;
			break;
	  }

  }
  
  if (error != -1) 
  {
	  xResult.xltype = xltypeErr;
	  xResult.val.err = error;
  }
  else 
  {
	  xResult.xltype = xltypeNum;
	  xResult.val.num = d;
  }

  return (LPXLOPER) &xResult;
}


///***************************************************************************
// fDance()
//
// Purpose:
//
//    This is an example of a lengthy operation. It calls the function xlAbort
//    occasionally. This yields the processor (allowing cooperative 
//    multitasking), and checks if the user has pressed ESC to abort this
//    operation. If so, it offers the user a chance to cancel the abort.
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

__declspec(dllexport) int WINAPI fDance(void)
{
  DWORD dtickStart;
  XLOPER xAbort, xConfirm;
  int boolSheet;
  int col=0;
  char rgch[32];

  //
  // Check what kind of sheet is active. If it is a worksheet or macro
  // sheet, this function will move the selection in a loop to show
  // activity. In any case, it will update the status bar with a countdown.
  //
  // Call xlSheetId; if that fails the current sheet is not a macro sheet or
  // worksheet. Next, get the time at which to start. Then start a while
  // loop that will run for one minute. During the while loop, check if the
  // user has pressed ESC. If true, confirm the abort. If the abort is
  // confirmed, clear the message bar and return; if the abort is not
  // confirmed, clear the abort state and continue. After checking for an
  // abort, move the active cell if on a worksheet or macro. Then
  // update the status bar with the time remaining.
  //
  // This block uses TempActiveCell(), which creates a temporary XLOPER.
  // The XLOPER contains a reference to a single cell on the active sheet.
  // This function is part of the framework library.
  //

  boolSheet = (Excel4(xlSheetId, 0, 0) == xlretSuccess);

  dtickStart = GetTickCount();

  while (GetTickCount() < dtickStart + 60000L) 
  {
    Excel(xlAbort, &xAbort, 0);
    if (xAbort.val.bool) 
	{
		Excel(xlcAlert, &xConfirm, 2,
	          TempStr(" Are you sure you want to cancel this operation?"),
	          TempNum(1));

	    if (xConfirm.val.bool) 
		{
		    Excel(xlcMessage, 0, 1, TempBool(0));
		    return 1;
	    }
	    else 
		{
		    Excel(xlAbort, 0, 1, TempBool(0));
	    }
    }

    if (boolSheet) 
	{
	    Excel(xlcSelect, 0, 1, TempActiveCell(0,(BYTE)col));
	    col = (col + 1) & 3;
    }

    wsprintf(rgch," 0:%lu",
	      (60000 + dtickStart - GetTickCount()) / 1000L);

    Excel(xlcMessage, 0, 2, TempBool(1), TempStr(rgch));

  }

  Excel(xlcMessage, 0, 1, TempBool(0));

  return 1;
}

///***************************************************************************
// fDialog()
//
// Purpose:
//
//         An example of how to create a Microsoft Excel
//         UDD (User Defined Dialog) from a DLL.
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

__declspec(dllexport) int WINAPI fDialog(void)
{
  int i, j;
  HANDLE hrgrgx;
  LPXLOPER prgrgx, px;
  XLOPER xDialog;
  XLOPER rgxList[5];
  XLOPER xList;
  XLOPER xDialogReturned;

  //
  // This block first allocates memory to hold the dialog box array. It then
  // fills this array with information from g_rgDialog[]. It replaces any
  // empty entries with NIL XLOPERs while filling the array. It then
  // creates the name "GENERIC_List1" to refer to an array containing the
  // list box values. The dialog box is then displayed. The dialog box is
  // redisplayed using the results from the last dialog box. Then the arrays
  // are freed and the name "GENERIC_List1" is deleted.
  //

  px = prgrgx = (LPXLOPER) GlobalLock(hrgrgx = 
    GlobalAlloc(GMEM_MOVEABLE, sizeof(XLOPER)* 7 * g_rgDialogRows));

  for (i = 0; i < g_rgDialogRows; i++) 
  {
	for (j = 0; j < g_rgDialogCols; j++)	
	{
	  if (g_rgDialog[i][j][0] == 0)
	    px->xltype = xltypeNil;
	  else 
	  {
		  px->xltype = xltypeStr;
		  px->val.str = g_rgDialog[i][j];
	  }
	  px++;
	}
  }

  xDialog.xltype = xltypeMulti;
  xDialog.val.array.lparray = prgrgx;
  xDialog.val.array.rows = g_rgDialogRows;
  xDialog.val.array.columns = 7;

  rgxList[0].val.str = (LPSTR) " Bake";
  rgxList[1].val.str = (LPSTR) " Broil";
  rgxList[2].val.str = (LPSTR) " Sizzle";
  rgxList[3].val.str = (LPSTR) " Fry";
  rgxList[4].val.str = (LPSTR) " Saute";

  for (i = 0; i < 5; i++) 
  {
    rgxList[i].xltype = xltypeStr;
    rgxList[i].val.str[0]=(BYTE)lstrlen(rgxList[i].val.str + 1);
  }

  xList.xltype = xltypeMulti;
  xList.val.array.lparray = (LPXLOPER) & rgxList;
  xList.val.array.rows = 5;
  xList.val.array.columns = 1;

  Excel(xlfSetName, 0, 2, TempStr(" GENERIC_List1"), (LPXLOPER) &xList);

  Excel(xlfDialogBox, &xDialogReturned, 1, (LPXLOPER) &xDialog);

  Excel(xlfDialogBox, 0, 1, (LPXLOPER) &xDialogReturned);

  Excel(xlFree, 0, 1, (LPXLOPER) &xDialogReturned);

  GlobalUnlock(hrgrgx);
  GlobalFree(hrgrgx);

  Excel(xlfSetName, 0, 1, TempStr(" GENERIC_List1"));

  return 1;
}


///***************************************************************************
// fExit()
//
// Purpose:
//
//  This is a user-initiated routine to exit GENERIC.XLL You may be tempted to
//  simply call UNREGISTER("GENERIC.XLL") in this function. Don't do it! It
//  will have the effect of forcefully unregistering all of the functions in
//  this DLL, even if they are registered somewhere else! Instead, unregister
//  the functions one at a time.
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

__declspec(dllexport) int WINAPI fExit(void)
{
  XLOPER  xDLL,                   // The name of this DLL //
          xFunc,                  // The name of the function //
          xRegId;                 // The registration ID //
  int     i;

  //
  // This code gets the DLL name. It then uses this along with information
  // from g_rgFuncs[] to obtain a REGISTER.ID() for each function. The
  // register ID is then used to unregister each function. Then the code
  // frees the DLL name and calls xlAutoClose.
  //

  // Make xFunc a string //
  xFunc.xltype = xltypeStr;

  Excel(xlGetName, &xDLL, 0);

  for (i = 0; i < g_rgWorksheetFuncsRows; i++) 
  {
	  xFunc.val.str = (LPSTR) (g_rgWorksheetFuncs[i][0]);
	  Excel(xlfRegisterId,&xRegId,2,(LPXLOPER)&xDLL,(LPXLOPER)&xFunc);
	  Excel(xlfUnregister, 0, 1, (LPXLOPER) &xRegId);
  }

  for (i = 0; i < g_rgCommandFuncsRows; i++) 
  {
	  xFunc.val.str = (LPSTR) (g_rgCommandFuncs[i][0]);
	  Excel(xlfRegisterId,&xRegId,2,(LPXLOPER)&xDLL,(LPXLOPER)&xFunc);
	  Excel(xlfUnregister, 0, 1, (LPXLOPER) &xRegId);
  }

  Excel(xlFree, 0, 1,  (LPXLOPER) &xDLL);

  return xlAutoClose();
}
