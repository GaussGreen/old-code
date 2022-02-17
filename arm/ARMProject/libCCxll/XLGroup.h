#ifndef XLGROUP_H
#define XLGROUP_H

#include <stdio.h>

#include <CCString.h>
#ifdef STL_WIN32
#pragma warning(disable:4786)
#include <vector>
using namespace std;
#else	// STL_WIN32
#include <vector.h>
#endif	// STL_WIN32

#include "CCxll.h"
#include "XLFunction.h"

#define MAX_NB_OF_ARGUMENTS		25
#define NB_FIRST_ROWS			9

class XLGroup
{
public:
	XLGroup () {};
	XLGroup (const char* n_name);
	~XLGroup () {};

	void addFunction (XLFunction& n_function)
	{
		functions.push_back (n_function);
	}
	void clear ()
	{
		functions.clear ();
	}
	int getNbFunctions () const
	{
		return functions.size ();
	}

	operator LPSTR**()
	{
		LPSTR** g_rgWorksheetFuncs = new LPSTR*[getNbFunctions ()];
		for(int i = 0; i < getNbFunctions (); i++)
		{
			g_rgWorksheetFuncs[i] = new LPSTR[MAX_NB_OF_ARGUMENTS + NB_FIRST_ROWS];
			
			g_rgWorksheetFuncs[i][0] = XL_StrC2StrPascal (functions[i].getName ());

			CCString codes ("R");
			for(int u = 0; u < functions[i].getNbArguments (); u++)
			{
				codes += "R";
			}
			g_rgWorksheetFuncs[i][1] = XL_StrC2StrPascal (codes);

			g_rgWorksheetFuncs[i][2] = XL_StrC2StrPascal (functions[i].getName ());

			CCString list ("");
			int u;
			for(u = 0; u < functions[i].getNbArguments () - 1; u++)
			{
				list += functions[i].getArgument (u).getName ();
				list += ",";
			}
			list += functions[i].getArgument (functions[i].getNbArguments () - 1).getName ();
			g_rgWorksheetFuncs[i][3] = XL_StrC2StrPascal (list);

			char type_str[2];
			//sprintf (type_str, "%d", functions[i].getType ());
			g_rgWorksheetFuncs[i][4] = XL_StrC2StrPascal (CCString (type_str));
						
			g_rgWorksheetFuncs[i][5] = XL_StrC2StrPascal (name);
			
			g_rgWorksheetFuncs[i][6] = XL_StrC2StrPascal (CCString (""));
			g_rgWorksheetFuncs[i][7] = XL_StrC2StrPascal (CCString (""));

			g_rgWorksheetFuncs[i][8] = XL_StrC2StrPascal (functions[i].getHelptext ());
		    int j;					
			for(j = 0; j < functions[i].getNbArguments (); j++)
			{
				g_rgWorksheetFuncs[i][j + NB_FIRST_ROWS] = XL_StrC2StrPascal (functions[i].getArgument (u).getHelptext ());
			}
			g_rgWorksheetFuncs[i][j + NB_FIRST_ROWS] = XL_StrC2StrPascal (functions[i].getResult_helptext ());
			for(j = functions[i].getNbArguments () + 1; j < MAX_NB_OF_ARGUMENTS + 1; j++)
			{
				g_rgWorksheetFuncs[i][j + NB_FIRST_ROWS] = XL_StrC2StrPascal (CCString (""));
			}
		}
		return g_rgWorksheetFuncs;
	}
	
private:
	CCString name;
	vector<XLFunction> functions;
};

#ifdef XLGroup_cpp

#endif	// XLGroup_cpp

#endif	// XLGROUP_H

// EOF %M%