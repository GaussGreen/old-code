#ifndef XLFUNCTION_H
#define XLFUNCTION_H

#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <xlcall.h>
#include <framewrk.h>
#ifdef __cplusplus
}
#endif

#include <CCString.h>
#ifdef STL_WIN32
#pragma warning(disable:4786)
#include <vector>
using namespace std;
#else	// STL_WIN32
#include <vector.h>
#endif	// STL_WIN32

#include "XLArgument.h"

class XLFunction
{
public:
	XLFunction () {};
	XLFunction (const char* n_name, int n_type, const char* n_helptext, const char* n_result_helptext, int n_result_type);
	~XLFunction () {};

	virtual LPXLOPER done (vector<LPXLOPER>& XL_arguments) { return NULL; };

	inline void addArgument (XLArgument& n_argument)
	{
		arguments.push_back (n_argument);
	}
	inline void clear ()
	{
		arguments.clear ();
	}
	inline int getNbArguments () const
	{
		return arguments.size ();
	}

	inline XLArgument getArgument (int index) const
	{
		return arguments[index];
	}

	inline CCString getName () const
	{
		return name;
	}

	inline int getType () const
	{
		return type;
	}

	inline CCString getHelptext () const
	{
		return helptext;
	}

	inline CCString getResult_helptext () const
	{
		return result_helptext;
	}

	inline int getResult_type () const
	{
		return result_type;
	}

private:
	CCString name;
	int type;
	CCString helptext;
	vector<XLArgument> arguments;
	CCString result_helptext;
	int result_type;
};

#endif	// XLFUNCTION_H