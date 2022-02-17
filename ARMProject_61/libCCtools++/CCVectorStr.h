/* -------------------------------------------------------------------------
 
   File: CCVectorStr.h
   Path: /home/julia/projets/dev/Common/libCCtools++/SCCS/s.CCVectorStr.h
   Description: header CCVectorStr
   Created: 99/07/27
   Author: Charles-Emmanuel MUSY
   Modified: 99/07/27 08:56:47
   Last maintained by: Charles-Emmanuel MUSY
   Revision: 1.2
 
   -------------------------------------------------------------------------
 
   Note:
 
   ------------------------------------------------------------------------- */
 
#ifndef CCVECTORSTR_H
#define CCVECTORSTR_H

#include <CCcommon.h>
SCCS_ID (CCVectorStr_h_SccsId, "@(#)CCVectorStr.h	1.2, modified 99/07/27");

#ifdef STL_WIN32
#include <vector>
#else
#include <vector.h>
#endif	// STL_WIN32

#include "CCString.h"

class CCVectorStr
{
public:
	CCVectorStr () : values ()
	{
	}

	CCVectorStr (int n_size) : values (n_size)
	{
	}

#ifdef STL_WIN32
	CCVectorStr (std::vector<CCString>& n_values) : values ()
#else
	CCVectorStr (vector<CCString>& n_values) : values ()
#endif	// STL_WIN32
	{
		set (n_values);
	}

#ifdef STL_WIN32
	inline void set (std::vector<CCString>& n_values)
#else
	inline void set (vector<CCString>& n_values)
#endif	// STL_WIN32
	{
		clear ();
		for(int i = 0; i < n_values.size (); i++)
		{
			values.push_back (n_values[i]);
		}
	}	

	inline CCString operator [](int index)
	{
		return values[index];
	}

	inline void push_back (CCString& value)
	{
		values.push_back (value);
	}

	inline void clear ()
	{
		values.clear ();
	}

	inline int size ()
	{
		return values.size ();
	}

private:
#ifdef STL_WIN32
	std::vector<CCString> values;
#else
	vector<CCString> values;
#endif	// STL_WIN32
};

#ifdef CCVectorStr_cc

#endif	// CCVectorStr_cc 

#endif	// CCVECTORSTR_H

// EOF CCVectorStr.h
