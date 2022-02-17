/* -------------------------------------------------------------------------
 
   File: CCVector.h
   Path: /home/julia/projets/dev/Common/libCCtools++/SCCS/s.CCVector.h
   Description: header CCVector
   Created: 99/05/11
   Author: Charles-Emmanuel MUSY
   Modified: 99/05/11 11:35:32
   Last maintained by: Charles-Emmanuel MUSY
   Revision: 1.3
 
   -------------------------------------------------------------------------
 
   Note:
 
   ------------------------------------------------------------------------- */
 
#ifndef CCVECTOR_H
#define CCVECTOR_H

#include <CCcommon.h>
SCCS_ID (CCVector_h_SccsId, "@(#)CCVector.h	1.3, modified 99/05/11");

#ifdef STL_WIN32
#include <vector>
#define VECTOR	std::vector
#else
#include <vector.h>
#define VECTOR	vector
#endif	// STL_WIN32

class CCVector
{
public:
	CCVector () : values ()
	{
	}

	CCVector (int n_size) : values (n_size)
	{
	}

	CCVector (VECTOR<double>& n_values) : values ()
	{
		set (n_values);
	}

	inline void set (VECTOR<double>& n_values)
	{
		clear ();
		for(int i = 0; i < n_values.size (); i++)
		{
			values.push_back (n_values[i]);
		}
	}	

	inline double operator [](int index)
	{
		return values[index];
	}

	inline void push_back (double value)
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
	VECTOR<double> values;
};

#ifdef CCVector_cc

#endif	// CCVector_cc 

#endif	// CCVECTOR_H

// EOF CCVector.h
