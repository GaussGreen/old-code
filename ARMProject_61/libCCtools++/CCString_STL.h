/* -------------------------------------------------------------------------

   File: CCString_STL.h
   Path: /home/julia/projets/dev/Common/libCCtools++/SCCS/s.CCString_STL.h
   Description: header CCString_STL
   Created: 98/06/02
   Author: Charles-Emmanuel MUSY
   Modified: 99/07/21 15:53:00
   Last maintained by: Charles-Emmanuel MUSY
   Revision: 1.6

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#ifndef CCSTRING_STL_H
#define CCSTRING_STL_H

#include <CCcommon.h>
SCCS_ID (CCString_STL_h_SccsId, "@(#)CCString_STL.h	1.6, modified 99/07/21");

#ifndef STL_WIN32
#include <hash_map.h>
#endif	// STL_WIN32

#include <CCString.h>

#ifndef STL_WIN32
__STL_TEMPLATE_NULL struct hash<CCString>
{
  size_t operator() (CCString string) const { return __stl_hash_string (string); }
};

__STL_TEMPLATE_NULL struct hash<CCString&>
{
  size_t operator() (const CCString& string) const { return __stl_hash_string (string); }
};
#endif	// STL_WIN32

struct CCeqstring
{
        bool operator () (const CCString &s1, const CCString &s2) const { return (s1 == s2 ? true : false); }
};

struct eqstr
{
	bool operator () (const char* s1, const char* s2) const { return (strcmp (s1, s2) == 0); }
};

struct ltstr
{
	bool operator () (const char* s1, const char* s2) const { return (strcmp (s1, s2) < 0); }
};

#endif 	// CCSTRING_STL_H

// EOF CCString_STL.h
