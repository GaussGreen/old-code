/***************************************************************
 * Module:	DR Library C++ Utilities
 * Submodule:	
 * File:	kplatdep.h
 * Function:	Standard Definition and Include files.
 * Author:	Andrew Chou (revised c. Daher).
 * Revision:	$Header$
 ***************************************************************/
#ifndef _kplatdep__h
#define _kplatdep__h


////////////////////////////////////////////////////////////////
//
// Microsoft Visual C++
//
////////////////////////////////////////////////////////////////
#if defined (_MSC_VER) || defined(_WIN32)


#include <exception>
#include <string>
#include <iostream>	
#include <fstream>
#include <strstream>
#include <sstream>

#include <algorithm>
#include <map>
#include <vector>
#include <utility>
#include <set>
//#include <list>
#include <iomanip>		/* setprecision and scientific format */

using namespace std;
#pragma warning(disable: 4786) // turns off annoying warning msg
#pragma warning(disable: 4800)	// turns off int to bool warning




// Defines the platform dependent include files
// On MSC, I use the stl libraries

#define IOSTREAM_H <iostream>
#define IOMANIP_H  <iomanip>
#define SSTREAM_H  <sstream>
#define STDARG_H   <cstdarg>
#define STDIO_H    <cstdio>
#define STRING_H   <cstring> 
#define FSTREAM_H  <fstream>
#define STDLIB_H  <cstdlib>

// Define some platform independent terms
#define OSTRINGSTREAM ostringstream
#define ISTRINGSTREAM istringstream

#include <cmath>
#include <ctime>
using namespace std;
#pragma warning(disable: 4786) // turns off annoying warning msg

const char PATH_DELIMITER[] = "\\";	// path delimiter for file names
const char HOMEPATH[] = "HOMEDRIVE";

// Solaris does not support default template args, so I'm stuck doing this.

#define KMap(a, b)		map < a, b >
#define KMultimap(a, b)		multimap < a, b >
#define KPair(a, b)		pair < a, b >
#define KVector(a)		vector < a >
#define KList(a)		list < a >
#define KSet(a)			set < a >

#define K_MEMBER_TEMPLATES

typedef exception EXCEPTION;
#define MYALLOC(a) allocator< a >


////////////////////////////////////////////////////////////////
//
// Solaris CC compilier Version 5.0
//
////////////////////////////////////////////////////////////////
#elif defined (__SUNPRO_CC)
#if (__SUNPRO_CC >= 0x0500)

const char PATH_DELIMITER[] = "/";
const char HOMEPATH[] = "HOME";


#include <string.h>
#include <iostream.h>	
#include <fstream.h>
#include <strstream.h>
#include <sstream>

#include <algorithm>
#include <map>
#include <utility>
#include <set>
#include <list>
#include <vector>
#include <iomanip>		/* setprecision format */

using namespace std;

#define KMap( a, b )		\
	std::map< 		\
		a, 		\
		b, 		\
		std::less<a>, 	\
		std::allocator<std::pair<const a, b> > >
#define KMultimap(a, b)	multimap < a, b, less < a > >
#define KPair(a, b)	pair < a, b >
#define KVector(a)	vector < a, std::allocator< a > >
#define KList(a)	list < a >


#define MYALLOC(a)	allocator_interface< allocator< a >, a >

#define KSet(a)		\
	set < a, less < a >, allocator_interface< allocator< a >, a > >



#include "kdrstring.h"

#undef K_MEMBER_TEMPLATES


////////////////////////////////////////////////////////////////
//
// Solaris CC compilier Version 4.0
//
////////////////////////////////////////////////////////////////
#elif (__SUNPRO_CC == 0x0420) 

# define __STL_NO_PROXY_ARROW_OPERATOR 1

#include "stl_config.h"

#include <iostream.h>	
#include <fstream.h>
#include <strstream.h>

#include <algorithm>
#include <map>
#include <utility>
#include <set>
#include <list>
#include <vector>
#include <iomanip>		/* setprecision format */

const char PATH_DELIMITER[] = "/";
const char HOMEPATH[] = "HOME";


#define KMap(a, b)		map < a, b, less < a > >
#define KMultimap(a, b)		multimap < a, b, less < a > >
#define KPair(a, b)		pair < a, b >
#define KVector(a)		vector < a >
#define KList(a)		list < a >
#define KSet(a)			set < a , less < a > >
#define MYALLOC(a)		alloc

#include "kdrstring.h"

#undef K_MEMBER_TEMPLATES


#endif	// __SUNPRO_CC
////////////////////////////////////////////////////////////////
//
// GCC Linux compilier
//
////////////////////////////////////////////////////////////////
#else


#include <iostream.h>	
#include <fstream.h>
#include <stdio.h>

#include <string>
#include <algorithm>
#include <map>
#include <utility>
#include <set>
#include <list>
#include <vector>
#include <strstream>
//#include <sstream>	// Problem on Linux gcc

using namespace std;

const char PATH_DELIMITER[] = "/";
const char HOMEPATH[] = "HOME";




#define KMap( a, b )		\
	std::map< 		\
		a, 		\
		b, 		\
		std::less<a>, 	\
		std::allocator<std::pair<const a, b> > >
#define KMultimap(a, b)	multimap < a, b, less < a > >
#define KPair(a, b)		pair < a, b >
#define KVector(a)		vector < a, std::allocator< a > >
#define KList(a)		list < a >
#define KSet(a)		\
	set < a, less < a >, allocator_interface< allocator< a >, a > >
#define MYALLOC(a)	allocator_interface< allocator< a >, a >

#define	setprecision(arg)	""	// Problem on Linux gcc


#endif




#endif // _kplatdep__h

