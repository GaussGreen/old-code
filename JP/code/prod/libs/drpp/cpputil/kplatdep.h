#ifndef __drplatdep__h
/**@#-*/
#define __drplatdep__h
/**@#+*/

#if defined (_MSC_VER)	// Microsoft Visual C++

// Defines the platform dependent include files
// On MSC, I use the stl libraries
#define IOSTREAM_H <iostream>
#define IOMANIP_H <iomanip>
#define SSTREAM_H <sstream>
#define STDARG_H <cstdarg>
#define STDIO_H <cstdio>
#define STRING_H  <cstring> 
#define FSTREAM_H <fstream>
#define STDLIB_H <cstdlib>

// Define some platform independent terms
#define OSTRINGSTREAM ostringstream
#define ISTRINGSTREAM istringstream

#include <exception>
#include <cmath>
#include <ctime>
using namespace std;
#pragma warning(disable: 4786) // turns off annoying warning msg

const char PATH_DELIMITER[] = "\\";	// path delimiter for file names
const char HOMEPATH[] = "HOMEDRIVE";

// Solaris does not support default template args, so I'm stuck doing this.

#define KMap(a, b) map < a, b >
#define KVector(a) vector < a >
#define KList(a) list < a >
#define KSet(a) set < a >

#define K_MEMBER_TEMPLATES

typedef std::exception EXCEPTION;
#define MYALLOC(a) allocator< a >

// #endif

/**@#-*/

#else // defined (__SUNPRO_CC) // Solaris CC compilier

#if 0
#define IOSTREAM_H <iostream.h>
#define IOMANIP_H <iomanip.h>
#define SSTREAM_H <strstream.h>
#define STDARG_H <stdarg.h>
#define STDIO_H <stdio.h>
#define STRING_H  <string.h> 
#define FSTREAM_H <fstream.h>
#define STDLIB_H <stdlib.h>
#endif

#define IOSTREAM_H <iostream>
#define IOMANIP_H <iomanip>
#define SSTREAM_H <strstream>
#define STDARG_H <cstdarg>
#define STDIO_H <cstdio>
#define STRING_H  <cstring> 
#define FSTREAM_H <fstream>
#define STDLIB_H <cstdlib>

#if 0
#include "stl_config.h"
#include "math.h"
#include "time.h"
#include "ctype.h"
#endif

// Define some platform independent terms
#define OSTRINGSTREAM ostrstream
#define ISTRINGSTREAM istrstream

#include <exception>
#include <cmath>
#include <ctime>
using namespace std;

const char PATH_DELIMITER[] = "/";
const char HOMEPATH[] = "HOME";

#define KMap(a, b) map < a, b >
#define KVector(a) vector < a >
#define KList(a) list < a >
#define KSet(a) set < a >

#define K_MEMBER_TEMPLATES
typedef std::exception EXCEPTION;

#define MYALLOC(a) allocator< a >

#if 0
#define KMap(a, b) map < a, b, less < a >, alloc >
#define KVector(a) vector < a, alloc >
#define KList(a) list < a, alloc >
#define KSet(a) set < a, less < a >, alloc >
#define MYALLOC(a) alloc

#undef K_MEMBER_TEMPLATES
#endif

#endif

/**@#+*/

#endif



