#ifndef __drplatdep__h
/**@#-*/
#define __drplatdep__h
/**@#+*/

#if defined (_MSC_VER)	// Microsoft Visual C++

#include "kplatform.h"
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
#define OSTRINGSTREAM std::ostringstream
#define ISTRINGSTREAM std::istringstream

#include <exception>
#include <cmath>
#include <ctime>
//using namespace std;

const char PATH_DELIMITER[] = "\\";	// path delimiter for file names
const char HOMEPATH[] = "HOMEDRIVE";

// Solaris does not support default template args, so I'm stuck doing this.
#define KMap(a, b) std::map < a, b >
#define KVector(a) std::vector < a >
#define KList(a) std::list < a >
#define KSet(a) std::set < a >
#define K_MEMBER_TEMPLATES

typedef exception EXCEPTION;
//#define MYALLOC(a) std::allocator< a >

#endif

/**@#-*/

#if defined (__SUNPRO_CC) // Solaris CC compilier

#define IOSTREAM_H <iostream.h>
#define IOMANIP_H <iomanip.h>
#define SSTREAM_H <strstream.h>
#define STDARG_H <stdarg.h>
#define STDIO_H <stdio.h>
#define STRING_H  <string.h> 
#define FSTREAM_H <fstream.h>
#define STDLIB_H <stdlib.h>

class EXCEPTION
{
public:
    EXCEPTION() {}
    virtual ~EXCEPTION() {}
    virtual const char* what() const = 0;
};

#define OSTRINGSTREAM ostrstream
#define ISTRINGSTREAM istrstream

//#include "stl_config.h"
#include "math.h"
#include "time.h"
#include "ctype.h"

const char PATH_DELIMITER[] = "/";
const char HOMEPATH[] = "HOME";
#define KMap(a, b) std::map < a, b, std::less< a >, std::allocator< std::pair< const a, b > > >
#define KVector(a) std::vector < a, std::allocator< a > >
#define KList(a) std::list < a >
#define KSet(a) std::set < a , std::less< a >, std::allocator< a > >
//#define MYALLOC(a) std::alloc

#undef K_MEMBER_TEMPLATES

#endif

/**@#+*/

#endif



