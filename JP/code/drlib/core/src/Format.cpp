//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Format.cpp
//
//   Description : Formatting class - helps create strings
//
//   Author      : Mark A Robson
//
//   Date        : 23 Jan 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/coreConfig.hpp"
#include "edginc/Format.hpp"

#include <cstdarg>
#include <cstdio>
#include <cmath>

#if defined( _MSC_VER )
#define vsnprintf _vsnprintf
#endif

CORE_BEGIN_NAMESPACE

using namespace std;

/** Creates a string from a double using %f format */
string Format::toString(double value){
    char buffer[50];
    if (fabs(value) < 1.e10) {
        sprintf(buffer, "%f", value);
    } else {
        sprintf(buffer, "%e", value);
    }
    return string(buffer);
}

/** Creates a string from an int using %d format */
string Format::toString(int    value){
    char buffer[50];
    sprintf(buffer, "%d", value);
    return string(buffer);
}

/** Creates a string from a size_t using %d format */
string Format::toString(size_t value) {
	return toString((int)value);
}

/** Creates a string ("true" or "false") from an bool */
string Format::toString(bool   value){
    return string(value? "true": "false");
}

/** Creates a string using sprintf conventions  - max string length
    is 512 */
string Format::toString(const char* format, ...){
    char buffer[512];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, 511, format, args);
    buffer[511] = 0;
    va_end(args);
    return string(buffer);
}

string Format::toUppercase(const string& src) {
    string usrc(src);
    transform(usrc.begin(), usrc.end(), usrc.begin(), (int(*)(int))toupper);

    return usrc;
}

CORE_END_NAMESPACE

