//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Format.hpp
//
//   Description : Formatting class - helps create strings
//
//   Author      : Mark A Robson
//
//   Date        : 23 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef FORMAT_HPP
#define FORMAT_HPP

#include <string>

CORE_BEGIN_NAMESPACE

using namespace std;


/** Set of utility methods for creating strings from doubles, ints etc */
class CORE_DLL Format{
public:
    /** Creates a string from a double using %f format */
    static string toString(double value);

    /** Creates a string from an int using %d format */
    static string toString(int    value);

    /** Creates a string from a size_t using %d format */
    static string toString(size_t value);

    /** Creates a string ("true" or "false") from an bool */
    static string toString(bool   value);

    /** Creates a string using sprintf conventions  - max string length
        is 512 */
    static string toString(const char* format, ...);

    /** Answers an all-uppercase version of the source string */
    static string toUppercase(const string& src);
};

CORE_END_NAMESPACE

#endif

