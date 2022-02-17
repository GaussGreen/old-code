/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/inttypes.h
// Purpose:     defines int8, int16, int32 and int64 types
// Created:     17.11.2002
// RCS-ID:      $Id: inttypes.h,v 1.5 2005/04/05 12:32:55 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_INTTYPES_H_
#define _ITO33_INTTYPES_H_

/**
    @file ito33/inttypes.h
    @brief Declare fixed-size integer types.

    This is a subset of C99 standard inttypes.h header which is not yet
    available on all platforms.
 */

// get the sizes of the fundamental types
#if defined(HAVE_ITO33_CONFIG_H)
    #include "ito33/config.h"
#elif defined(_WIN32)
    #define SIZEOF_SHORT 2
    #define SIZEOF_INT 4
    #define SIZEOF_LONG 4

    #if defined(_WIN64)
      #define SIZEOF_SIZE_T 8
      #define SIZEOF_VOID_P 8
    #else // _WIN32
      #define SIZEOF_SIZE_T 4
      #define SIZEOF_VOID_P 4
    #endif // _WIN64/_WIN32
#endif

#if !defined(SIZEOF_SHORT) || \
    !defined(SIZEOF_INT) || \
    !defined(SIZEOF_LONG)
    #error "Unknown int types sizes."
#endif


namespace ito33
{

// define int8
typedef signed char int8;
typedef unsigned char uint8;

// define int16 if possible
#if SIZEOF_SHORT == 2
    typedef short int16;
    typedef unsigned short uint16;
#elif SIZEOF_INT == 2
    typedef int int16;
    typedef unsigned int uint16;
#endif

// define int32
#if SIZEOF_INT == 4
    typedef int int32;
    typedef unsigned int uint32;
#elif SIZEOF_LONG == 4
    typedef long int32;
    typedef unsigned long uint32;
#endif

// and finally define int64
#if SIZEOF_LONG == 8
    typedef long int64;
    typedef unsigned long uint64;
#elif defined(_MSC_VER)
    typedef __int64 int64;
    typedef unsigned __int64 uint64;
#elif defined(SIZEOF_LONG_LONG)
    typedef long long int64;
    typedef unsigned long long uint64;
#endif

} // namespace ito33

#endif // _ITO33_INTTYPES_H_

