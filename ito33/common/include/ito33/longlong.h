/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/longlong.h
// Purpose:     defines LongLong typedef: 64 bit int type
// Author:      Vadim Zeitlin
// Created:     25.05.03
// RCS-ID:      $Id: longlong.h,v 1.5 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/longlong.h
    @brief  Defines LongLong and ULongLong typedefs.

    This header defines LongLong type which represents a 64 bit integer type.
    Such type is often needed but unfortunately there is no standard name for
    it and, in fact, not even all compilers/platforms support it so we may have
    to write a class to emulate 64 bit ints ourselves (using 2 32 bit ints) if
    this is ever needed. For now we just provide portable names for the types
    defined by the compiler.

    LongLong represents a (signed) 64 bit type, ULongLong the unsigned
    version.

    ITO33_LL_FMT is the format spec to be used with printf() for printing
    LongLongs (it isn't always possible to print unsigned long longs like
    this so we don't define ITO33_ULL_FMT for now). Finally,
    ITO33_LL_SUFFIX is used by ITO33_LL below.
 */

#ifndef _ITO33_LONGLONG_H_
#define _ITO33_LONGLONG_H_

namespace ito33
{
#if defined(HAVE_LONGLONG)
  typedef long long LongLong;
  typedef unsigned long long ULongLong;

  #define ITO33_LL_SUFFIX ll
  #define ITO33_LL_FMT "ll"
#elif defined(_WIN32)
  #if defined(_MSC_VER)
    /// LongLong is a signed 64 bit integer type
    typedef __int64 LongLong;

    /// ULongLong is a signed 64 bit integer type
    typedef unsigned __int64 ULongLong;

    /// the suffix used with 64 bit constants, only for ITO33_LL
    #define ITO33_LL_SUFFIX i64

    /**
        Format string which should be used with @c printf() to print long
        longs.

        Usage example follows, note how the string concatenation is used
        to construct the full format string:
        \code
            LongLong ll = ITO33_LL(12345678901234567);
            printf("64 bit value = %" ITO33_LL_FMT "d\n", ll);
        \endcode
     */
    #define ITO33_LL_FMT "I64"
  #else // other compiler (at least Borland and CodeWarrior could be added)
    #error "don't know how to define LongLong with this compiler."
  #endif
#else // unknown OS
  typedef long long LongLong;
  typedef unsigned long long ULongLong;

  #error "don't know how to define LongLong on this platform."
#endif // OS

/// helper for ITO33_LL_HELPER_2
#define ITO33_LL_HELPER_1(x, s) x ## s

/// helper for ITO33_LL
#define ITO33_LL_HELPER_2(x, s) ITO33_LL_HELPER_1(x, s)

/**
  This macro allows defining long long constants in a portable way.
 */
#define ITO33_LL(x) ITO33_LL_HELPER_2(x, ITO33_LL_SUFFIX)

} // namespace ito33

#endif // _ITO33_LONGLONG_H_

