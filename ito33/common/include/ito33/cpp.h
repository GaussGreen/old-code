///////////////////////////////////////////////////////////////////////////////
// Name:        ito33/cpp.h
// Purpose:     Various preprocessor helpers
// Author:      Vadim Zeitlin
// Modified by:
// Created:     03.07.2003
// RCS-ID:      $Id: cpp.h,v 1.2 2004/11/12 17:02:18 zhang Exp $
// Copyright:   (c) 2003 ITO33-Solutions
///////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_CPP_H_
#define _ITO33_CPP_H_

/// Helper for CONCAT, don't use directly.
#define ITO33_CONCAT_HELPER(x, y) x ## y

/**
    Concatenates two tokens into a single identifier.

    @param x first token
    @param y second token
    @return concatenation of x and y, i.e. xy
 */
#define ITO33_CONCAT(x, y) ITO33_CONCAT_HELPER(x, y)

/**
    Helper for ITO33_MAKE_UNIQUE_NAME, don't use directly.

    Note that for normal compilers it is enough to use __LINE__ here but VC++
    doesn't handle it properly when "Edit and Continue" is on (/ZI switch) and
    instead introduces, starting in version 7.0, another macro which does work
    as expected.
 */
#if defined(_MSC_VER) && (_MSC_VER >= 1300)
    #define ITO33_CONCAT_LINE(name) ITO33_CONCAT(name, __COUNTER__)
#else
    #define ITO33_CONCAT_LINE(name) ITO33_CONCAT(name, __LINE__)
#endif

/**
    Helper macro to be able to define unique/anonymous objects: this works by
    appending the current line number to the given identifier to reduce the
    probability of the conflict (it may still happen if this is used in the
    headers, hence you should avoid doing it or provide unique prefixes then)

    @param name the identifier used as base for the "unique" name
    @return name with the line number appended
 */
#define ITO33_MAKE_UNIQUE_NAME(name) ITO33_CONCAT_LINE(name)

#endif // _ITO33_CPP_H_

