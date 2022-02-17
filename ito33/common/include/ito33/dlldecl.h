/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/dlldecl.h
// Purpose:     declares ITO33_DLLDECL macro
// Author:      Vadim Zeitlin
// Created:     2005-02-17
// RCS-ID:      $Id: dlldecl.h,v 1.5 2006/08/19 18:41:24 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/dlldecl.h
    @brief Declaration of macros used to import or export classes from DLL.

    The macros here are only used under Windows and expand to nothing under
    other platforms.
 */

#ifndef _ITO33_DLLDECL_H_
#define _ITO33_DLLDECL_H_

/**
    @def ITO33_DLLEXPORT

    Macro which expands into compiler-specific construct used to export from
    DLL. Under platforms other than Windows expands to nothing.
 */

/**
    @def ITO33_DLLIMPORT

    Macro which expands into compiler-specific construct used to import from
    DLL. Under platforms other than Windows expands to nothing.
 */

#ifdef _WIN32
  // this is, of course, not standard, but nowadays most Windows compilers do
  // support this syntax
  #define ITO33_DLLEXPORT __declspec(dllexport)
  #define ITO33_DLLIMPORT __declspec(dllimport)
#else // !_WIN32
  #define ITO33_DLLEXPORT
  #define ITO33_DLLIMPORT
#endif


/**
    @def ITO33_DLLDECL

    Macro which should be used in declaration of all public classes, functions
    and variables.

    This macro expands to something which tells the compiler to export the
    declaration from DLL when building the DLL and to something which tells it
    to import the declaration from DLL otherwise. This ensures that the same
    header can be used both for building the DLL itself and the programs using
    it.

    Note that @c ITO33_USING_DLL @b must be defined when using ITO33 DLLs from
    user code.
 */

// there are 3 different cases: building DLL, using DLL, or not DLL at all
#if defined(ITO33_BUILDING_DLL)
  #define ITO33_DLLDECL ITO33_DLLEXPORT
#elif defined(ITO33_USING_DLL)
  #define ITO33_DLLDECL ITO33_DLLIMPORT
#else
  #define ITO33_DLLDECL
#endif

/*
  Disable the warning on foo class needs to have dll-interface to be used by 
  clients of class Bar. It seems harmless.
*/
#ifdef _MSC_VER
  #pragma warning(disable:4251)
#endif // VC++

#endif // _ITO33_DLLDECL_H_
