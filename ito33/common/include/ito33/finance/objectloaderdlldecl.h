/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/objectloaderdlldecl.h
// Purpose:     declares ITO33_OBJECTLOADER_DLLDECL macro
// Author:      ZHANG Yunzhi
// Created:     2006-08-23
// RCS-ID:      $Id: objectloaderdlldecl.h,v 1.1 2006/08/24 10:03:59 zhang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/objectloaderdlldecl.h
    @brief declares ITO33_OBJECTLOADER_DLLDECL macro macros used to import
           or export classes from DLL.
 */

#ifndef _ITO33_FINANCE_OBJECTLOADERDLLDECL_H_
#define _ITO33_FINANCE_OBJECTLOADERDLLDECL_H_

#include "ito33/dlldecl.h"

/**
    @def ITO33_OBJECTLOADER_DLLDECL

    Macro which should be used in declaration of all public classes, functions
    and variables.

    This macro expands to something which tells the compiler to export the
    declaration from DLL when building the DLL and to something which tells it
    to import the declaration from DLL otherwise. This ensures that the same
    header can be used both for building the DLL itself and the programs using
    it.

    Note that @c ITO33_USING_OBJECTLOADER_DLL @b must be defined when using
    objectloader.dll from
    user code.
 */

// there are 3 different cases: building DLL, using DLL, or not DLL at all
#if defined(ITO33_BUILDING_OBJECTLOADER_DLL)
  #define ITO33_OBJECTLOADER_DLLDECL ITO33_DLLEXPORT
#elif defined(ITO33_USING_OBJECTLOADER_DLL)
  #define ITO33_OBJECTLOADER_DLLDECL ITO33_DLLIMPORT
#else
  #define ITO33_OBJECTLOADER_DLLDECL
#endif

/*
  Disable the warning on foo class needs to have dll-interface to be used by 
  clients of class Bar. It seems harmless.
*/
#ifdef _MSC_VER
  #pragma warning(disable:4251)
#endif // VC++

#endif // _ITO33_FINANCE_OBJECTLOADERDLLDECL_H_
