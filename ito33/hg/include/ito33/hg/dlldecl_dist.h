/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/dlldecl.h
// Purpose:     declares ITO33_HG_DLLDECL macro
// Created:     2005/04/06
// RCS-ID:      $Id: dlldecl_dist.h,v 1.2 2006/08/19 23:46:11 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/dlldecl.h
    @brief Declaration of macros used to import classes from DLL.

    The macros here are only used under Windows and expand to nothing under
    other platforms.
 */

#ifndef _ITO33_HG_DLLDECL_H_
#define _ITO33_HG_DLLDECL_H_

#include "ito33/dlldecl.h"

/**
    @def ITO33_HG_DLLDECL

    Macro which expands into compiler-specific construct used to import from
    DLL. Under platforms other than Windows expands to nothing.
 */
#ifdef _WIN32
  #define ITO33_HG_DLLDECL __declspec(dllimport)
#else // !_WIN32
  #define ITO33_HG_DLLDECL
#endif

/*
  Disable the warning on foo class needs to have dll-interface to be used by 
  clients of class Bar. It seems harmless.
*/
#ifdef _MSC_VER
  #pragma warning(disable:4251)
#endif // VC++

#endif // _ITO33_HG_DLLDECL_H_
