/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/dlldecl.h
// Purpose:     declares ITO33_IHG_DLLDECL macro
// Author:      Vadim Zeitlin
// Created:     2005-02-18
// RCS-ID:      $Id: dlldecl.h,v 1.2 2005/03/31 10:36:12 pedro Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/dlldecl.h
    @brief Declaration of macros used to import or export classes from DLL.

    @sa ito33/dlldecl.h
 */

#ifndef _ITO33_IHG_DLLDECL_H_
#define _ITO33_IHG_DLLDECL_H_

#include "ito33/dlldecl.h"

/**
    @def ITO33_IHG_DLLDECL

    Macro which should be used in declaration of all public classes, functions
    and variables in IHG.

    @sa ITO33_DLLDECL
 */

// there are 3 different cases: building DLL, using DLL, or not DLL at all
#if defined(ITO33_IHG_BUILDING_DLL)
  #define ITO33_IHG_DLLDECL ITO33_DLLEXPORT
#elif defined(ITO33_IHG_USING_DLL)
  #define ITO33_IHG_DLLDECL ITO33_DLLIMPORT
#else
  #define ITO33_IHG_DLLDECL
#endif

#endif // _ITO33_IHG_DLLDECL_H_

