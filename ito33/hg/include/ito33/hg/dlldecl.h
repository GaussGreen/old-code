/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/dlldecl.h
// Purpose:     declares ITO33_HG_DLLDECL macro
// Created:     2005/04/06
// RCS-ID:      $Id: dlldecl.h,v 1.1 2005/04/06 14:11:43 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/dlldecl.h
    @brief Declaration of macros used to import or export classes from DLL.

    @sa ito33/dlldecl.h
 */

#ifndef _ITO33_HG_DLLDECL_H_
#define _ITO33_HG_DLLDECL_H_

#include "ito33/dlldecl.h"

/**
    @def ITO33_HG_DLLDECL

    Macro which should be used in declaration of all public classes, functions
    and variables in HG.

    @sa ITO33_DLLDECL
 */

// there are 3 different cases: building DLL, using DLL, or not DLL at all
#if defined(ITO33_HG_BUILDING_DLL)
  #define ITO33_HG_DLLDECL ITO33_DLLEXPORT
#elif defined(ITO33_HG_USING_DLL)
  #define ITO33_HG_DLLDECL ITO33_DLLIMPORT
#else
  #define ITO33_HG_DLLDECL
#endif

#endif // _ITO33_HG_DLLDECL_H_

