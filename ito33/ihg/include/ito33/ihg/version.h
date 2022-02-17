/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/version.h
// Purpose:     file containing version information for IHG product
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: version.h,v 1.39 2006/02/28 10:24:36 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_IHG_VERSION_H_
#define _ITO33_IHG_VERSION_H_

#include "ito33/finance/version.h"

// ----------------------------------------------------------------------------
// version info: Just the same version as in finance
// ----------------------------------------------------------------------------

// major version
#define ITO33_IHG_VERSION_MAJOR ITO33_ITO_VERSION_MAJOR 

// minor version
#define ITO33_IHG_VERSION_MINOR ITO33_ITO_VERSION_MINOR

// build number
#define ITO33_IHG_VERSION_BUILD ITO33_ITO_VERSION_BUILD

// ----------------------------------------------------------------------------
// other version-related helpers, normally don't need to be ever changed
// ----------------------------------------------------------------------------

#ifndef INNO_SETUP_RUNNING

// we must define this before including the header below
#define ITO33_REALLY_NEED_VERSION_H
#include "../../../../common/include/ito33/version.h"

// these are used in the .rc files
#define ITO33_IHG_VERSION_STRING \
  ITO33_MAKE_VERSION_STRING(ITO33_IHG_VERSION_MAJOR, \
                            ITO33_IHG_VERSION_MINOR)

#define ITO33_IHG_SERIES_DOT_STRING \
  ITO33_MAKE_SERIES_DOT_STRING(ITO33_IHG_VERSION_MAJOR, \
                                ITO33_IHG_VERSION_MINOR)

#define ITO33_IHG_VERSION_DOT_STRING \
  ITO33_MAKE_VERSION_DOT_STRING(ITO33_IHG_VERSION_MAJOR, \
                                ITO33_IHG_VERSION_MINOR, \
                                ITO33_IHG_VERSION_BUILD)


// check if the current version is equal or greater than the given one
#define ITO33_IHG_CHECK_VERSION(x, y, z) \
    ITO33_COMPARE_VERSIONS(IHG, x, y, z)

#endif // INNO_SETUP_RUNNING


#endif // _ITO33_IHG_VERSION_H_
