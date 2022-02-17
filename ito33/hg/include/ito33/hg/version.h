/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/version.h
// Purpose:     file containing version information for HG product
// Created:     2005/01/13
// RCS-ID:      $Id: version.h,v 1.6 2006/02/28 10:24:56 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_HG_VERSION_H_
#define _ITO33_HG_VERSION_H_

#include "ito33/finance/version.h"

// ----------------------------------------------------------------------------
// version info: Just the same version as in finance
// ----------------------------------------------------------------------------

// major version
#define ITO33_HG_VERSION_MAJOR ITO33_ITO_VERSION_MAJOR 

// minor version
#define ITO33_HG_VERSION_MINOR ITO33_ITO_VERSION_MINOR 

// build number
#define ITO33_HG_VERSION_BUILD ITO33_ITO_VERSION_BUILD 

// ----------------------------------------------------------------------------
// other version-related helpers, normally don't need to be ever changed
// ----------------------------------------------------------------------------

#ifndef INNO_SETUP_RUNNING

// we must define this before including the header below
#define ITO33_REALLY_NEED_VERSION_H
#include "../../../../common/include/ito33/version.h"

// these are used in the .rc files
#define ITO33_HG_VERSION_STRING \
  ITO33_MAKE_VERSION_STRING(ITO33_HG_VERSION_MAJOR, \
                            ITO33_HG_VERSION_MINOR)

#define ITO33_HG_SERIES_DOT_STRING \
  ITO33_MAKE_SERIES_DOT_STRING(ITO33_HG_VERSION_MAJOR, \
                               ITO33_HG_VERSION_MINOR)

#define ITO33_HG_VERSION_DOT_STRING \
  ITO33_MAKE_VERSION_DOT_STRING(ITO33_HG_VERSION_MAJOR, \
                                ITO33_HG_VERSION_MINOR, \
                                ITO33_HG_VERSION_BUILD)

// check if the current version is equal or greater than the given one
#define ITO33_HG_CHECK_VERSION(x, y, z) \
        ITO33_COMPARE_VERSIONS(HG, x, y, z)

#endif // INNO_SETUP_RUNNING


#endif // #ifndef _ITO33_HG_VERSION_H_

