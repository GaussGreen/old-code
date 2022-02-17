/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/version.h
// Purpose:     file containing version information for basic ito library
// Author:      ZHANG Yunzhi
// Created:     2004-07-02
// RCS-ID:      $Id: version.h,v 1.23 2006/03/01 17:42:59 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_FINANCE_VERSION_H_
#define _ITO33_FINANCE_VERSION_H_

// ----------------------------------------------------------------------------
// version info: must be updated manually each time when the version changes
// ----------------------------------------------------------------------------

// major version
#define ITO33_ITO_VERSION_MAJOR 3

// minor version
#define ITO33_ITO_VERSION_MINOR 0

// build number
#define ITO33_ITO_VERSION_BUILD 0

// ----------------------------------------------------------------------------
// other version-related helpers, normally don't need to be ever changed
// ----------------------------------------------------------------------------

#ifndef INNO_SETUP_RUNNING

// we must define this before including the header below
#define ITO33_REALLY_NEED_VERSION_H
#include "ito33/version.h"

#define ITO33_ITO_SERIES_DOT_STRING \
  ITO33_MAKE_SERIES_DOT_STRING(ITO33_ITO_VERSION_MAJOR, \
                                ITO33_ITO_VERSION_MINOR)

#define ITO33_ITO_VERSION_DOT_STRING \
  ITO33_MAKE_VERSION_DOT_STRING(ITO33_ITO_VERSION_MAJOR, \
                                ITO33_ITO_VERSION_MINOR, \
                                ITO33_ITO_VERSION_BUILD)

#endif // INNO_SETUP_RUNNING

#endif // _ITO33_FINANCE_VERSION_H_
