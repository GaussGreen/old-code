/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/showversion.h
// Purpose:     file containing version information for IHG product
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: showversion.h,v 1.4 2004/11/10 15:54:57 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_IHG_SHOWVERSION_H_
#define _ITO33_IHG_SHOWVERSION_H_

#include "ito33/string.h"
#include "ito33/ihg/version.h"

namespace ito33
{
  namespace ihg
  {

  /**
    Gets the version

    @return version
    */
  inline std::string ShowVersion()
  {
    return ITO33_IHG_VERSION_DOT_STRING;
  }

  }
}

#endif // _ITO33_IHG_VERSION_H_

