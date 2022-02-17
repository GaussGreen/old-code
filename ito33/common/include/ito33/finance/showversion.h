/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/showversion.h
// Purpose:     show version function
// Author:      ZHANG Yunzhi
// Created:     2004-05-04
// RCS-ID:      $Id: showversion.h,v 1.5 2004/11/10 15:53:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_FINANCE_SHOWVERSION_H_
#define _ITO33_FINANCE_SHOWVERSION_H_

#include "ito33/string.h"
#include "ito33/finance/version.h"

namespace ito33
{
  namespace finance
  {


  /**
    Gets the version

    @return version
    */
  inline std::string ShowVersion()
  {
    return ITO33_ITO_VERSION_DOT_STRING;
  }

  }
}

#endif // _ITO33_FINANCE_SHOWVERSION_H_

