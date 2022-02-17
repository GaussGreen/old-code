/////////////////////////////////////////////////////////////////////////////
// Name:        win32/net.cpp
// Purpose:     implements NetworkInitializer
// Author:      Vadim Zeitlin
// Created:     2004-11-20
// RCS-ID:      $Id: net.cpp,v 1.2 2004/11/23 13:03:53 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"
#include "ito33/error.h"
#include "ito33/exception.h"
#include "ito33/gettext.h"

#include "ito33/net.h"

#include <winsock2.h>

extern const ito33::Error ITO33_SYS_ERROR;

using namespace ito33;

// ============================================================================
// NetworkInitializer implementation
// ============================================================================

NetworkInitializer::NetworkInitializer()
{
    const WORD wVerNeeded = MAKEWORD(1, 1);
    WSADATA data;
    std::string msg;

    int rc = ::WSAStartup(wVerNeeded, &data);
    if ( rc != 0 )
    {
      msg = String::Printf(TRANS("error %ld"), rc);
    }
    else if ( data.wVersion != wVerNeeded )
    {
      msg = String::Printf(TRANS("version %u available but %u requested"),
                           data.wVersion, wVerNeeded);
    }

    if ( !msg.empty() )
    {
      throw EXCEPTION_MSG
            (
              ITO33_SYS_ERROR,
              String::Printf
              (
                TRANS("Windows network initialization failed (%s)"),
                msg.c_str()
              )
            );
    }
}

NetworkInitializer::~NetworkInitializer()
{
  if ( ::WSACleanup() != 0 )
  {
    // don't throw from dtor but still at least log the error somewhere
    ::OutputDebugString(
        String::Printf("WSACleanup failed: error code %d\r\n",
                       ::WSAGetLastError()).c_str()
      );
  }
}

