/////////////////////////////////////////////////////////////////////////////
// Name:        win32/exception.cpp
// Purpose:     implementation of Win32::Exception class
// Author:      Vadim Zeitlin
// Created:     24.12.02
// RCS-ID:      $Id: exception.cpp,v 1.11 2004/11/23 13:03:53 wang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"

#include "ito33/debug.h"
#include "ito33/error.h"

#include "ito33/gettext.h"

#include "ito33/win32/exception.h"
#include "ito33/win32/winwrap.h"

#include <ctype.h>

extern const ito33::Error ITO33_SYS_ERROR;

using namespace ito33;

// ============================================================================
// Win32::Exception implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctor
// ----------------------------------------------------------------------------

Win32::Exception::Exception(const char *apiname,
              int error,
              const char *filename,
              size_t line,
              const char *function)
        : ito33::Exception(ITO33_SYS_ERROR, "", filename, line, function),
         m_apiname(apiname)
{
  // use the provided error if it is non 0, otherwise get the error ourselves
  // from the system
  m_lastError = error ? error : static_cast<unsigned long>(::GetLastError());
}

// ----------------------------------------------------------------------------
// message formatting
// ----------------------------------------------------------------------------

std::string
Win32::Exception::GetWin32Message() const
{
  std::string msg;

  // ask the system for the message to be put in lpMsgBuf
  LPVOID lpMsgBuf;
  if ( ::FormatMessage
     (
      FORMAT_MESSAGE_ALLOCATE_BUFFER |    // flags: alloc lpMsgBuf
      FORMAT_MESSAGE_FROM_SYSTEM |        //        get sys error msg
      FORMAT_MESSAGE_MAX_WIDTH_MASK,      //        no "\r\n"s please
      NULL,                               // source (ignored)
      m_lastError,                        // message id
      0,                                  // lang id (use best found)
      (LPTSTR)&lpMsgBuf,                  // output buffer
      0,                                  // min size of buffer to alloc
      NULL                                // arguments (none)
     ) == 0 || !lpMsgBuf ) {
    msg = "Unknown error";
  }
  else  { // FormatMessage() ok
    // copy it to our buffer and free memory
    msg = reinterpret_cast<const char *>(lpMsgBuf);

    ::LocalFree(lpMsgBuf);

    // normally this can't ever happen but check nevertheless
    if ( !msg.empty() ) {
      // returned string is capitalized which doesn't look good
      msg[0] = (char)(tolower(msg[0]));
    }
    else {
      FAIL( "::FormatMessage() returned empty message?" );
    }
  }

  return msg;
}

// ----------------------------------------------------------------------------
// overridden Exception functions
// ----------------------------------------------------------------------------

std::string Win32::Exception::GetFullMessage() const
{
  std::string msg;

  msg = String::Printf(TRANS("A Windows function %s() failed "),
            m_apiname);
  msg += GetLocation();
  msg += '\n';
  msg += String::Printf(TRANS("with error code %#lx (%s)"),
             m_lastError, GetWin32Message().c_str());

  return msg;
}

