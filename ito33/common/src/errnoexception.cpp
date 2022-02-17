/////////////////////////////////////////////////////////////////////////////
// Name:        errnoexception.cpp
// Purpose:     implementation of ErrnoException class
// Author:      Vadim Zeitlin
// Created:     19.02.03
// RCS-ID:      $Id: errnoexception.cpp,v 1.6 2006/05/27 20:05:06 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  // 'strerror' was declared deprecated in VC8. But we can still use it in this
  // file
  #define _CRT_SECURE_NO_DEPRECATE 1
#endif

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"

#include "ito33/debug.h"
#include "ito33/error.h"
#include "ito33/gettext.h"

#include "ito33/errnoexception.h"

#include <errno.h>
#include <string.h>

extern const ito33::Error ITO33_SYS_ERROR;

using namespace ito33;

// ============================================================================
// ErrnoException implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctor
// ----------------------------------------------------------------------------

ErrnoException::ErrnoException(const char *apiname,
              int error,
              const char *filename,
              size_t line,
              const char *function)
        : ito33::Exception(ITO33_SYS_ERROR, "", filename, line, function),
         m_apiname(apiname)
{
  // use the provided errno if it is non 0, otherwise get it ourselves from
  // the system
  m_errno = error ? error : errno;
}

// ----------------------------------------------------------------------------
// message formatting
// ----------------------------------------------------------------------------

std::string
ErrnoException::GetCMessage() const
{
  std::string msg(strerror(m_errno));

  return msg;
}

// ----------------------------------------------------------------------------
// overridden Exception functions
// ----------------------------------------------------------------------------

std::string ErrnoException::GetFullMessage() const
{
  std::string msg;

  msg = String::Printf(TRANS("A standard function %s() failed "),
            m_apiname);
  msg += GetLocation();
  msg += '\n';
  msg += String::Printf(TRANS("with error code %d (%s)"),
             m_errno, GetCMessage().c_str());

  return msg;
}
