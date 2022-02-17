/////////////////////////////////////////////////////////////////////////////
// Name:        utils/exception_base.cpp
// Purpose:     implements ito33::Exception class
// Author:      Vadim Zeitlin
// Created:     12.02.03
// RCS-ID:      $Id: exception_base.cpp,v 1.8 2005/03/28 12:58:56 zeitlin Exp $
// Copyright:   (c) 2002-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/exception.h"

#include "ito33/gettext.h"

#include "ito33/crashrpt.h"

using namespace ito33;

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// Exception ctor and dtor
// ----------------------------------------------------------------------------

Exception::Exception(int errorCode,
          const std::string& msg,
          const char *filename,
          size_t line,
          const char *function)
    : m_errorCode(errorCode),
     m_errorMsg(msg),
     m_filename(filename),
     m_line(line),
     m_function(function)
{
  ITO33_GENERATE_BACKTRACE();
}

Exception::~Exception()
{
}

// ----------------------------------------------------------------------------
// GetFullMessage() helper function
// ----------------------------------------------------------------------------

std::string ito33::Exception::GetLocation() const
{
  std::string msg;

  if ( !m_function.empty() )
  {
    msg += String::Printf(TRANS("in %s() "), m_function.c_str());
  }

  msg += String::Printf(TRANS("at line %lu of file %s"),
             (unsigned long)m_line, m_filename.c_str());

  return msg;
}

// ----------------------------------------------------------------------------
// GetFullMessage() itself
// ----------------------------------------------------------------------------

std::string ito33::Exception::GetFullMessage() const
{
  std::string msg;

  msg = String::Printf(TRANS("Error %d occured "), m_errorCode);

  msg += GetLocation();

  if ( !m_errorMsg.empty() )
  {
    msg += ":\n\n";
    msg += m_errorMsg;
  }

  return msg;
}

