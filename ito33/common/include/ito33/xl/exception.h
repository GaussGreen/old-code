/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/XL/exception.h
// Purpose:     exception corresponding to Excel-related errors
// Author:      Vadim Zeitlin
// Created:     2006-04-14
// RCS-ID:      $Id: exception.h,v 1.2 2006/04/14 21:08:41 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/XL/exception.h
    @brief Exception thrown when an Excel function fails.
 */

#ifndef _ITO33_XL_EXCEPTION_H_
#define _ITO33_XL_EXCEPTION_H_

#include "ito33/error.h"
#include "ito33/exception.h"

/// The error code used for Excel-related errors.
extern const ito33::Error ITO33_EXCEL_ERROR;

namespace ito33
{

namespace XL
{

/**
    Exception thrown when a call to Excel fails.
 */
class Exception : public ::ito33::Exception
{
public:
  /**
      Ctor for the Exception object.

      Use the standard EXCEPTION macro to create Exception objects, this
      frees you from having to type __FILE__, __LINE__ and __FUNCTION__

      The @a xlError should be one of xlretXXX constants (not xlretSuccess).
    */
  Exception(int xlError,
            const std::string& msg,
            const char *filename,
            size_t line,
            const char *function)
    : ito33::Exception(ITO33_EXCEL_ERROR, msg, filename, line, function),
      m_xlError(xlError)
  {
  }

  /**
      Return the Excel error code.
   */
  int GetXLCode() const { return m_xlError; }

  /**
      Override the base class version to add information about Excel error code.
   */
  virtual std::string GetFullMessage() const
  {
    return ito33::Exception::GetFullMessage() +
              String::Printf("\nExcel error %x", m_xlError);
  }

private:
  int m_xlError;
};


} // namespace XL

} // namespace ito33

#endif // _ITO33_XL_EXCEPTION_H_

