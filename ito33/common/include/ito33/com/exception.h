/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/exception.h
// Purpose:     declares an exception class for COM errors
// Author:      Vadim Zeitlin
// Created:     13.01.03
// RCS-ID:      $Id: exception.h,v 1.9 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002,2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/exception.h
    @brief An exception class for COM errors.

    This file declares the COM::Exception class and COM_EXCEPTION macro for
    throwing the exceptions of this type.
 */

#ifndef _ITO33_COM_EXCEPTION_H_
#define _ITO33_COM_EXCEPTION_H_

#include "ito33/win32/exception.h"

#include "ito33/win32/winwrap.h"

#include <ole2.h>                   // for HRESULT

namespace ito33
{

namespace COM
{

/**
    Exception class for COM errors.

    COM::Exception stores the HRESULT of the failed function in addition to its
    name and the last Win32 error code (which are stored by the base class).
    Usually this HRESULT will be the same as the Win32 error code but not
    always, so this class should be used for all Windows API functions which
    explicitly return HRESULT as they might not set the last error.
 */
class Exception : public ito33::Win32::Exception
{
public:
  /**
      Ctor stores the HRESULT of the function which failed which can be
      accessed later by the code in the catch clause.

      In addition to this, it takes the function name which is needed for
      the base class.

      @param apiname the name of the API function which failed (this should
                     be a static string normally, as we don't copy it but
                     just hold on to the pointer)
      @param hr the OLE/COM error code of the function which failed
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  Exception(const char *apiname,
       HRESULT hr,
       const char *filename,
       size_t line,
       const char *function)
    : ito33::Win32::Exception(apiname, hr, filename, line, function),
     m_hr(hr)
  {
  }

  /// returns the HRESULT of the function which failed
  HRESULT GetHresult() const { return m_hr; }

  /// override the base class version to return the correct system error code
  virtual unsigned long GetSystemErrorCode() const { return GetHresult(); }

private:
  // COM error code
  HRESULT m_hr;
};

/**
    Use this macro to create a COM::Exception() object with minimal hassle.

    @param apiname the name of the API which has failed
    @param hr the error code
 */
#define COM_EXCEPTION(apiname, hr) \
  ::ito33::COM::Exception(apiname, hr, __FILE__, __LINE__, __FUNCTION__)


} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_EXCEPTION_H_



