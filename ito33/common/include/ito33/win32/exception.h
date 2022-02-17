/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/exception.h
// Purpose:     exception thrown when a Win32 API fails
// Author:      Vadim Zeitlin
// Created:     24.12.02
// RCS-ID:      $Id: exception.h,v 1.10 2004/10/05 09:13:40 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/win32/exception.h
    @brief Class for exceptions used in Win32 code.

    This class is used when a Win32 API (unexpectedly) fails.
 */

#ifndef _ITO33_WIN32_EXCEPTION_H_
#define _ITO33_WIN32_EXCEPTION_H_

#include "ito33/exception.h"        // base class

namespace ito33
{

namespace Win32
{

/**
    Exception stores the Win32 last error and can return the associated error
    message.

    Exceptions of this class should be always thrown when a Win32 function
    unexpectedly fails.
 */
class Exception : public ito33::Exception
{
public:
  /**
      Ctor creates the exception object with the predefined ITO33_SYS_ERROR
      code and automatically remembers the Win32 last error code.

      The caller must only specify the name of the function which failed and
      the location of the exception.

      @param apiname the name of the API function which failed (this should
                     be a static string normally, as we don't copy it but
                     just hold on to the pointer)
      @param error the Win32 error code, if it is 0 ::GetLastError() is used
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  Exception(const char *apiname,
       int error,
       const char *filename,
       size_t line,
       const char *function);

  /// return the Win32 error code
  unsigned long GetWin32Error() const { return m_lastError; }

  /// returns the name of the API which failed
  const char *GetWin32Function() const { return m_apiname; }

  /// returns the error message corresponding to GetWin32Error()
  std::string GetWin32Message() const;

  /// override base class function to show the Win32 function name
  virtual std::string GetFullMessage() const;

  /// and to return the error code
  virtual unsigned long GetSystemErrorCode() const { return GetWin32Error(); }

private:
  // Win32 last error
  unsigned long m_lastError;

  // the name of the API which failed
  const char *m_apiname;
};

/**
    Use this macro to create a Win32::Exception() object without having to
    specify __FILE__, __LINE__ and __FUNCTION__.

    This macro must be used with the functions which return their error code
    and don't use GetLastError(). Otherwise, WIN32_EXCEPTION may be used for
    even less typing.

    @param apiname the name of the API which has failed
    @param error the Win32 error return by apiname
 */
#define WIN32_EXCEPTION_ERR(apiname, error) \
  ::ito33::Win32::Exception(apiname, error, __FILE__, __LINE__, __FUNCTION__)

/**
    Use this macro to create a Win32::Exception() object with minimal hassle.

    @param apiname the name of the API which has failed
 */
#define WIN32_EXCEPTION(apiname) WIN32_EXCEPTION_ERR(apiname, 0)

} // namespace Win32

} // namespace ito33

#endif // _ITO33_WIN32_EXCEPTION_H_

