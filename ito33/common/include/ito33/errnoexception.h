/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/errnoexception.h
// Purpose:     exception thrown when a C library function fails
// Author:      Vadim Zeitlin
// Created:     19.02.03
// RCS-ID:      $Id: errnoexception.h,v 1.5 2005/04/02 14:03:08 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/errnoexception.h
    @brief Class for exceptions used together with C library functions.

    This class is used when a C library API (unexpectedly) fails.

    @sa Win32::Exception
 */

#ifndef _ITO33_ERRNOEXCEPTION_H_
#define _ITO33_ERRNOEXCEPTION_H_

#include "ito33/exception.h"        // base class

namespace ito33
{

/**
    Exception stores errno and can return the associated error message.

    Exceptions of this class should be always thrown when a C library function
    unexpectedly fails.
 */
class ErrnoException : public ito33::Exception
{
public:
  /**
      Ctor creates the exception object with the predefined ITO33_SYS_ERROR
      code and remembers the errno or the code given.

      The caller must only specify the name of the function which failed and
      the location of the exception.

      @param apiname the name of the API function which failed (this should
                     be a static string normally, as we don't copy it but
                     just hold on to the pointer)
      @param error the errno error code, if it is 0 errno is used
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  ErrnoException(const char *apiname,
       int error,
       const char *filename,
       size_t line,
       const char *function);

  /// return the errno error code
  unsigned long GetErrno() const { return m_errno; }

  /// returns the name of the API which failed
  const char *GetCFunction() const { return m_apiname; }

  /// returns the error message corresponding to GetErrno()
  std::string GetCMessage() const;

  /// override base class function to show the function name and error msg
  virtual std::string GetFullMessage() const;

private:
  // the error code
  int m_errno;

  // the name of the API which failed
  const char *m_apiname;
};

/**
    Use this macro to create a ErrnoException() object without having to
    specify __FILE__, __LINE__ and __FUNCTION__.

    This macro must be used with the functions which return their error code
    and don't use GetLastError(). Otherwise, ERRNO_EXCEPTION may be used for
    even less typing.

    @param apiname the name of the API which has failed
    @param error the errno error return by apiname
 */
#define ERRNO_EXCEPTION_ERR(apiname, error) \
  ::ito33::ErrnoException(apiname, error, __FILE__, __LINE__, __FUNCTION__)

/**
    Use this macro to create a ErrnoException() object with minimal hassle.

    @param apiname the name of the API which has failed
 */
#define ERRNO_EXCEPTION(apiname) ERRNO_EXCEPTION_ERR(apiname, 0)

} // namespace ito33

#endif // _ITO33_ERRNOEXCEPTION_H_

