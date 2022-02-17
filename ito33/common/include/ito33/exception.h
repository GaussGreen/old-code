/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/exception.h
// Purpose:     declares the base class of all ITO33 exceptions hierarchy
// Author:      Vadim Zeitlin
// Created:     16.12.02
// RCS-ID:      $Id: exception.h,v 1.18 2006/06/15 16:36:48 wang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/exception.h
    @brief Base class of all exceptions used in ITO33 code.

    All exceptions thrown by ITO33 code should derive from the Exception class
    defined here. This base class is trivial and only carries information about
    the error code and also the location of the exception.
 */

#ifndef _ITO33_EXCEPTION_H_
#define _ITO33_EXCEPTION_H_

#include "ito33/string.h"

namespace ito33
{

/**
    Exception class is the base of ITO33 exceptions hierarchy.
 */
class ITO33_DLLDECL Exception
{
public:
  /**
      Ctor creates the exception with the given error code and message.

      Use EXCEPTION() or EXCEPTION_MSG() macros to avoid having to write
      __FILE__, __LINE__ and __FUNCTION__ when constructing an exception
      object.

      @param errorCode the error code
      @param msg the human readable error message (may be empty)
      @param filename the name of the file where the exception occured
             (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
             (__FUNCTION__ is unfortunately not yet supported by all
              compilers so this argument is left empty for them)
   */
  Exception(int errorCode,
       const std::string& msg,
       const char *filename,
       size_t line,
       const char *function);

  /// virtual dtor as for any base class
  virtual ~Exception();


  /// @name Accessors
  //@{

  /// get the exception error code
  int GetErrorCode() const;

  /// get the human readable error message
  const char *GetErrorMessage() const;

  /// get the filename
  const char *GetFileName() const;

  /// get the line number
  size_t GetLineNumber() const;

  /// get the function name (may be empty)
  const char *GetFunctionName() const;

  //@}

  /**
      Returns the full error message combining all available information.

      This function is virtual and is meant to be overridden by the derived
      classes which may wish to add more details. We only show the error
      code, message (if any) and the location of the exception.
   */
  virtual std::string GetFullMessage() const;

  /**
      Returns the system code if we carry information about it.

      This is really just a hook for derived classes. The base class always
      returns 0 from here which means that it doesn't have information about
      the system error associated with this exception.
   */
  virtual unsigned long GetSystemErrorCode() const;

protected:
  /// part of GetFullMessage() returning the exception location
  std::string GetLocation() const;

private:
  // the information carried by the exception
  int m_errorCode;
  std::string m_errorMsg;

  // the exception location
  std::string m_filename;
  size_t m_line;
  std::string m_function;
};

/**
    EXCEPTION macro allows to create an Exception object without having to
    manually specifiy the exception location.

    This macro expands into an object of class Exception which is intentionally
    left without scope qualification -- this allows to use it inside the code
    of some class Foo to throw exceptions of class Foo::Exception.

    Note that not all compilers support __FUNCTION__ macro used below (only gcc
    and VC++ 7 are known to support it right now, in fact), and so the function
    field of the exception object risks to be empty in some cases.

    @param code the error code
 */
#define EXCEPTION(code) EXCEPTION_MSG(code, code.GetMessage())

/**
    EXCEPTION_MSG is similar to EXCEPTION but also allows to specify the error
    message and not only error code.

    The error message may be shown to the user so it must be human readable (in
    particular translated) and "user understandable".

    @param code the error code
    @param msg the error message
 */
#define EXCEPTION_MSG(code, msg) \
  Exception(code, msg, __FILE__, __LINE__, __FUNCTION__)

/**
    Catches any exception and logs it. This macro is meant for use in
    destructors, where it's necessary to ensure that no exception is left
    uncaught.

    @code
    try
    {
      ...do something that may throw...
    }
    CATCH_AND_LOG_EXCEPTION(MY_TRACE)
    @endcode

    @param logfunc Function, object or macro that is called to log an
                   exception. This must support printf-like syntax, i.e. it
                   must be possible to write @a logfunc("foo %s...", bar, ...).
 */
#define CATCH_AND_LOG_EXCEPTION(logfunc)                                    \
  catch ( const ito33::Exception& e )                                       \
  {                                                                         \
    logfunc("*** Exception caught and ignored in %s: %s ***",               \
            __FUNCTION__, e.GetFullMessage().c_str());                      \
  }                                                                         \
  catch ( const std::exception& e )                                         \
  {                                                                         \
    logfunc("*** Exception caught and ignored in %s: %s ***",               \
            __FUNCTION__, e.what());                                        \
  }                                                                         \
  catch ( ... )                                                             \
  {                                                                         \
    logfunc("*** Exception caught and ignored in %s ***",                   \
            __FUNCTION__);                                                  \
  }


// ============================================================================
// implementation only from now on
// ============================================================================

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

inline
int Exception::GetErrorCode() const
{
  return m_errorCode;
}

inline
const char *Exception::GetErrorMessage() const
{
  return m_errorMsg.c_str();
}

inline
const char *Exception::GetFileName() const
{
  return m_filename.c_str();
}

inline
size_t Exception::GetLineNumber() const
{
  return m_line;
}

inline
const char *Exception::GetFunctionName() const
{
  return m_function.c_str();
}

inline
unsigned long Exception::GetSystemErrorCode() const
{
  return 0;
}

} // namespace ito33

#endif // _ITO33_EXCEPTION_H_

