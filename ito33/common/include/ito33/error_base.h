/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/error_base.h
// Purpose:     error codes returned by ITO33 functions
// Author:      ZHANG Yunzhi
// Created:     2006-02-09
// RCS-ID:      $Id: error_base.h,v 1.2 2006/06/13 11:48:23 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/error_base.h
  @brief This file contains the Error class declaration and some common error
         codes.

  We represent the error codes by class objects instead of the enum elements to
  make adding new error codes less painful (with an enum approach everything
  must be recompiled whenever a new error code is added even if is only used in
  one place).
 */
#ifndef _ITO33_ERROR_BASE_H_
#define _ITO33_ERROR_BASE_H_

#include "ito33/vector.h"
#include "ito33/common.h"
#include "ito33/gettext.h"

namespace ito33
{

/**
    Base Error class. 

    This class must not be used directly. To define a specific Error
    class, you should use DECLARE/IMPLEMENT_GIVEN_ERROR_CLASS macros
    below.
 */
class ErrorBase
{
public:

  /// Implicit conversion to integer.
  operator int() const { return m_code; }

protected:
  /**
      The (unique) error code.

      The error code is > 0 for errors, < 0 for warnings and == 0 only for
      ITO33_OK which indicates absence of any error.
   */
  int m_code;
};

//_____________________________________________________________________________
// to generate a humain readable description of a real Error class
//_____________________________________________________________________________
#ifdef DOXYGEN

/// prototype of a specific Error class
class Error
{
public:
  /// Start of error codes.
  enum { BASE = 1 };

  /// Creates a new error code with associated error message
  Error(const char* msg);

  /**
      Gets the error message corresponding to this error code.

      Error message is translated to the current language, if any.

      The returned pointer is static and shouldn't be deleted nor freed.
   */
  virtual const char *GetMessage() const;
}                       

#endif // #ifdef DOXYGEN

////////////////////////////////////////////////////////////////////////////////
// MORE GENERAL MACRO

/**
    Declares an error class.

    @param Error the name of error class
    @param offset the start of error code

    Example to declare MyError class in the name space myspace with error codes
    starting from 1000:
    @code
    namespace myspace
    {
      ITO33_DECLARE_GIVEN_ERROR_CLASS(MyError, 1000);
    }
    @endcode
 */
#define ITO33_DECLARE_GIVEN_ERROR_CLASS(Error, offset)                        \
class Error : public ito33::ErrorBase                                         \
{                                                                             \
public:                                                                       \
  Error(const char* msg);                                                     \
  const char *GetMessage() const;                                             \
                                                                              \
private:                                                                      \
  enum { BASE = offset };                                                     \
                                                                              \
  static std::vector<const char*> gs_messages;                                \
}                                                                             \

/**
    Declares the @b Error class. 

    @param offset the start of error code

    Example to declare myspace::Error class with error codes starting from 101
    @code
    namespace myspace
    {
      ITO33_DECLARE_ERROR_CLASS(101);
    }
    @endcode
 */
#define ITO33_DECLARE_ERROR_CLASS(offset) \
        ITO33_DECLARE_GIVEN_ERROR_CLASS(Error, offset)

/**
    Implments the given error class.

    Example to implement myspace::MyError class in source file.
    @code
    namespace myspace
    {
      ITO33_IMPLEMENT_GIVEN_ERROR_CLASS(MyError);
    }

    extern const myspace::MyError
      Err_First("First error."),
      ...
    @endcode

    @param Error the name of error class
 */
#define ITO33_IMPLEMENT_GIVEN_ERROR_CLASS(Error)                              \
std::vector<const char*> Error::gs_messages;                                  \
                                                                              \
Error::Error(const char* msg)                                                 \
{                                                                             \
  static int s_code = Error::BASE;                                            \
  m_code = s_code++;                                                          \
  gs_messages.push_back(msg);                                                 \
}                                                                             \
                                                                              \
const char* Error::GetMessage() const                                         \
{                                                                             \
  return ito33::GetTranslation(gs_messages[m_code - Error::BASE]);            \
} struct Dummy


/**
    Implements the @b Error class

    Example to implement myspace::Error class in source file.
    @code
    namespace myspace
    {
      ITO33_IMPLEMENT_ERROR_CLASS;
    }
    @endcode
 */
#define ITO33_IMPLEMENT_ERROR_CLASS ITO33_IMPLEMENT_GIVEN_ERROR_CLASS(Error)

} // namespace ito33

#endif /* _ITO33_ERROR_BASE_H_ */
