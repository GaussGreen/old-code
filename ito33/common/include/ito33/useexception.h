/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/useexception.h
// Purpose:     the only file to be included for using exception.
// Author:      ZHANG Yunzhi
// Created:     03/10/17
// RCS-ID:      $Id: useexception.h,v 1.8 2006/06/15 12:09:23 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/useexception.h
    @brief  the only file to be included for using exception. It offers
    as well some useful macros.

    It includes all related header files.
  */
#ifndef _ITO33_USEEXCEPTION_H_
#define _ITO33_USEEXCEPTION_H_

#include "ito33/exception.h"
#include "ito33/gettext.h"
#include "ito33/error.h"


/**
    CHECK_PTR macro checks "ptr" which can be either a standard pointer or a
    smart one. When ptr is 0, the macro throws an exception with given error
    "err" and given message "msg"
    
    The macro just returns "ptr". So it can be used in the following way

        p2 = CHECK_PTR_MSG(p1, ITO33_BAD_PARAM, "Setting invalid p1.")
    or
        return CHECK_PTR_MSG(p1, ITO33_BAD_DATA, "Accessing invalid data.")
 */
#define CHECK_PTR_MSG(ptr, err, msg) \
      (ptr ? ptr : (throw EXCEPTION_MSG(err, TRANS(msg)), ptr))

/**
  It is similar to CHECK_PTR_MSG macro. But the user doesn't have to specify the
  error message: the default message of "err" is used.
  */
#define CHECK_PTR(ptr, err) \
      (ptr ? ptr : (throw EXCEPTION_MSG(err, err.GetMessage()), ptr))

////////////////////////////////////////////////////////////////////////////////
// general macros similar to EXCEPTION_MSG
////////////////////////////////////////////////////////////////////////////////

/**
    EXCEPTION_MSG_1 is similar to EXCEPTION_MSG, but the error message is
    printf()-like format which takes one argument.

    Example,
    @code
      throw EXCEPTION_MSG_1(1, "Bad value of foo: %f.", 0.03);
    @endcode

    @param code the error code
    @param msg the error message which is a printf()-like format
    @param arg the only argument of the format (msg)
 */
#define EXCEPTION_MSG_1(code, msg, arg) \
  EXCEPTION_MSG(code, ito33::String::Printf(msg, arg).c_str())

/**
    EXCEPTION_MSG_2 is similar to EXCEPTION_MSG_1, but the error message is
    printf()-like format which takes two arguments.

    Example,
    @code
      throw EXCEPTION_MSG_2(2, "Bad value of %s: %f.", "foo", 0.03);
    @endcode

    @param code the error code
    @param msg the error message which is a printf()-like format
    @param arg1 the first argument of the format (msg)
    @param arg2 the second argument of the format (msg)
 */
#define EXCEPTION_MSG_2(code, msg, arg1, arg2) \
  EXCEPTION_MSG(code, ito33::String::Printf(msg, arg1, arg2).c_str())

/**
    EXCEPTION_MSG_3 is similar to EXCEPTION_MSG_1, but the error message is
    printf()-like format which takes three arguments.

    Example,
    @code
      throw EXCEPTION_MSG_3(3, "Bad date: %02d/%02d/%02d.", 72, 2, 30);
    @endcode

    @param code the error code
    @param msg the error message which is a printf()-like format
    @param arg1 the first argument of the format (msg)
    @param arg2 the second argument of the format (msg)
    @param arg3 the third argument of the format (msg)
 */
#define EXCEPTION_MSG_3(code, msg, arg1, arg2, arg3) \
  EXCEPTION_MSG(code, ito33::String::Printf(msg, arg1, arg2, arg3).c_str())

////////////////////////////////////////////////////////////////////////////////
// helper macros which should not be used directly
////////////////////////////////////////////////////////////////////////////////

/**
    CHECK_COND_MSG_GENERAL should not be used directly. Instead you should
    use CHECK_COND_MSG below or CHECK_MSG macros.
    
    The macro checks the condition "cond".
    When the condition is not satisfied, the macro throws an exception with 
    given error "err" and given message "msg"
 */
#define CHECK_COND_MSG_GENERAL(ns, cond, err, msg)                             \
      {                                                                        \
        if(!(cond)) throw ns::EXCEPTION_MSG(err, TRANS(msg));                  \
      } struct Dummy


/**
    CHECK_COND_GENERAL should not be used directly.

    It is similar to CHECK_COND_MSG_GENERAL macro. But the user doesn't have to
    specify the error message: the default message of "err" is used.
 */
#define CHECK_COND_GENERAL(ns, cond, err)                                      \
      {                                                                        \
        if(!(cond)) throw ns::EXCEPTION_MSG(err, err.GetMessage());            \
      } struct Dummy

/**
    CHECK_COND_1_GENERAL should not be used directly.

    It is similar to CHECK_COND_GENERAL macro. But the default message of err
    is a printf()-like format which takes 1 argument
 */
#define CHECK_COND_1_GENERAL(ns, cond, err, arg)                               \
      {                                                                        \
        if(!(cond)) throw ns::EXCEPTION_MSG_1(err, err.GetMessage(), arg);     \
      } struct Dummy

/**
    CHECK_COND_2_GENERAL should not be used directly.

    It is similar to CHECK_COND_GENERAL macro. But the default message of err
    is a printf()-like format which takes 2 arguments
 */
#define CHECK_COND_2_GENERAL(ns, cond, err, a1, a2)                            \
      {                                                                        \
        if(!(cond)) throw ns::EXCEPTION_MSG_2(err, err.GetMessage(), a1, a2);  \
      } struct Dummy

/**
    CHECK_COND_3_GENERAL should not be used directly.

    It is similar to CHECK_COND_GENERAL macro. But the default message of err
    is a printf()-like format which takes 3 arguments
 */
#define CHECK_COND_3_GENERAL(ns, cond, err, a1, a2, a3)                        \
      {                                                                        \
        if(!(cond))                                                            \
          throw ns::EXCEPTION_MSG_3(err, err.GetMessage(), a1, a2, a3);        \
      } struct Dummy

////////////////////////////////////////////////////////////////////////////////
// macros for checking condition and throwing exactly a ito33::Exception
////////////////////////////////////////////////////////////////////////////////

/**
    CHECK_COND macro checks the condition "cond" and throws exactly a
    ito33::Exception.

    When the condition is not satisfied, the macro throws an exception with 
    given error "err" and given message "msg".

    Note that msg will be translated by the macro.
 */
#define ITO_CHECK_COND_MSG(cond, err, msg) \
        CHECK_COND_MSG_GENERAL(ito33, cond, err, msg) 

/**
    It is similar to CHECK_COND_MSG macro. But the user doesn't have to specify the
    error message: the default message of "err" is used.
  */
#define ITO_CHECK_COND(cond, err) CHECK_COND_GENERAL(ito33, cond, err)

/**
    It is similar to CHECK_COND macro. But the error message of err
    should a printf()-like format which takes 1 argument arg
 */
#define ITO_CHECK_COND_1(cond, err, arg) \
        CHECK_COND_1_GENERAL(ito33, cond, err, arg)

/**
    It is similar to CHECK_COND macro. But the error message of err
    should a printf()-like format which takes 2 arguments arg1 and arg2
 */
#define ITO_CHECK_COND_2(cond, err, arg1, arg2) \
        CHECK_COND_2_GENERAL(ito33, cond, err, arg1, arg2)

/**
    It is similar to CHECK_COND macro. But the error message of err
    should a printf()-like format which takes 3 arguments
 */
#define ITO_CHECK_COND_3(cond, err, arg1, arg2, arg3) \
        CHECK_COND_3_GENERAL(ito33, cond, err, arg1, arg2, arg3)

////////////////////////////////////////////////////////////////////////////////
// macros for checking condition and throwing Exception of the current scope.
// Use these macros where there won't be any confusion.
////////////////////////////////////////////////////////////////////////////////

/**
    CHECK_COND_MSG is similar to ITO_CHECK_COND_MSG but it does not throw 
    exactly ito33::Exception. It just throws Exception, so the Exception in
    the scope where the macro is used.

    The macro checks the condition "cond".
    When the condition is not satisfied, the macro throws an exception with 
    given error "err" and given message "msg"
 */
#define CHECK_COND_MSG(cond, err, msg)                                         \
      {                                                                        \
        if(!(cond)) throw EXCEPTION_MSG(err, TRANS(msg));                      \
      } struct Dummy


/**
    CHECK_COND is similar to ITO_CHECK_COND but it does not throw 
    exactly ito33::Exception. It just throws Exception, so the Exception in
    the scope where the macro is used.

    The macro checks the condition "cond".
    When the condition is not satisfied, the macro throws an exception with 
    given error "err" and the default message of "err".
 */
#define CHECK_COND(cond, err)                                                  \
      {                                                                        \
        if(!(cond)) throw EXCEPTION_MSG(err, err.GetMessage());                \
      } struct Dummy

/**
    CHECK_COND_1 is similar to ITO_CHECK_COND_1 but it does not throw 
    exactly ito33::Exception. It just throws Exception, so the Exception in
    the scope where the macro is used.

    The macro checks the condition "cond".
    When the condition is not satisfied, the macro throws an exception with 
    given error code "err" and the default message of "err" which is 
    a printf()-like format which takes 1 argument.
 */
#define CHECK_COND_1(cond, err, arg)                                           \
      {                                                                        \
        if(!(cond)) throw EXCEPTION_MSG_1(err, err.GetMessage(), arg);         \
      } struct Dummy

/**
    CHECK_COND_2 is similar to ITO_CHECK_COND_2 but it does not throw 
    exactly ito33::Exception. It just throws Exception, so the Exception in
    the scope where the macro is used.
 */
#define CHECK_COND_2(cond, err, a1, a2)                                        \
      {                                                                        \
        if(!(cond)) throw EXCEPTION_MSG_2(err, err.GetMessage(), a1, a2);      \
      } struct Dummy

/**
    CHECK_COND_3 is similar to ITO_CHECK_COND_3 but it does not throw 
    exactly ito33::Exception. It just throws Exception, so the Exception in
    the scope where the macro is used.
 */
#define CHECK_COND_3(ns, cond, err, a1, a2, a3)                                \
      {                                                                        \
        if(!(cond))                                                            \
          throw EXCEPTION_MSG_3(err, err.GetMessage(), a1, a2, a3);            \
      } struct Dummy

#endif
