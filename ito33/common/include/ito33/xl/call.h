/////////////////////////////////////////////////////////////////////////////
// Name:        to33/ito33/XL/call.h
// Purpose:     wrappers around Excel4() function
// Author:      Vadim Zeitlin
// Created:     2006-04-14
// RCS-ID:      $Id: call.h,v 1.4 2006/05/05 12:18:05 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/XL/call.h
    @brief Wrappers around Excel SDK Excel4() function.

    This header declares a low-level XL::Call() function which is a direct
    wrapper of Excel4() and also several convenience functions.
 */

#ifndef _ITO33_XL_CALL_H_
#define _ITO33_XL_CALL_H_

/**
    Maximal number of arguments accepted by Excel4 call.

    It must be less than 30 because of Excel4 constraints.
 */
#ifndef XL_MAX_ARGUMENTS
  #define XL_MAX_ARGUMENTS 30
#endif

namespace ito33
{

namespace XL
{

/// Check if the return code is xlretSuccess and throw if it isn't
inline void ThrowIfFailed(int xlFunc, int xlError)
{
  if ( xlError != xlretSuccess )
  {
    throw EXCEPTION_MSG
          (
            xlError,
            String::Printf("Call to Excel function %x failed", xlFunc)
          );
  }
}

#ifdef DOXYGEN

/**
    Call Excel safely: unlike the standard Excel4() function, this one is
    overloaded for all possible number of parameters and so avoids the need to
    specify their number in the call itself (and suppresses the possibility of
    mistake).

    This function throws if the Excel functions fails, use CallNoThrow() if
    it's ok for the function to fail.

    @param xlFunc the Excel function to call, one of xlfXXX or xlcYYY
    @param rc the return value of the call, may be @c NULL if not needed
    @param arg1 the first argument
    ...
    @param argN the last argument
 */
inline void
Call(int xlFunc, XLOPER *rc, const Oper& arg1, ... const Oper& argN);

/**
    Same as Call() but simply returns the error code from the Excel function
    which can be non 0 (indicating an error).
 */
inline int
CallNoThrow(int xlFunc, XLOPER *rc, const Oper& arg1, ... const Oper& argN);

#else // !DOXYGEN

#define MAKE_CALL_n(z, n, unused)                                             \
inline int                                                                    \
CallNoThrow(int xlFunc, XLOPER *rc                                            \
            BOOST_PP_ENUM_TRAILING_PARAMS(n, const Oper& arg))                \
{                                                                             \
  return Excel4(xlFunc, rc, n                                                 \
                BOOST_PP_ENUM_TRAILING_PARAMS(n, (const XLOPER *)arg));       \
}                                                                             \
                                                                              \
inline void                                                                   \
Call(int xlFunc, XLOPER *rc                                                   \
     BOOST_PP_ENUM_TRAILING_PARAMS(n, const Oper& arg))                       \
{                                                                             \
  ThrowIfFailed(xlFunc, CallNoThrow(xlFunc, rc                                \
                                    BOOST_PP_ENUM_TRAILING_PARAMS(n, arg)));  \
}

BOOST_PP_REPEAT_FROM_TO(0, XL_MAX_ARGUMENTS, MAKE_CALL_n, ~)

#undef MAKE_CALL_n

#endif // DOXYGEN/!DOXYGEN


/**
    Alert kinds.

    The values of the elements in this enum correspond to the constants used by
    Excel, don't change them.
 */
enum AlertKind
{
  /// Shows a question dialog box with "Ok" and "Cancel" buttons
  Alert_Question = 1,

  /// Shows a dialog box with "Ok" button only and information icon
  Alert_Info,

  /// Shows a dialog box with "Ok" button only and error icon
  Alert_Error
};

/**
    Generic Excel alert function, for internal use only.
 */
inline void Alert(const std::string& msg, AlertKind kind, XLOPER *rc = NULL)
{
  XL::Call(xlcAlert, rc, msg, kind);
}

/**
    Show an Excel dialog with a question.

    @param msg the message to show to the user
    @return true if the user chose the "Ok" button or false if "Cancel"
 */
inline bool Question(const std::string& msg)
{
  Oper rc;
  Alert(msg, Alert_Question, ByRef(rc));
  return rc == TRUE;
}

/**
    Show an Excel information dialog with the given text.
 */
inline void Info(const std::string& msg)
{
  Alert(msg, Alert_Info);
}


/**
    Show an Excel error dialog with the given text.
 */
inline void Error(const std::string& msg)
{
  Alert(msg, Alert_Error);
}

} // namespace XL

} // namespace ito33

#endif // _ITO33_XL_CALL_H_

