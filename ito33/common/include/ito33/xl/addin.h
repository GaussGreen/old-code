/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ito33/XL/addin.h
// Purpose:     everything you need to write Excel add-ins painlessly
// Author:      Vadim Zeitlin
// Created:     2006-03-21
// RCS-ID:      $Id: addin.h,v 1.36 2006/07/18 22:15:49 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ito33/XL/addin.h
    @brief Stuff needed to write Excel add-ins in C++.

    This header declares XL::AddIn class which should be used for implementing
    an XL add-in or an XLL.
 */

#ifndef _ITO33_XL_ADDIN_H_
#define _ITO33_XL_ADDIN_H_

#include "ito33/log.h"

#include "ito33/XL/oper.h"

#include <boost/preprocessor/repetition.hpp>

#ifdef _MSC_VER
  #pragma comment(lib, "xlcall32.lib")
#endif

/**
    Maximal number of parameters for an Excel-callable function.
 */
#ifndef XL_MAX_PARAMETERS
  #define XL_MAX_PARAMETERS 11
#endif

extern const ito33::Error ITO33_UNEXPECTED;


/**
    All functions called by Excel must be declared with this macro in front of
    its return value.

    For example:
      @code
        XLLEXPORT double MyExcelFunction() { ... }
      @endcode

    Note that, very strangely, stdcall doesn't seem to be needed here.
 */
#define XLLEXPORT extern "C" ITO33_DLLEXPORT

namespace ito33
{

namespace XL
{

/**
    Represents the entire add-in and provides methods for registering its
    functions.

    Each Excel add-in project must derive a class from this one and use
    XL_IMPLEMENT_ADDIN() macro somewhere after the derived class declaration.
 */
class AddIn
{
public:
  /**
      Get the unique instance of this class.

      Create() must have been called previously. If it hadn't, or if the global
      add-in object has been already destroyed, an exception is thrown.
   */
  static AddIn& Get()
  {
    if ( !ms_instance )
    {
      throw ::ito33::EXCEPTION_MSG(ITO33_UNEXPECTED, "No global add-in object.");
    }

    return *ms_instance;
  }

  /**
      Get the user-visible name of this add-in.

      This name must be specified as the ctor parameter.
   */
  const std::string& GetName() const { return m_name; }


  /**
      Return the error message corresponding to the currently propagating
      exception.

      This function may be only called from a catch clause, i.e. while handling
      an exception. It recognizes the standard exception classes and reports
      them using XL::Error() function. Handling of other, custom, exceptions
      may be added by overriding DoGetExceptionMessage() virtual function.

      This method can be called even when there is no global add-in instance
      but in this case DoGetExceptionMessage() is not used.
   */
  static std::string GetExceptionMessage();


  /**
      Called by the library to register all the functions exported by the
      add-in.

      Normally this is implemented by just calling Register() for all the
      functions we export.
   */
  virtual void RegisterAllFunctions() = 0;

  /**
      Called by the library when XL is closed or the addin is unloaded.

      Some cleanup needs to be done here and not in the dtor because doing some
      things, such as shutting down threads, from the dtor results in a
      deadlock inside the CRT because dtor is called from inside DllMain()
      which holds locks also used by the thread termination code.

      Notice however that the addin shouldn't attempt to unregister the
      functions it registered here as it doesn't work anyhow and it could
      actually be harmful because Excel sometimes can call xlAutoClose and then
      continue to use the add-in (this happens when the user starts closing
      Excel but then aborts it), so any cleanup done here must be reversible.
   */
  virtual void OnClose() { }

#ifdef DOXYGEN
  /**
      Register a function with Excel.

      @param func the function pointer (used for type deduction, not called here)
      @param name name of the function, supposed to be the same in Excel and in
                  the DLL; use XL_REGISTER_N macro to avoid having to specify it
      @param desc the description of the function shown in Excel
      @param arg1 the name of the first parameter shown in Excel
      @param desc1 the description of the first parameter shown in Excel
      ...
      @param argN the name of the last parameter shown in Excel
      @param descN the description of the last parameter shown in Excel
   */
  template <typename T1, ..., TN>
  void Register(T (*func)(T1, ..., TN),
                const char *name,
                const char *desc,
                const char *arg1,
                const char *desc1,
                ...,
                const char *argN,
                const char *descN);
#else // !DOXYGEN
  #define MAKE_ARG_DESC_n(z, n, unused)                                       \
    const char *arg##n,                                                       \
    const char *desc##n

  #define MAKE_PARAM_n(z, n, unused) arg##n, desc##n

  #define MAKE_REGISTER_n(z, n, unused)                                       \
  template <typename T BOOST_PP_ENUM_TRAILING_PARAMS(n, typename T)>          \
  void Register(T (*)(BOOST_PP_ENUM_PARAMS(n, T)),                            \
                const char *name,                                             \
                const char *desc                                              \
                BOOST_PP_ENUM_TRAILING(n, MAKE_ARG_DESC_n, ~))                \
  {                                                                           \
    DoRegister(name,                                                          \
               Concat(Traits<T>::Sig,                                         \
                      n                                                       \
                      BOOST_PP_ENUM_TRAILING_BINARY_PARAMS                    \
                      (                                                       \
                        n,                                                    \
                        Traits<T,                                             \
                        >::Sig BOOST_PP_INTERCEPT                             \
                      )),                                                     \
               desc,                                                          \
               n BOOST_PP_ENUM_TRAILING(n, MAKE_PARAM_n, ~));                 \
  }

  BOOST_PP_REPEAT_FROM_TO(0, XL_MAX_PARAMETERS, MAKE_REGISTER_n, ~)

  #undef MAKE_REGISTER_n
  #undef MAKE_ARG_DESC_n
  #undef MAKE_PARAM_n
#endif // DOXYGEN/!DOXYGEN


  /**
      Implementation of the real DllMain().

      @sa XL_IMPLEMENT_DLLMAIN
   */
  static BOOL DllMain(HMODULE hModule, DWORD dwReason, void *lpReserved);

protected:
  /**
      Default constructor.

      There should be only one AddIn object in the program. When this unique
      object is created, it sets itself as the unique AddIn instance in the
      program and can be later retrieved using Get().

      @param name the add-in name to be shown in Excel add-ins manager dialog
      @param category the category for the functions of this plugin (must not
                      be @c NULL nor empty)
   */
  AddIn(const char *name, const char *category)
    : m_name(name),
      m_xloCategory(category)
  {
    ms_instance = this;
    Call(xlGetName, XL::ByRef(m_xloDllName));
  }

  /**
      Destructor unregisters this object as the global add-in instance.

      The add-in is only destroyed by the library itself when the XLL is
      unloaded.
   */
  virtual ~AddIn()
  {
    ms_instance = NULL;
  }

private:
  /**
      Create the global add-in object.

      This method must be implemented by the add-in implementation and should
      return a heap-allocated instance of AddIn-derived class.

      It is called from the library itself when the XLL is loaded by Excel and
      must not be called directly.

      @sa XL_IMPLEMENT_ADDIN()
   */
  static AddIn *Create();

  /**
      May be overridden to report custom exceptions.

      Return a non-empty string if the exception was recognized by this
      function or an empty string otherwise (this is the default behaviour).

      The messages returned by this function are used in GetExceptionMessage()
      and hence by ReportException() and ReturnException().
   */
  virtual std::string DoGetExceptionMessage() { return std::string(); }


  // helper function concatenating n chars into a string
  static std::string Concat(char ch, int n, /* n "char" parameters */ ...);

  // do call xlfRegister
  //
  // throws if the registration failed
  void DoRegister(const char *name,
                  const std::string& sig,
                  const char *desc,
                  int n,
                  /* n pairs of "const char *" params */
                  ...);


  // the unique plugin instance
  static AddIn *ms_instance;


  // the name to be shown to the user
  const std::string m_name;

  // the category for the functions defined in this plugin
  const Oper m_xloCategory;

  // the full path to this DLL
  Oper m_xloDllName;


  // it creates and destroys us
  friend class AddInManager;

  NO_COPY_CLASS(AddIn);
};

/**
    Define a standard DLL entry function for an XLL.

    This macro must be used in the global scope in one of the main XLL project
    files (otherwise the standard dummy DllMain() may get linked in instead).
 */
#define XL_IMPLEMENT_DLLMAIN()                                                \
  BOOL APIENTRY DllMain(void *hModule, DWORD dwReason, void *lpReserved)      \
  {                                                                           \
    return ito33::XL::AddIn::DllMain((HMODULE)hModule, dwReason, lpReserved); \
  }

/**
    This macro must be used in the add-in implementation file.

    Notice that the @a name class must have a default ctor for this to work. If
    there is no default ctor you must implement the static AddIn::Create()
    method manually instead of using this macro.

    The macro expansion also includes a GetTheName() global function which can
    be used to access the global add-in object cast to the correct type.

    The macro must be used in the global scope as it uses XL_IMPLEMENT_DLLMAIN.

    @param name name of the AddIn-derived class
 */
#define XL_IMPLEMENT_ADDIN(name)                                              \
  XL_IMPLEMENT_DLLMAIN()                                                      \
  ::ito33::XL::AddIn *::ito33::XL::AddIn::Create() { return new name; }       \
  name& GetThe##name() { return static_cast<name&>(::ito33::XL::AddIn::Get()); }


/**
    Report the current exception using XL::Error() with the message returned by
    XL::AddIn::GetExceptionMessage().

    Notice that XL::Error() function can only be used from the commands, not
    functions, so use ReturnException() instead from the functions!

    As XL::AddIn::GetExceptionMessage(), this function can only be called
    during exception handling, i.e. from a catch clause.
 */
extern void ReportException();

/**
    Return an Oper containing the error message from the currently propagating
    exception as string value.

    ReturnException() can be used in the exported functions, unlike
    ReportException(). It must also be only called from a catch clause.
 */
extern XLOPER *ReturnException();


#ifdef DOXYGEN

/**
    @def XL_REGISTER_N

    Macro to register the add-in functions.

    The benefit of using the macro instead of calling AddIn::Register()
    directly is that it allows to avoid specifying both the function pointer
    and the function name in the common case when they are the same.

    @param func the name of the function to register
    @param desc the description of the function shown in Excel
    @param arg1 the name of the first parameter shown in Excel
    @param desc1 the description of the first parameter shown in Excel
    ...
    @param argN the name of the last parameter shown in Excel
    @param descN the description of the last parameter shown in Excel
 */
#define XL_REGISTER_N(func, desc, arg1, desc1, ..., argN, descN)

#else // !DOXYGEN

#define XL_REGISTER_0(func, desc)                                             \
  ito33::XL::AddIn::Get().Register(&func, #func, desc)

#define XL_REGISTER_1(func, desc, arg1, desc1)                                \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1)

#define XL_REGISTER_2(func, desc, arg1, desc1, arg2, desc2)                   \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1, arg2, desc2)

#define XL_REGISTER_3(func, desc, arg1, desc1, arg2, desc2, arg3, desc3)      \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3)

#define XL_REGISTER_4(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,      \
                      arg4, desc4)                                            \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4)

#define XL_REGISTER_5(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,      \
                      arg4, desc4, arg5, desc5)                                            \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4,     \
                                   arg5, desc5)

#define XL_REGISTER_6(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,      \
                      arg4, desc4, arg5, desc5, arg6, desc6)                  \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4,     \
                                   arg5, desc5, arg6, desc6)

#define XL_REGISTER_7(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,      \
                      arg4, desc4, arg5, desc5, arg6, desc6, arg7, desc7)     \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4,     \
                                   arg5, desc5, arg6, desc6, arg7, desc7)

#define XL_REGISTER_8(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,      \
                      arg4, desc4, arg5, desc5, arg6, desc6, arg7, desc7,     \
                      arg8, desc8)                                            \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4,     \
                                   arg5, desc5, arg6, desc6, arg7, desc7,     \
                                   arg8, desc8)

#define XL_REGISTER_9(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,      \
                      arg4, desc4, arg5, desc5, arg6, desc6, arg7, desc7,     \
                      arg8, desc8, arg9, desc9)                               \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4,     \
                                   arg5, desc5, arg6, desc6, arg7, desc7,     \
                                   arg8, desc8, arg9, desc9)

#define XL_REGISTER_10(func, desc, arg1, desc1, arg2, desc2, arg3, desc3,     \
                       arg4, desc4, arg5, desc5, arg6, desc6, arg7, desc7,    \
                       arg8, desc8, arg9, desc9, arg10, desc10)               \
  ito33::XL::AddIn::Get().Register(&func, #func, desc, arg1, desc1,           \
                                   arg2, desc2, arg3, desc3, arg4, desc4,     \
                                   arg5, desc5, arg6, desc6, arg7, desc7,     \
                                   arg8, desc8, arg9, desc9, arg10, desc10)

#endif // DOXYGEN/!DOXYGEN

/**
    Log category for the XLL messages.
 */
ITO33_DECLARE_LOG_CATEGORY(LogXLL);

/**
    Macro which can used for logging with "XLL" facility.
 */
#define XLL_TRACE    ITO33_TRACE_CATEGORY(::ito33::XL::LogXLL)

} // namespace XL

} // namespace ito33

#endif // _ITO33_XL_ADDIN_H_

