/////////////////////////////////////////////////////////////////////////////
// Name:        src/XL/addin.cpp
// Purpose:     implementation of Excel add-in support classes
// Author:      Vadim Zeitlin
// Created:     2006-04-11
// RCS-ID:      $Id: addin.cpp,v 1.12 2006/07/19 14:05:55 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/XL/addin.h"

#include "ito33/win32/module.h"

namespace ito33
{

namespace XL
{

// ----------------------------------------------------------------------------
// logging
// ----------------------------------------------------------------------------

ITO33_DEFINE_LOG_CATEGORY(LogXLL, "logxll");

// ----------------------------------------------------------------------------
// AddInManager: the class which can (being its friend) create/destroy AddIns
// ----------------------------------------------------------------------------

class AddInManager
{
public:
  // call this to create the global add-in object (safe to call multiple times)
  static void Initialize();

  // call this to destroy the global add-in object (safe to call even if
  // Initialize() hadn't been called but shouldn't be called more than once)
  static void Cleanup();

private:
  // initially false, set to true when Initialize() is called
  static bool ms_initialized;
};

// ============================================================================
// AddInManager implementation
// ============================================================================

bool AddInManager::ms_initialized = false;

/* static */
void AddInManager::Initialize()
{
  if ( !ms_initialized )
  {
    AddIn::Create();

    ms_initialized = true;
  }
}

/* static */
void AddInManager::Cleanup()
{
  if ( ms_initialized )
  {
    ms_initialized = false;

    // this will automatically reset AddIn::ms_instance to NULL
    delete &(AddIn::Get());
  }
}

// ============================================================================
// AddIn implementation
// ============================================================================

AddIn *AddIn::ms_instance = NULL;

/* static */
std::string AddIn::GetExceptionMessage()
{
  std::string msg;
  if ( ms_instance )
  {
    msg = ms_instance->DoGetExceptionMessage();
  }

  if ( msg.empty() )
  {
    try
    {
      throw;
    }
    catch ( XL::Exception& e )
    {
      msg = e.GetErrorMessage();
    }
    catch ( ito33::Exception& e )
    {
      msg = e.GetErrorMessage();
    }
    catch ( std::exception& e )
    {
      msg = e.what();
    }
    catch ( ... )
    {
      msg = "An unexpected error occured.";
    }
  }

  // this error message will be returned to Excel and there is a hard limit on
  // the Excel string length which is enforced by throwing an exception if it's
  // exceeded -- ensure that, whatever happens, we don't throw from the
  // exception handler from which we're usually called
  if ( msg.length() > 255 )
  {
    // so that the parts before and after "..." will have 126 characters each
    // and the total string length will be exactly 255
    static const size_t maxHalfLen = 126;

    msg = std::string(msg.begin(), msg.begin() + maxHalfLen + 1) +
          "..." +
          std::string(msg.end() - maxHalfLen, msg.end());
  }

  return msg;
}

/* static */
std::string AddIn::Concat(char ch, int n, /* n "char" parameters */ ...)
{
  std::string s;
  s.reserve(n + 1);
  s = ch;

  va_list argptr;
  va_start(argptr, n);

  for ( int i = 0; i < n; i++ )
  {
    s += va_arg(argptr, char);
  }

  va_end(argptr);

  return s;
}

void AddIn::DoRegister(const char *name,
                       const std::string& sig,
                       const char *desc,
                       int n,
                       /* n pairs of "const char *" params */
                       ...)
{
  XLL_TRACE("Registering %s(), signature \"%s\"", name, sig.c_str());

  enum
  {
    Arg_DllName,
    Arg_DllFunc,
    Arg_Signature,
    Arg_ExcelFunc,
    Arg_ArgNames,
    Arg_FuncType,
    Arg_FuncCategory,
    Arg_Unused,
    Arg_HelpTopic,
    Arg_Desc,
    Arg_FirstParamDesc
  };

  // prepare the arguments
  const XLOPER *opers[Arg_FirstParamDesc + XL_MAX_PARAMETERS];
  opers[Arg_DllName] = m_xloDllName.ByVal();

  std::auto_ptr<Oper> xloName(new Oper(name));
  opers[Arg_DllFunc] =
  opers[Arg_ExcelFunc] = xloName->ByVal();

  // append '#' to give the function macro sheet permissions: this is needed in
  // order to be able to call many built-in Excel functions and doesn't seem to
  // do any harm otherwise, so just do it for all functions
  std::auto_ptr<Oper> xloSig(new Oper(sig + "#"));
  opers[Arg_Signature] = xloSig->ByVal();

  static const Oper s_xloKindFunc(1); // this is a function, not a command
  opers[Arg_FuncType] = s_xloKindFunc.ByVal();

  opers[Arg_FuncCategory] = m_xloCategory.ByVal();
  opers[Arg_Unused] =
  opers[Arg_HelpTopic] = Oper::Missing().ByVal();

  std::auto_ptr<Oper> xloDesc(new Oper(desc));
  opers[Arg_Desc] = xloDesc->ByVal();

  va_list argptr;
  va_start(argptr, n);

  std::string argNames;
  std::auto_ptr<Oper> argDescs[XL_MAX_PARAMETERS];
  for ( int i = 0; i < n; i++ )
  {
    if ( i )
      argNames += ",";
    argNames += va_arg(argptr, char *);

    // due to a well-known bug in XL the last parameter description is
    // truncated and we need to add 2 filler characters to prevent this from
    // happening
    const char * const d = va_arg(argptr, char *);
    argDescs[i].reset(new Oper(i < n - 1 ? d
                                         : (std::string(d) + ". ").c_str()));
    opers[Arg_FirstParamDesc + i] = argDescs[i]->ByVal();
  }

  va_end(argptr);

  std::auto_ptr<Oper> xloArgNames(new Oper(argNames));
  opers[Arg_ArgNames] = xloArgNames->ByVal();


  // do register
  XL::Oper xloRC;
  Excel4v(xlfRegister, XL::ByRef(xloRC), Arg_FirstParamDesc + n,
          const_cast<XLOPER **>(&opers[0]));

  if ( xloRC.IsError() )
  {
    throw EXCEPTION_MSG
          (
            xloRC.GetErrorCode(),
            String::Printf("Failed to register %s().", name)
          );
  }
}

// ============================================================================
// global functions implementation
// ============================================================================

void ReportException()
{
  Error(AddIn::GetExceptionMessage());
}

XLOPER *ReturnException()
{
  return Oper::Return(AddIn::GetExceptionMessage());
}

// ============================================================================
// standard XLL entry functions implementation
// ============================================================================

XLLEXPORT int xlAutoOpen()
{
  XLL_TRACE("xlAutoOpen() called");

  try
  {
    // we may be called before xlAddInManagerInfo when the add-in is being
    // loaded by Excel on startup so make sure we're initialized
    AddInManager::Initialize();

    // register all exported functions
    AddIn::Get().RegisterAllFunctions();

    return TRUE;
  }
  catch ( ... )
  {
    ReportException();

    return FALSE;
  }
}

XLLEXPORT XLOPER *xlAddInManagerInfo(Oper *xloAction)
{
  XLL_TRACE("xloAction(%s) called", xloAction->Dump().c_str());

  try
  {
    // we may be called before xlAutoOpen when the add-in is added from the
    // Excel "Tools|Add ins" menu and not loaded on startup
    AddInManager::Initialize();

    // we need to return the descriptive string if we're passed 1
    if ( *xloAction != 1 )
      return Oper::ErrValue();

    return Oper::Return(AddIn::Get().GetName());
  }
  catch ( ... )
  {
    ReportException();

    return NULL;
  }
}

XLLEXPORT void xlAutoFree(XL::Oper *xlo)
{
  try
  {
    delete xlo;
  }
  catch ( ... )
  {
    ReportException();
  }
}

XLLEXPORT int xlAutoClose()
{
  XLL_TRACE("xlAutoClose() called");

  try
  {
    AddIn::Get().OnClose();
  }
  catch ( ... )
  {
    ReportException();
  }

  return TRUE;
}

// ============================================================================
// DLL entry point
// ============================================================================

BOOL AddIn::DllMain(HMODULE hModule, DWORD dwReason, void* /* lpReserved */)
{
  using namespace ito33;
  using namespace ito33::XL;

  static log::DebugSink *s_sink = NULL;

  switch ( dwReason )
  {
    case DLL_PROCESS_ATTACH:
      {
        // store the DLL handle: we may need it for GetModuleFileName() &c
        Win32::Module::SetHandle(hModule);

        s_sink = new log::DebugSink;
        s_sink->SubscribeAll();

        XLL_TRACE("XLL \"%s\" loaded", Win32::Module::GetFileName().c_str());
      }
      break;

    case DLL_PROCESS_DETACH:
      {
        // now we surely won't be used any more so it's safe to do it
        AddInManager::Cleanup();

        XLL_TRACE("XLL unloaded");

        delete s_sink;
        s_sink = NULL;
      }
      break;
  }

  return TRUE;
}

} // namespace XL

} // namespace ito33

