/////////////////////////////////////////////////////////////////////////////
// Name:        src/unix/dynlib.cpp
// Purpose:     implementation of DynamicLibrary class for Unix with dlopen()
// Author:      Vadim Zeitlin
// Created:     26.07.03
// RCS-ID:      $Id: dynlib.cpp,v 1.9 2006/02/24 14:45:33 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/error.h"

#include "ito33/errnoexception.h"
#include "ito33/gettext.h"
#include "ito33/dynlib.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// DLException: dlxxx() functions provide an error message which we want to
//              show to the user as part of this exception
// ----------------------------------------------------------------------------

class DLException : public ErrnoException
{
public:
  DLException(const char *apiname,
        const char *filename,
        size_t line,
        const char *function);

  // override the base class function to show m_strDlError too
  virtual std::string GetFullMessage() const;

private:
  // the error message from dlerror()
  std::string m_strDlError;
};

// ============================================================================
// DLException implementation
// ============================================================================

DLException::DLException(const char *apiname,
            const char *filename,
            size_t line,
            const char *function)
     : ErrnoException(apiname, 0, filename, line, function),
      m_strDlError(dlerror())
{
}

std::string DLException::GetFullMessage() const
{
  // we short circuit the base class version because we can't really append
  // to the standard message -- we have to replace it
  std::string msg;
  msg = String::Printf(TRANS("A dynamic library function %s() failed "),
            GetCFunction());
  msg += GetLocation();
  msg += '\n';

  msg += String::Printf(TRANS("with error code %d (%s)"),
             GetErrno(), m_strDlError.c_str());

  return msg;
}

// ============================================================================
// DynamicLibrary implementation
// ============================================================================

// ----------------------------------------------------------------------------
// loading
// ----------------------------------------------------------------------------

bool DynamicLibrary::Load(const char *libname, int flags)
{
  // first unload any library we currently have loaded
  Unload();

  // next adjust the library name
  m_name = libname;
  if ( !(flags & Load_Verbatim) )
  {
    static const char *DLL_SUFFIX = ".so";

    // add the suffix if it's not already present
    if ( m_name.length() < 3 ||
        String::Stricmp(std::string(m_name.end() - 3,
                      m_name.end()).c_str(),
                DLL_SUFFIX) != 0 )
    {
      m_name += DLL_SUFFIX;
    }

    // also add the "lib" prefix if not already there as all Unix
    // libraries start with it
    std::string prefix("lib");
    if ( String::Strnicmp(m_name, prefix, prefix.length()) != 0 )
    {
      m_name = prefix + m_name;
    }
  }

  // also adjust the flags
  int flagsDL = RTLD_NOW; // default as Load_Now == 0
  if ( flags & Load_Lazy )
   flagsDL = RTLD_LAZY;   // can't be combined with RTLD_NOW. don't use |=
  if ( flags & Load_Global )
   flagsDL |= RTLD_GLOBAL;

  // finally do load the library
  m_handle = dlopen(m_name.c_str(), flagsDL);

  return m_handle != 0;
}

// ----------------------------------------------------------------------------
// unloading
// ----------------------------------------------------------------------------

/* static */
bool DynamicLibrary::Unload(DllHandle_t handle)
{
  if ( IsOk(handle) )
  {
    if ( !dlclose(handle) )
    {
      // TODO: log the error somewhere
      return false;
    }
  }

  return true;
}

// ----------------------------------------------------------------------------
// symbol access
// ----------------------------------------------------------------------------

void *DynamicLibrary::HasSymbol(const char *name) const
{
  CHECK( IsOk(), NULL, "Library not loaded, cannot resolve symbols" );
  CHECK( name, NULL, "DynamicLibrary: NULL symbol name" );

  return dlsym(m_handle, name);
}

// ----------------------------------------------------------------------------
// error handling
// ----------------------------------------------------------------------------

std::string DynamicLibrary::GetFullPath() const
{
  // TODO: for Linux we can grep for the library name in /proc/self/maps,
  //       but there is no portable way to find shared object path in general
  return m_name;
}

void DynamicLibrary::Throw(const char *funcname)
{
  throw DLException(funcname, __FILE__, __LINE__, __FUNCTION__);
}

void DynamicLibrary::ThrowOnLoad()
{
  Throw("dlopen");
}

void DynamicLibrary::ThrowOnGetSymbol()
{
  Throw("dlsym");
}


