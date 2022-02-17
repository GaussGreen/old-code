/////////////////////////////////////////////////////////////////////////////
// Name:        src/win32/dynlib.cpp
// Purpose:     implementation of DynamicLibrary class
// Author:      Vadim Zeitlin
// Created:     21.03.03
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

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/win32/exception.h"
#include "ito33/dynlib.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// local classes
// ----------------------------------------------------------------------------

// used in DynamicLibrary::Load() but can't be declared there as it is used as
// template parameter and so must have external linkage <sigh>
struct IsBackslashPredicate
{
  inline bool operator()(char ch) const { return ch == '\\'; }
};

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
    static const char *DLL_SUFFIX = ".dll";

    // add the suffix if it's not already present
    if ( m_name.length() < 4 ||
        String::Stricmp(m_name.c_str() + m_name.length() - 4,
                DLL_SUFFIX) != 0 )
    {
      m_name += DLL_SUFFIX;
    }
  }

  // MSDN says that the slashes are not allowed in the path passed to it so
  // replace them with backslashes
  std::replace_if(m_name.begin(), m_name.end(),
          IsBackslashPredicate(), '\\');

  // finally do load the library
  m_handle = LoadLibrary(m_name.c_str());

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
    if ( !FreeLibrary(handle) )
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

  return GetProcAddress(m_handle, name);
}

// ----------------------------------------------------------------------------
// error handling
// ----------------------------------------------------------------------------

std::string DynamicLibrary::GetFullPath() const
{
  char path[MAX_PATH];
  if ( !::GetModuleFileName(m_handle, path, SIZEOF(path)) )
  {
    // TODO: log the error somewhere
    return m_name;
  }

  return path;
}

void DynamicLibrary::Throw(const char *funcname)
{
  throw WIN32_EXCEPTION(funcname);
}

void DynamicLibrary::ThrowOnLoad()
{
  Throw("LoadLibrary");
}

void DynamicLibrary::ThrowOnGetSymbol()
{
  Throw("GetProcAddress");
}

