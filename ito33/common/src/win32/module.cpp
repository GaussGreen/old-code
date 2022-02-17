/////////////////////////////////////////////////////////////////////////////
// Name:        win32/module.cpp
// Purpose:     implementation of global (module-level) Win32 functions
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: module.cpp,v 1.5 2004/10/05 09:13:47 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"

#include "ito33/win32/exception.h"
#include "ito33/win32/module.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// global data
// ----------------------------------------------------------------------------

// this variable is not protected by critical section because normally it is
// set only once, in the very beginning of the program (and hence before any
// threads other than the main one were created), and is only read afterwards
static HMODULE gs_hModule = NULL;

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// Module::Set/GetHandle()
// ----------------------------------------------------------------------------

inline HMODULE ito33::Win32::Module::GetHandle()
{
  ASSERT_MSG( gs_hModule, "SetHandle() must have been called first!" );

  return gs_hModule;
}

void ito33::Win32::Module::SetHandle(HMODULE hModule)
{
  ASSERT_MSG( !gs_hModule, "SetHandle() should only be called once" );
  ASSERT_MSG( hModule, "HMODULE in SetHandle() must not be NULL" );

  gs_hModule = hModule;
}

// ----------------------------------------------------------------------------
// Module::GetFileName()
// ----------------------------------------------------------------------------

std::string ito33::Win32::Module::GetFileName()
{
  // this should be enough -- and is the best we can do as
  // GetModuleFileName() doesn't return the size of the needed buffer
  char buf[4096];
  if ( !::GetModuleFileName(GetHandle(), buf, SIZEOF(buf)) )
  {
    throw WIN32_EXCEPTION("GetModuleFileName");
  }

  return buf;
}

