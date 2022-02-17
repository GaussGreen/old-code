/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/comdll.cpp
// Purpose:     standard COM DLL functions
// Author:      Vadim Zeitlin
// Created:     2006-05-15 (extracted from src/com/dllmain_impl.cpp)
// RCS-ID:      $Id: comdll.cpp,v 1.1 2006/05/15 15:01:23 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file src/com/comdll.cpp
    @brief Implementation of standard COM DLL functions.

    COM in-process server must define several global functions which are more
    or less the same for all COM servers and so they are implemented here once
    and for all.

    The customization is achieved by building the linked list of CoClassInfo
    structures in the main program using DEFINE_COM_COCLASS macro.
 */

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"
#include "ito33/gettext.h"

#include "ito33/win32/regkey.h"

#include "ito33/com/ptr.h"
#include "ito33/com/uuid.h"
#include "ito33/com/coclass.h"
#include "ito33/com/unknown_impl.h"
#include "ito33/com/register.h"

using namespace ito33;

// ============================================================================
// implementation of standard COM DLL functions
// ============================================================================

// ----------------------------------------------------------------------------
// this function creates a new COM object supporting the specified CLSID
// ----------------------------------------------------------------------------

STDAPI
DllGetClassObject(REFCLSID clsid, REFIID iid, void** ppObj)
{
  try
  {
    // find whether we support this CLSID
    const COM::CoClassInfo *coclass = COM::CoClassInfo::Find(clsid);

    if ( !coclass )
    {
      // no coclass for this CLSID
      return CLASS_E_CLASSNOTAVAILABLE;
    }

    // now ask the corresponding object about the interface requested
    COM::Ptr<IUnknown> pObject(coclass->CreateInstance());

    return pObject->QueryInterface(iid, ppObj);
  }
  catch ( const std::bad_alloc& )
  {
    return E_OUTOFMEMORY;
  }
  catch ( ... )
  {
    return E_UNEXPECTED;
  }
}

// ----------------------------------------------------------------------------
// should return S_OK only when the DLL can be safely unloaded
// ----------------------------------------------------------------------------

STDAPI
DllCanUnloadNow()
{
  // we can be safely unloaded only when there are no more objects created by
  // this DLL alive
  return COM::GetNumberOfActiveObjects() ? S_FALSE : S_OK;
}

// ----------------------------------------------------------------------------
// this function should create all the registry entries needed for a COM server
// ----------------------------------------------------------------------------

STDAPI
DllRegisterServer()
{
  try
  {
    COM::RegisterAll(COM::Server_InProcess);

    return S_OK;
  }
  catch ( const std::bad_alloc& )
  {
    return E_OUTOFMEMORY;
  }
  catch ( const ito33::Exception& e )
  {
    ::MessageBox
      (
       NULL,
       e.GetFullMessage().c_str(),
       TRANS("COM Server Registration Error"),
       MB_ICONERROR | MB_OK
      );

    return E_FAIL;
  }
  catch ( ... )
  {
    return E_UNEXPECTED;
  }
}

// ----------------------------------------------------------------------------
// this function should delete all entries created by DllRegisterServer()
// ----------------------------------------------------------------------------

STDAPI
DllUnregisterServer()
{
    return COM::UnregisterAll() ? S_OK : E_UNEXPECTED;
}

