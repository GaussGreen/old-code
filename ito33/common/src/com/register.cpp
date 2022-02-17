/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/register.cpp
// Purpose:     implements COM servers registration functions
// Author:      Vadim Zeitlin
// Created:     19.01.03
// RCS-ID:      $Id: register.cpp,v 1.10 2004/11/13 10:29:32 zeitlin Exp $
// Copyright:   (c) 2002-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"

#include "ito33/win32/regkey.h"
#include "ito33/win32/module.h"

#include "ito33/com/exception.h"
#include "ito33/com/typelib.h"
#include "ito33/com/coclass.h"
#include "ito33/com/register.h"

using ito33::Win32::RegKey;

// ============================================================================
// RegKey implementation
// ============================================================================

namespace ito33
{

namespace COM
{

// ----------------------------------------------------------------------------
// COM servers [un]registration
// ----------------------------------------------------------------------------

void
RegisterCOMServer(const std::string& clsidWithoutBraces,
                  const std::string& progid,
                  int version,
                  const std::string& description,
                  ServerKind kind)
{
  // coclass registration
  // --------------------

  const std::string clsid = '{' + clsidWithoutBraces + '}';

  const std::string
    progidVer = String::Printf("%s.%d", progid.c_str(), version);

  // version independent part
  RegKey key(HKEY_CLASSES_ROOT, progid);
  key = description;
  RegKey(key, "Clsid") = clsid;
  RegKey(key, "CurVer") = progidVer;

  // version dependent part
  RegKey keyVer(HKEY_CLASSES_ROOT, progidVer);
  keyVer = description;
  RegKey(keyVer, "Clsid") = clsid;

  // CLSID
  const std::string fullpath = Win32::Module::GetFileName();

  RegKey keyClsid(HKEY_CLASSES_ROOT, "CLSID\\" + clsid);
  RegKey keyServer(keyClsid, kind == Server_InProcess ? "InprocServer32"
                                                      : "LocalServer32");
  keyServer = fullpath;
  keyServer.Set("ThreadingModel", "Apartment");
  RegKey(keyClsid, "ProgID") = progidVer;
  RegKey(keyClsid, "VersionIndependentProgID") = progid;

  // interface registration (TODO)

  // typelib registration
  // --------------------

  // notice that we assume here that the type library is embedded into the
  // EXE or DLL itself as a resource, not in separate .TLB file
  COM::TypeLibrary().Register();
}

void
UnregisterCOMServer(const std::string& clsidWithoutBraces,
                    const std::string& progid,
                    int version)
{
  // undo the work of RegisterCOMServer() but be ready to handle the keys
  // which don't exist any longer

  // version independent part
  if ( RegKey::Exists(HKEY_CLASSES_ROOT, progid) )
  {
    RegKey key(HKEY_CLASSES_ROOT, progid);
    RegKey::Delete(key, "Clsid");
    RegKey::Delete(key, "CurVer");
    RegKey::Delete(HKEY_CLASSES_ROOT, progid);
  }

  // version dependent part
  const std::string
    progidVer = String::Printf("%s.%d", progid.c_str(), version);
  if ( RegKey::Exists(HKEY_CLASSES_ROOT, progidVer) )
  {
    RegKey keyVer(HKEY_CLASSES_ROOT, progidVer);
    RegKey::Delete(keyVer, "Clsid");
    RegKey::Delete(HKEY_CLASSES_ROOT, progidVer);
  }

  const std::string clsid = "CLSID\\{" + clsidWithoutBraces + '}';
  if ( RegKey::Exists(HKEY_CLASSES_ROOT, clsid) )
  {
    RegKey keyClsid(HKEY_CLASSES_ROOT, clsid);
    if ( RegKey::Exists(keyClsid, "InprocServer32") )
      RegKey::Delete(keyClsid, "InprocServer32");
    if ( RegKey::Exists(keyClsid, "LocalServer32") )
      RegKey::Delete(keyClsid, "LocalServer32");
    RegKey::Delete(keyClsid, "ProgID");
    RegKey::Delete(keyClsid, "VersionIndependentProgID");
    RegKey::Delete(HKEY_CLASSES_ROOT, clsid);
  }

  COM::TypeLibrary().Unregister();
}

// ----------------------------------------------------------------------------
// [un]registration of all coclasses
// ----------------------------------------------------------------------------

void RegisterAll(ServerKind kind)
{
  // register all coclasses
  for ( const COM::CoClassInfo *coclass = COM::CoClassInfo::GetFirst();
        coclass;
        coclass = coclass->GetNext() )
  {
    RegisterCOMServer
    (
     COM::Uuid(coclass->GetClsid()).AsString(),
     coclass->GetProgid(),
     coclass->GetVersion(),
     "",  // no description for now (TODO)
     kind
    );
  }
}

bool UnregisterAll()
{
  try
  {
    // unregister all coclasses
    for ( const COM::CoClassInfo *coclass = COM::CoClassInfo::GetFirst();
          coclass;
          coclass = coclass->GetNext() )
    {
      UnregisterCOMServer
      (
       COM::Uuid(coclass->GetClsid()).AsString(),
       coclass->GetProgid(),
       coclass->GetVersion()
      );
    }

    return true;
  }
  catch ( ... )
  {
    // we don't give any error messages here, it's probably not that
    // important that we failed to unregister (if this ever happens, that
    // is...)
    return false;
  }
}

} // namespace Win32

} // namespace ito33
