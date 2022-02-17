/////////////////////////////////////////////////////////////////////////////
// Name:        com/typelib.cpp
// Purpose:     implementation of TypeLibrary
// Author:      Vadim Zeitlin
// Created:     25.03.03
// RCS-ID:      $Id: typelib.cpp,v 1.6 2004/11/13 10:29:14 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/string.h"

#include "ito33/win32/module.h"

#include "ito33/com/exception.h"
#include "ito33/com/typelib.h"

using namespace ito33::COM;


TypeLibrary::TypeLibrary()
{
  ITypeLib *pTypeLib;

  HRESULT hr = ::LoadTypeLib
                 (
                  String::MB2WC(Win32::Module::GetFileName()),
                  &pTypeLib
                 );

  if ( FAILED(hr) )
    throw COM_EXCEPTION("LoadTypeLib", hr);

  // for some reason using a simple assignment doesn't work in VC++ 6 here
  Ptr<ITypeLib>::operator=(pTypeLib);
}

void TypeLibrary::Register()
{
  HRESULT hr = ::RegisterTypeLib
                 (
                  Get(),
                  String::MB2WC(Win32::Module::GetFileName()),
                  NULL /* help dir */
                 );

  if ( FAILED(hr) )
    throw COM_EXCEPTION("RegisterTypeLib", hr);
}

bool TypeLibrary::Unregister()
{
  TLIBATTR *pAttrs = NULL;
  HRESULT hr = Get()->GetLibAttr(&pAttrs);
  if ( FAILED(hr) || !pAttrs )
    return false;

  hr = ::UnRegisterTypeLib
         (
            pAttrs->guid,
            pAttrs->wMajorVerNum,
            pAttrs->wMinorVerNum,
            pAttrs->lcid,
            pAttrs->syskind
         );

  Get()->ReleaseTLibAttr(pAttrs);

  return SUCCEEDED(hr);
}

