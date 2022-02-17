/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/dispatch.cpp
// Purpose:     implementation of IDispatch support
// Author:      Vadim Zeitlin
// Created:     25.03.03
// RCS-ID:      $Id: dispatch.cpp,v 1.8 2004/10/05 09:13:42 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/com/exception.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/typelib.h"

using namespace ito33::COM;

// ----------------------------------------------------------------------------
// SupportIDispatch
// ----------------------------------------------------------------------------

HRESULT
SupportIDispatch::DoQueryInterface(const IID& iid, void **ppObj)
{
  if ( iid != IID_IDispatch )
    return E_NOINTERFACE;

  if ( !ppObj )
    return E_POINTER;

  IDispatch *pDispatch = GetDispatch();
  if ( !pDispatch )
    return E_NOINTERFACE;

  pDispatch->AddRef();
  *ppObj = static_cast<void *>(pDispatch);

  return S_OK;
}

IDispatch *SupportIDispatch::GetDispatch()
{
  if ( !m_pDispatch )
  {
    try
    {
      TypeLibrary typeLib;

      IUnknown *pUnknown;
      HRESULT hr = ::CreateStdDispatch
             (
              m_pThis,          // outer unknown
              m_pThis,          // the object to expose
              typeLib.GetTypeInfo(m_iidThis).Get(),
              &pUnknown
             );

      if ( FAILED(hr) )
      {
        throw COM_EXCEPTION("CreateStdDispatch", hr);   
      }

      hr = pUnknown->QueryInterface
             (
              IID_IDispatch,
              reinterpret_cast<void **>(&m_pDispatch)
             );

      if ( FAILED(hr) )
      {
        throw COM_EXCEPTION("QueryInterface(IID_IDispatch)", hr);
      }
    }
    catch ( ... )
    {
      FAIL( "Failed to create the standard IDispatch implementation" );
    }
  }

  return m_pDispatch;
}

