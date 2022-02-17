/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/coinit.h
// Purpose:     COM::CoInitializer class
// Author:      Vadim Zeitlin
// Created:     Jan 20, 2004
// RCS-ID:      $Id: coinit.h,v 1.5 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/coinit.h
    @brief Helper class for COM initialization.

    Some very simple helper classes wrapping low-level COM functions.
 */

#ifndef _ITO33_COM_COINIT_H_
#define _ITO33_COM_COINIT_H_

#include "ito33/win32/winwrap.h"
#include <objbase.h>

#include "ito33/com/exception.h"
#include "ito33/com/ptr.h"

namespace ito33
{

namespace COM
{

/**
  A tiny helper class to initialize and uninitialize COM for the application.

  Before using COM ::CoInitialize() must always be called and
  ::CoUninitialize() has to be called exactly once for each successfull
  CoInitialize() call. This begs for a simple class encapsulating these class
  in its ctor and dtor respectively -- and this is exactly what CoInitializer
  is.

  Usually an object of this class should be created in main()/WinMain() before
  making any COM calls.
 */
class CoInitializer
{
public:
  /**
      Ctor initializes COM and possibly returns whether this was done for the
      first time or not.

      If COM can't be initialized, an exception is thrown.

      @param pFirstTime if this pointer is non @c NULL, true is returned in it
                        if this is the first time CoInitialize() is called in
                        this thread and false otherwise, i.e. if COM had been
                        already (successfully) initialzied in this thread
   */
  CoInitializer(bool *pFirstTime = NULL);

  /**
      Dtor uninitializes COM.
   */
  ~CoInitializer();
};


/**
    Class for registerting the (main) application class object and revoking its
    registration automatically.

    This class registers the class object in its ctor and revokes it in its
    dtor. As all out of process servers must do the former on startup and the
    latter on shutdown, the usual usage of this class is to simply create a
    local object of its type in WinMain().
 */
class ClassObjectRegistrar
{
public:
  /**
      Registers the class object in an out of process server with COM.

      We hold on to the pUnknown passed to us and will release it in dtor.

      @param clsid the CLSID to be registered
      @param pUnknown the object implementing the clsid
   */
  ClassObjectRegistrar(const CLSID& clsid, const COM::Ptr<IUnknown>& pUnknown);

  /**
      Unregisters the class object registered by our ctor.

      After this, the clients can't access the class object any more.

      The dtor also releases the object passed to the ctor which normally
      results in its destruction.
   */
  ~ClassObjectRegistrar();

private:
  // our registration token
  DWORD m_dwToken;
};

// ============================================================================
// inline functions implementation
// ============================================================================

// ----------------------------------------------------------------------------
// CoInitializer
// ----------------------------------------------------------------------------

CoInitializer::CoInitializer(bool *pFirstTime)
{
  HRESULT hr = ::CoInitialize(NULL /* reserved */);
  if ( FAILED(hr) )
  {
    throw COM_EXCEPTION("CoInitialize", hr);
  }

  if ( pFirstTime )
  {
    // S_FALSE means we succeeded but it had been already done before
    ASSERT_MSG( hr == S_OK || hr == S_FALSE,
                  "unexpected CoInitialize() return value" );

    *pFirstTime = hr == S_OK;
  }
}

CoInitializer::~CoInitializer()
{
  // no error return
  ::CoUninitialize();
}

// ----------------------------------------------------------------------------
// ClassObjectRegistrar
// ----------------------------------------------------------------------------

ClassObjectRegistrar::ClassObjectRegistrar(const CLSID& clsid,
                                           const COM::Ptr<IUnknown>& pUnknown)
{
  HRESULT hr = ::CoRegisterClassObject
                 (
                    clsid,
                    pUnknown.Get(),
                    CLSCTX_SERVER,
                    REGCLS_MULTIPLEUSE,
                    &m_dwToken
                 );
  if ( hr != S_OK )
  {
    // we have to release it manually here as CoRegisterClassObject() failed,
    // if it succeeds it does this itself when we revoke the class object later
    pUnknown->Release();

    throw COM_EXCEPTION("CoRegisterClassObject", hr);
  }
}

ClassObjectRegistrar::~ClassObjectRegistrar()
{
#ifndef _NDEBUG
  HRESULT hr =
#endif // !_NDEBUG

  ::CoRevokeClassObject(m_dwToken);

#ifndef _NDEBUG
  if ( hr != S_OK )
  {
    // dtor shouldn't throw so we don't let the exception propagate but
    // we should still log it
    std::string
      msg = COM_EXCEPTION("CoRevokeClassObject", hr).GetErrorMessage();
    ::OutputDebugString((msg + "\r\n").c_str());
  }
#endif // !_NDEBUG
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_COINIT_H_

