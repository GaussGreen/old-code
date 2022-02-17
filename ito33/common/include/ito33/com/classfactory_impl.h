/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/classfactory_impl.h
// Purpose:     provides a trivial implementation of IClassFactory interface
// Author:      Vadim Zeitlin
// Created:     08.01.03
// RCS-ID:      $Id: classfactory_impl.h,v 1.12 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/classfactory_impl.h
    @brief Provides a simple IClassFactory implementation.

    This implementation can only create the objects of a single type which must
    be specified as the template parameter.
 */

#ifndef _ITO33_COM_CLASSFACTORY_H_
#define _ITO33_COM_CLASSFACTORY_H_

#include "ito33/com/ptr.h"
#include "ito33/com/traits.h"
#include "ito33/com/unknown_impl.h"

DEFINE_COM_TRAITS(IClassFactory, IUnknown);

namespace ito33
{

namespace COM
{


/**
    Template class factory implementation class.

    This is the simplest possible IClassFactory implementation: it doesn't
    support aggregation and can create the objects of a single type only. The
    type of the objects to create is the template parameter. This class must
    have a default ctor and inherit from IUnknown, i.e. be a COM interface
    implementation.
 */
template <class T>
class ImplementClassFactory : public ImplementUnknown<IClassFactory>
{
public:
  /// save the name of the class we create so that it can be accessed later
  typedef T ClassName;

  /// @name IClassFactory methods
  //@{

  /// creates an instance of the object with the given IID
  STDMETHODIMP CreateInstance(IUnknown* pUnknownOuter,
                const IID& iid,
                void** ppObj);

  /// locks or unlocks the server
  STDMETHODIMP LockServer(BOOL bLock);

  //@}
};

// ----------------------------------------------------------------------------
// template implementation
// ----------------------------------------------------------------------------

template <class T>
STDMETHODIMP
ImplementClassFactory<T>::CreateInstance(IUnknown* pUnknownOuter,
                    const IID& iid,
                    void** ppObj)
{
  CHECK( ppObj, E_UNEXPECTED, "NULL ppObj in IClassFactory::CreateInstance" );

  *ppObj = NULL;

  if ( pUnknownOuter )
  {
    // this implementation doesn't support COM aggregation
    return CLASS_E_NOAGGREGATION;
  }

  // forward the job of deciding whether we support this interface or not to
  // the real object
  try
  {
    COM::Ptr<ClassName> pObject(new ClassName);

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

template <class T>
inline
STDMETHODIMP ImplementClassFactory<T>::LockServer(BOOL bLock)
{
  Server *server = Server::Get();
  if ( !server )
    return E_FAIL;

  if ( bLock )
    server->Lock();
  else
    server->Unlock();

  return S_OK;
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_CLASSFACTORY_H_

