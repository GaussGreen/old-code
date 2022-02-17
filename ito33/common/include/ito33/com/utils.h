/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/utils.h
// Purpose:     various COM utility classes and functions
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: utils.h,v 1.7 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/utils.h
    @brief Various COM utility classes and functions.
 */

#ifndef _ITO33_COM_UTILS_H_
#define _ITO33_COM_UTILS_H_

#include "ito33/debug.h"

#include "ito33/com/traits.h"
#include "ito33/com/ptr.h"

namespace ito33
{

namespace COM
{

/**
    UuidOf returns the IID (interface id) of the interface of the given
    pointer.

    This function replaces the (unportable) MSVC++ @c uuidof operator and
    returns the UUID (IID, GUID, CLSID or whatever) for the given pointer. For
    this to work, the COM::Traits<> specialization must be provided.

    Note that it is possible to use this function either with a pointer and
    then the template argument will be deduced automatically or by explicitly
    specifying the template specialization to use -- and then you don't have to
    specify the argument at all as it is not used anyhow and has a default
    value.
 */
template <class Iface>
inline
const UUID& UuidOf(const Iface * /* iface */ = NULL)
{
  return ::ito33::COM::Traits<Iface>::Uuid();
}

/**
    Type-safely cast COM pointer to another type.

    @param pUnknown any COM pointer which must not be @c NULL
    @return pointer cast to target IID (to be released by caller) or @c NULL
 */
template <class T>
inline
T *CastTo(IUnknown *pUnknown)
{
  T *pIface;
  if ( FAILED(pUnknown->QueryInterface
                        (
                          Traits<T>::Uuid(),
                          reinterpret_cast<void **>(&pIface)
                        )) )
  {
    pIface = NULL;
  }

  return pIface;
}

/**
    @name Helper functions for implementing COM accessors.
 */
//@{

/// returns the value of the given parameter
template <typename T>
inline
HRESULT ReturnValue(const T& param, T *pResult)
{
    CHECK( pResult, E_POINTER, "NULL output pointer" );

    *pResult = param;

    return S_OK;
}

/// special overload for booleans as they need to be translated to VARIANT_BOOL
inline
HRESULT ReturnValue(const bool& param, VARIANT_BOOL *pResult)
{
    CHECK( pResult, E_POINTER, "NULL output pointer" );

    *pResult = param ? VARIANT_TRUE : VARIANT_FALSE;

    return S_OK;
}

//@}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_UTILS_H_

