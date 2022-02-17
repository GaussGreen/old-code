///////////////////////////////////////////////////////////////////////////////
// Name:        hellocom/number_impl.cpp
// Purpose:     INumber and INumbers implementation
// Author:      Vadim Zeitlin
// Created:     28.03.03
// RCS-ID:      $Id: number_impl.cpp,v 1.5 2004/10/05 09:13:50 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

#include "number_impl.h"
#include "sign_impl.h"

#include "ito33/com/enum.h"

#undef COM_APPNAME
#define COM_APPNAME "Ito33HelloLib"
DEFINE_COM_COCLASS(Number, 1);

// ----------------------------------------------------------------------------
// NumberImpl
// ----------------------------------------------------------------------------

STDMETHODIMP NumberImpl::get_Sign(ISign **ppSign)
{
  *ppSign = new SignImpl(m_value);

  return S_OK;
}

// ----------------------------------------------------------------------------
// NumbersImpl
// ----------------------------------------------------------------------------

NumbersImpl::NumbersImpl()
{
  m_points.push_back(17);
  m_points.push_back(9);
  m_points.push_back(51);
}

STDMETHODIMP NumbersImpl::get_Count(long *pCount)
{
  *pCount = m_points.size();

  return S_OK;
}

STDMETHODIMP NumbersImpl::get_Item(long n, long *pNumber)
{
  if ( n < 0 || (size_t)n >= m_points.size() )
    return E_INVALIDARG;

  *pNumber = m_points[n];

  return S_OK;
}

STDMETHODIMP NumbersImpl::get__NewEnum(IUnknown **retval)
{
  COM::StdEnum<Points> *enumPoints = new COM::StdEnum<Points>(m_points);
  *retval = enumPoints;

  return S_OK;
}

