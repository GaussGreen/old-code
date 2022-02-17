/////////////////////////////////////////////////////////////////////////////
// Name:        hellocom/sign_impl.h
// Purpose:     declaration of class imlpementing ISign interface
// Author:      Vadim Zeitlin
// Created:     22.01.03
// RCS-ID:      $Id: sign_impl.h,v 1.4 2006/01/03 17:13:54 zhang Exp $
// Copyright:   (c) 2002-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _SIGN_IMPL_H_
#define _SIGN_IMPL_H_

#include "hello.h"

#include "ito33/com/unknown_impl.h"

DEFINE_COM_TRAITS(ISign, IUnknown);

// normally there should be no "using" in header but this is just a test...
using namespace ito33;

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

class SignImpl : public COM::ImplementUnknown<ISign>
{
public:
  SignImpl(int sign) : m_sign(sign > 0 ? 0 : -1) { }

  STDMETHODIMP get_IsPositive(BOOL *sign) { *sign = m_sign; return S_OK; }

private:
  int m_sign;
};

#endif // _SIGN_IMPL_H_

