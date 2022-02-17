/////////////////////////////////////////////////////////////////////////////
// Name:        hellocom/hello_impl.h
// Purpose:     declaration of class imlpementing IHello interface
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: hello_impl.h,v 1.15 2005/09/22 19:13:43 zeitlin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _HELLO_IMPL_H_
#define _HELLO_IMPL_H_

#include "hello.h"

#include "ito33/date.h"

#include "ito33/com/ptr.h"
#include "ito33/com/traits.h"
#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"

DEFINE_COM_TRAITS(IHello, IUnknown);

// normally there should be no "using" in header but this is just a test...
using namespace ito33;

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

class Ito33HelloImpl : public COM::ImplementCoClass
               <
                IHello,
                Ito33HelloImpl,
                TYPE_LIST_2(ISupportErrorInfo, IDispatch)
               >
{
public:
  Ito33HelloImpl() { }

  // IHello methods
  STDMETHODIMP SayHello(BSTR message);
  STDMETHODIMP FindAverage(SAFEARRAY **data);
  STDMETHODIMP get_Values(SAFEARRAY **data);
  STDMETHODIMP Fail();

  STDMETHODIMP get_Date(DATE *pDate);
  STDMETHODIMP put_Date(DATE date);
  STDMETHODIMP MaxDate(SAFEARRAY **dates, DATE *pDate);

  STDMETHODIMP putref_NumberRef(INumber *ppNumber);
  STDMETHODIMP put_Number(INumber *pNumber);
  STDMETHODIMP CreateNumber(INumber **ppNumber);
  STDMETHODIMP get_Numbers(INumbers **ppNumbers);

private:
  COM::Ptr<INumber> m_pNumber;
  ito33::Date m_date;

  NO_COPY_CLASS(Ito33HelloImpl);
};

#endif // _HELLO_IMPL_H_

