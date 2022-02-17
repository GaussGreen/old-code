///////////////////////////////////////////////////////////////////////////////
// Name:        hellocom/number_impl.h
// Purpose:     declaration of class imlpementing INumber interface
// Author:      Vadim Zeitlin
// Created:     22.01.03
// RCS-ID:      $Id: number_impl.h,v 1.9 2006/01/03 17:13:54 zhang Exp $
// Copyright:   (c) 2002-2003 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

#ifndef _NUMBER_IMPL_H_
#define _NUMBER_IMPL_H_

#include "hello.h"

#include "ito33/com/traits.h"
#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

DEFINE_COM_TRAITS(INumber, IUnknown);
DEFINE_COM_TRAITS(INumbers, IUnknown);

// normally there should be no "using" in header but this is just a test...
using namespace ito33;

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

class NumberImpl : public COM::ImplementCoClass<INumber, NumberImpl>
{
public:
  NumberImpl() : m_value(0) { }

  STDMETHODIMP put_Value(int value) { m_value = value; return S_OK; }
  STDMETHODIMP get_Value(int *value) { *value = m_value; return S_OK; }

  STDMETHODIMP get_Sign(ISign **ppSign);

private:
  int m_value;

  NO_COPY_CLASS(NumberImpl);
};

class NumbersImpl : public COM::ImplementCustomUnknown
                <
                  INumbers,
                  TYPE_LIST_2(ISupportErrorInfo, IDispatch)
                >
{
public:
  NumbersImpl();

  STDMETHODIMP get_Count(long *pCount);
  STDMETHODIMP get_Item(long n, long *pNumber);
  STDMETHODIMP get__NewEnum(IUnknown **ppUnknown);

private:
  typedef std::vector<long> Points;
  Points m_points;
};

#endif // _NUMBER_IMPL_H_

