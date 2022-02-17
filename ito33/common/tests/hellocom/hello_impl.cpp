/////////////////////////////////////////////////////////////////////////////
// Name:        hellocom/hello_impl.cpp
// Purpose:     implementation of IHello interface
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: hello_impl.cpp,v 1.16 2005/09/26 15:31:36 zeitlin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "hello_impl.h"
#include "number_impl.h"

#include "ito33/com/coclass.h"
#include "ito33/com/safearray.h"

extern const ito33::Error ITO33_BAD_PARAM, ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::COM;

#undef COM_APPNAME
#define COM_APPNAME "Ito33HelloLib"
DEFINE_COM_COCLASS(Ito33Hello, 1);

// ----------------------------------------------------------------------------
// various IHello methods
// ----------------------------------------------------------------------------

STDMETHODIMP Ito33HelloImpl::SayHello(BSTR message)
{
  wchar_t buf[512];
  wsprintfW(buf, L"Hello, %.100s!", message);

  MessageBoxW(NULL, buf, L"Hello COM Server", MB_OK | MB_ICONINFORMATION);

  return S_OK;
}

STDMETHODIMP Ito33HelloImpl::FindAverage(SAFEARRAY **ppData)
{
  try
  {
    SafeArrayAccessor<int> data(*ppData);

    const size_t size = data.GetCount();
    if ( size == 0 )
      throw ito33::EXCEPTION_MSG(ITO33_BAD_PARAM, "Array can't be empty");

    char buf[512];

    int sum = 0;
    for ( size_t n = 0; n < size; ++n )
    {
      sum += data[n];
    }

    sprintf(buf, "Average of %d elements is %g", size, ((float)sum)/size);

    MessageBoxA(NULL, buf, "Hello COM Server", MB_OK | MB_ICONINFORMATION);

    return S_OK;
  }
  catch ( ito33::Exception& e )
  {
    return SetErrorInfo(e);
  }
  catch ( ... )
  {
    return E_UNEXPECTED;
  }
}

STDMETHODIMP Ito33HelloImpl::Fail()
{
  SetErrorInfo(ito33::EXCEPTION_MSG
            (
              ITO33_UNEXPECTED,
              "This function always fails -- as it just did."
            ));

  return E_FAIL;
}

STDMETHODIMP Ito33HelloImpl::get_Date(DATE *pDate)
{
  SYSTEMTIME st;
  ::ZeroMemory(&st, sizeof(st));
  st.wYear = static_cast<WORD>(m_date.GetYear());
  st.wMonth = static_cast<WORD>(m_date.GetMonth());
  st.wDayOfWeek = static_cast<WORD>(m_date.GetDayOfWeek());
  st.wDay = (WORD)m_date.GetDay();

  if ( !::SystemTimeToVariantTime(&st, pDate) )
    return E_FAIL;

  return S_OK;
}

STDMETHODIMP Ito33HelloImpl::put_Date(DATE date)
{
  SYSTEMTIME st;
  if ( !::VariantTimeToSystemTime(date, &st) )
    return E_FAIL;

  m_date.Set(st.wYear, static_cast<Date::Month>(st.wMonth), st.wDay);

  return S_OK;
}

STDMETHODIMP Ito33HelloImpl::MaxDate(SAFEARRAY **ppDates, DATE *pDate)
{
  try
  {
    SafeArrayAccessor<DATE> pDates(*ppDates);

    const size_t size = pDates.GetCount();
    if ( size == 0 )
      throw ito33::EXCEPTION_MSG(ITO33_BAD_PARAM, "Array can't be empty");

    *pDate = 0;
    for ( size_t n = 0; n < size; ++n )
    {
      if ( pDates[n] > *pDate )
        *pDate = pDates[n];
    }

    return S_OK;
  }
  catch ( ito33::Exception& e )
  {
    return SetErrorInfo(e);
  }
  catch ( ... )
  {
    return E_UNEXPECTED;
  }
}

STDMETHODIMP Ito33HelloImpl::get_Values(SAFEARRAY **data)
{
  try
  {
    SafeArray<int> array(3);

    SafeArrayAccessor<int> ptr(array.Get());
    ptr[0] = 17;
    ptr[1] = 7;
    ptr[2] = 9;

    *data = array.Get();

    return S_OK;
  }
  catch ( ito33::Exception& e )
  {
    return SetErrorInfo(e);
  }
  catch ( ... )
  {
    return E_UNEXPECTED;
  }
}

STDMETHODIMP Ito33HelloImpl::CreateNumber(INumber **ppNumber)
{
  *ppNumber = new NumberImpl;

  return S_OK;
}

STDMETHODIMP Ito33HelloImpl::get_Numbers(INumbers **ppNumbers)
{
  *ppNumbers = new NumbersImpl;

  return S_OK;
}

STDMETHODIMP Ito33HelloImpl::putref_NumberRef(INumber *pNumber)
{
  CHECK( pNumber, E_INVALIDARG, "NULL parameter in Hello::putref_NumberRef" );

  // shouldn't Release() the [in] parameter!
  pNumber->AddRef();
  m_pNumber = pNumber;

  m_pNumber->put_Value(1);

  return S_OK;
}

STDMETHODIMP Ito33HelloImpl::put_Number(INumber *pNumber)
{
  pNumber->AddRef();
  m_pNumber = pNumber;

  m_pNumber->put_Value(2);

  return S_OK;
}

