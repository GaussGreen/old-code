/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/callback_impl.cpp
// Purpose:     helper for calling ICallback::GetValues() from C++
// Author:      Vadim Zeitlin
// Created:     2004-06-29
// RCS-ID:      $Id: callback_impl.cpp,v 1.2 2004/10/04 18:04:07 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/com/safearray.h"

// MIDL-generated header defining ICallback
#include "callback.h"

using namespace ito33;

void
GetValuesFromCallback(ICallback *pCallback,
                      double dTime,
                      const double *pdS,
                      double *pdValues,
                      size_t nNbS)
{
  CHECK_VOID( pCallback, "NULL pCallback parameter");

  COM::SafeArrayPtr paSpots(COM::ConvertToSafeArray(nNbS, pdS));
  SAFEARRAY *paValues;
  HRESULT hr = pCallback->GetValues(dTime, paSpots, &paValues);
  if ( FAILED(hr) )
  {
    throw COM_EXCEPTION("ICallback::GetValues", hr);
  }

  // make sure we destroy the returned SAFEARRAY
  COM::SafeArrayPtr ptrValues(paValues);

  COM::SafeArray<double> aValues(paValues);
  COM::SafeArrayAccessor<double> values(aValues);
  for ( size_t n = 0; n < nNbS; n++ )
  {
    *pdValues++ = values[n];
  }
}
