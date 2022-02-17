/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/volatilitycallback_impl.cpp
// Purpose:     implementation of VolatilityCallbackImpl
// Author:      Vadim Zeitlin
// Created:     2004-06-29
// RCS-ID:      $Id: volatilitycallback_impl.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "com/volatility_impl.h"
#include "com/volatilitycallback_impl.h"
#include "com/callback_impl.h"

#include "ito33/com/c2a.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// classes
// ----------------------------------------------------------------------------

// COMVolatilityCallback implements VolatilityCallback using the provided
// ICallback pointer
class COMVolatilityCallback : public ihg::VolatilityCallBackBase
{
public:
  COMVolatilityCallback(const COM::Ptr<ICallback>& pCallback,
                        double dShift = 0.)
    : ihg::VolatilityCallBackBase(dShift),
      m_pCallback(pCallback)
  {
  }

  // only for calling ICallback methods directly, don't use it otherwise
  ICallback *GetCallback() const { return m_pCallback.Get(); }

protected:
  virtual VolatilityCallBackBase *Clone(double dShift) const
  {
    return new COMVolatilityCallback(m_pCallback, dShift);
  }

  virtual void DoGetValues(double dTime,
                           const double *pdS,
                           double *pdVols,
                           size_t nNbS) const
  {
    return GetValuesFromCallback(m_pCallback.Get(), dTime, pdS, pdVols, nNbS);
  }


  COM::Ptr<ICallback> m_pCallback;
};

// ============================================================================
// VolatilityCallbackImpl implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctor
// ----------------------------------------------------------------------------

VolatilityCallbackImpl::VolatilityCallbackImpl(ICallback *pCallback)
{
  if ( !pCallback )
  {
    throw;
  }

  COM::Ptr<ICallback> ptr;
  ptr.Assign(pCallback);

  m_pImpl = make_ptr( new COMVolatilityCallback(ptr) );
}

// ----------------------------------------------------------------------------
// interface methods
// ----------------------------------------------------------------------------

STDMETHODIMP
VolatilityCallbackImpl::GetValues(double dTime,
                                  SAFEARRAY(double) * aSpots,
                                  SAFEARRAY(double) * rc)
{
  return static_pointer_cast<COMVolatilityCallback>(m_pImpl)->
             GetCallback()->GetValues(dTime, aSpots, rc);
}

/*
STDMETHODIMP
VolatilityCallbackImpl::Perturb(double dShift, IVolatility **rc)
{
  try
  {
      *rc = C2A::COM::Translate<IVolatility *>::To(m_pImpl->Perturb(dShift));

      return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
*/

