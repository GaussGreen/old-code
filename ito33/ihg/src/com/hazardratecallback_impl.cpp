/////////////////////////////////////////////////////////////////////////////
// Name:        src/com/hazardratecallback_impl.cpp
// Purpose:     implementation of HazardRateCallbackImpl
// Author:      Vadim Zeitlin
// Created:     2004-06-29
// RCS-ID:      $Id: hazardratecallback_impl.cpp,v 1.6 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "com/hazardrate_impl.h"
#include "com/hazardratecallback_impl.h"
#include "com/callback_impl.h"

#include "ito33/com/c2a.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// classes
// ----------------------------------------------------------------------------

// COMHazardRateCallback implements HazardRateCallback using the provided
// ICallback pointer
class COMHazardRateCallback : public ihg::HazardRate
{
public:
  COMHazardRateCallback(const COM::Ptr<ICallback>& pCallback)
    : m_pCallback(pCallback)
  {
  }

  virtual void GetHazardRates(double dTime,
                              const double *pdS,
                              double *pdValues,
                              size_t nNbS) const
  {
    return GetValuesFromCallback(m_pCallback.Get(), dTime, pdS, pdValues, nNbS);
  }

  virtual void Dump(ito33::XML::Tag& /* tagParent */) const
  {
  }

  void Visit(ito33::ihg::HazardRateVisitor &) const
  {
  }

  // only for calling ICallback methods directly, don't use it otherwise
  ICallback *GetCallback() const { return m_pCallback.Get(); }

private:
  COM::Ptr<ICallback> m_pCallback;
};

// ============================================================================
// HazardRateCallbackImpl implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctor
// ----------------------------------------------------------------------------

HazardRateCallbackImpl::HazardRateCallbackImpl(ICallback *pCallback)
{
  if ( !pCallback )
  {
    throw;
  }

  COM::Ptr<ICallback> ptr;
  ptr.Assign(pCallback);

  m_pImpl = make_ptr( new COMHazardRateCallback(ptr) );
}

// ----------------------------------------------------------------------------
// interface methods
// ----------------------------------------------------------------------------

STDMETHODIMP
HazardRateCallbackImpl::GetValues(double dTime,
                                  SAFEARRAY(double) * aSpots,
                                  SAFEARRAY(double) * rc)
{
    return static_pointer_cast<COMHazardRateCallback>(m_pImpl)->
              GetCallback()->GetValues(dTime, aSpots, rc);
}

