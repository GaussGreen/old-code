/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/com/hazardratecallback.h
// Purpose:     hazardrate callbacks for COM
// Author:      Vadim Zeitlin
// Created:     2004-06-29
// RCS-ID:      $Id: hazardratecallback_impl.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_COM_HAZARDRATECALLBACK_IMPL_H_
#define _IHG_COM_HAZARDRATECALLBACK_IMPL_H_

#include "ito33/sharedptr.h"

#include "ito33/com/ptr.h"
#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"
#include "ito33/com/unknown_impl.h"

#include "ito33/ihg/hazardrate.h"

// include the MIDL-generated headers defining the interfaces
#include "callback.h"
#include "hazardrate.h"

// ----------------------------------------------------------------------------
// IHazardRate implementation using the user-provided ICallback
// ----------------------------------------------------------------------------

class HazardRateCallbackImpl : public ito33::COM::ImplementCustomUnknown
                                      <
                                          IHazardRate,
                                          TYPE_LIST_2
                                          (
                                            ISupportErrorInfo,
                                            IDispatch
                                          )
                                      >
{
public:
  HazardRateCallbackImpl(ICallback *pCallback);

  STDMETHODIMP GetValues(double dTime, SAFEARRAY(double) * aSpots, SAFEARRAY(double) * rc);

  // for internal use only
  const ito33::shared_ptr<ito33::ihg::HazardRate>& GetImpl() const { return m_pImpl; }

private:
  ito33::shared_ptr<ito33::ihg::HazardRate> m_pImpl;

  HazardRateCallbackImpl(const HazardRateCallbackImpl&);
  HazardRateCallbackImpl& operator=(const HazardRateCallbackImpl&);
};

#endif // _IHG_COM_HAZARDRATECALLBACK_IMPL_H_
