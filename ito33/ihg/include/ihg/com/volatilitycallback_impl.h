/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/com/volatilitycallback.h
// Purpose:     volatility callbacks for COM
// Author:      Vadim Zeitlin
// Created:     2004-06-29
// RCS-ID:      $Id: volatilitycallback_impl.h,v 1.4 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_COM_VOLATILITYCALLBACK_IMPL_H_
#define _IHG_COM_VOLATILITYCALLBACK_IMPL_H_

#include "ito33/sharedptr.h"

#include "ito33/com/ptr.h"
#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"
#include "ito33/com/unknown_impl.h"

#include "ito33/ihg/volatilitycallback.h"

// include the MIDL-generated headers defining the interfaces
#include "callback.h"
#include "volatility.h"

// ----------------------------------------------------------------------------
// IVolatility implementation using the user-provided ICallback
// ----------------------------------------------------------------------------

class VolatilityCallbackImpl : public ito33::COM::ImplementCustomUnknown
                                      <
                                          IVolatility,
                                          TYPE_LIST_2
                                          (
                                            ISupportErrorInfo,
                                            IDispatch
                                          )
                                      >
{
public:
  // NB: not the same as VolatilityCallBack from ito33/ihg/volatilitycallback.h
  typedef ito33::ihg::VolatilityCallBackBase VolatilityCallback;

  VolatilityCallbackImpl(ICallback *pCallback);

  STDMETHODIMP GetValues(double dTime, SAFEARRAY(double) * aSpots, SAFEARRAY(double) * rc);
  // STDMETHODIMP Perturb(double dShift, IVolatility ** rc);


  // for internal use only
  const ito33::shared_ptr<VolatilityCallback>& GetImpl() const { return m_pImpl; }

private:
  ito33::shared_ptr<VolatilityCallback> m_pImpl;

  VolatilityCallbackImpl(const VolatilityCallbackImpl&);
  VolatilityCallbackImpl& operator=(const VolatilityCallbackImpl&);
};

#endif // _IHG_COM_VOLATILITYCALLBACK_IMPL_H_
