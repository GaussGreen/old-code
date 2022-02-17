/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/com/ihg_impl.cpp
// Purpose:     declaration and implementation of IIHG interface
// Author:      Vadim Zeitlin
// Created:     13.12.03
// RCS-ID:      $Id: ihg_impl.cpp,v 1.38 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2002-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/ihg/showversion.h"
#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"
#include "ito33/com/unknown_impl.h"
#include "ito33/com/c2a.h"

#include "com/hazardrate_impl.h"
#include "com/hazardratetimeonly_impl.h"
#include "com/hazardratepower_impl.h"
#include "com/spotcomponent_impl.h"
#include "com/hrspotcomponentpower_impl.h"
#include "com/hazardratecallback_impl.h"
#include "com/hazardratecombo_impl.h"
#include "com/volatility_impl.h"
#include "com/volatilitycallback_impl.h"
#include "com/volatilityflat_impl.h"
#include "com/volatilitypower_impl.h"
#include "com/volatilitytanh_impl.h"
#include "com/volatilitytimeonly_impl.h"
#include "com/volatilitycombo_impl.h"
#include "com/parametrization_hrwithspotcomponentpower_impl.h"
#include "com/parametrization_volpower_impl.h"
#include "com/parametrization_voltanh_impl.h"
#include "com/parametrization_volwithtimecomponent_impl.h"

#include "ihg.h"

using namespace ito33;
using std::string;

DEFINE_COM_IFACE(GlobalIHG, IUnknown);

// ----------------------------------------------------------------------------
// IHGImpl declaration
// ----------------------------------------------------------------------------

class GlobalIHGImpl : public COM::ImplementCoClass
                             <
                               IGlobalIHG,
                               GlobalIHGImpl,
                               TYPE_LIST_2(ISupportErrorInfo, IDispatch)
                             >
{
public:

  GlobalIHGImpl() { }

  STDMETHODIMP get_Version(BSTR *version);

  // IIHG methods
  STDMETHODIMP NewVolatilityFlat(double dVol, IVolatilityFlat **ppVol);

  STDMETHODIMP NewVolatilityPower(double dAlpha, 
                                  double dBeta,
                                  double dS0,
                                  IVolatilityPower **ppVol);

  STDMETHODIMP NewVolatilityTanh(double dLeft, 
                                 double dRight,
                                 double dScale,
                                 double dS0,
                                 IVolatilityTanh **ppVol);

  STDMETHODIMP NewVolatilityTimeOnly(SAFEARRAY(DATE) *dates,
                                     SAFEARRAY(double) *values,
                                     IVolatilityTimeOnly **ppVol);

  STDMETHODIMP NewVolatilityCombo(ISpotComponent *pSpotComponent,
                                  SAFEARRAY(DATE) *dates,
                                  SAFEARRAY(double) *values,
                                  IVolatilityCombo **ppVol);

  STDMETHODIMP NewVolatilityCallback(ICallback *pCallback, IVolatility **ppVol);

#if 0
  STDMETHODIMP NewHazardRateFlat(double dValue, IHazardRateFlat **ppHR);
#endif
  STDMETHODIMP NewHazardRateCallback(ICallback *pCallback, IHazardRate **ppHR);
  STDMETHODIMP NewHazardRateTimeOnly(SAFEARRAY(DATE) *dates,
                                     SAFEARRAY(double) *values,
                                     IHazardRateTimeOnly **ppHR);

  STDMETHODIMP NewHazardRatePower(double dAlpha,
                                  double dBeta,
                                  double dS0,
                                  IHazardRatePower **ppHR);

  STDMETHODIMP NewHRSpotComponentPower(double dBeta,
                                       double dS0,
                                       IHRSpotComponentPower **ppSpotcomponent);

  STDMETHODIMP NewHazardRateCombo(ISpotComponent *pSpotComponent,
                                  SAFEARRAY(DATE) *dates,
                                  SAFEARRAY(double) *values,
                                  IHazardRateCombo **ppHR);

  /*
  STDMETHODIMP NewImpliedCDSSpreads(DATE contractingDate,
                                    DATE firstPaymentDate,
                                    DayCountConvention dcc,
                                    Frequency freq,
                                    double dRecoveryRate,
                                    IImpliedCDSSpreads **ppImpliedCDSSpreads);
  */

  STDMETHODIMP NewParametrizationHRWithSpotComponentPower(IVolatility *pVol,
                                                          IParametrizationHRWithSpotComponentPower **ppParametrization);

  STDMETHODIMP NewParametrizationVolPower(IHazardRate *pHR,
                                          IParametrizationVolPower **ppParametrizationVolPower); 

  STDMETHODIMP NewParametrizationVolTanh(IHazardRate *pHR,
                                         IParametrizationVolTanh **ppParametrizationVolTanh); 

  STDMETHODIMP NewParametrizationVolWithTimeComponent(IHazardRate *pHR,
                                                      IParametrizationVolWithTimeComponent **ppParametrization); 

private:

  GlobalIHGImpl(const GlobalIHGImpl&);
  GlobalIHGImpl& operator=(const GlobalIHGImpl&);
};


// ============================================================================
// IIHG implementation
// ============================================================================

DEFINE_COM_COCLASS_FULL(CLSID_GlobalIHG, "IHG.GlobalIHG", 1, GlobalIHGImpl);

STDMETHODIMP 
GlobalIHGImpl::get_Version(BSTR *version)
{
  try
  {
    *version = C2A::COM::Translate<string>::To( ito33::ihg::ShowVersion() );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

// ----------------------------------------------------------------------------
// Volatility creation
// ----------------------------------------------------------------------------

STDMETHODIMP
GlobalIHGImpl::NewVolatilityFlat(double dVol, IVolatilityFlat **ppVol)
{
  try
  {
    *ppVol = new VolatilityFlatImpl(shared_ptr<ihg::VolatilityFlat>(
                new ihg::VolatilityFlat(dVol)));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
GlobalIHGImpl::NewVolatilityPower(double dAlpha,
                            double dBeta,
                            double dS0,
                            IVolatilityPower **ppVol)
{
  try
  {
    *ppVol = new VolatilityPowerImpl(shared_ptr<ihg::VolatilityPower>(
                new ihg::VolatilityPower(dAlpha, dBeta, dS0)));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalIHGImpl::NewVolatilityTanh(double dLeft, 
                           double dRight,
                           double dScale,
                           double dS0,
                           IVolatilityTanh **ppVol)
{
  try
  {
    *ppVol = new VolatilityTanhImpl(shared_ptr<ihg::VolatilityTanh>(
                new ihg::VolatilityTanh(dLeft, dRight, dScale, dS0)));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
GlobalIHGImpl::NewVolatilityTimeOnly(SAFEARRAY(DATE) *dates,
                                     SAFEARRAY(double) *values,
                                     IVolatilityTimeOnly **ppVol)
{
  try
  {
    *ppVol = new VolatilityTimeOnlyImpl(shared_ptr<ihg::VolatilityTimeOnly>(
               new ihg::VolatilityTimeOnly(C2A::COM::ToVector<Date, DATE>(dates),
                                          C2A::COM::ToVector<double, double>(values))
             ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalIHGImpl::NewVolatilityCombo(ISpotComponent *pSpotComponent,
                                  SAFEARRAY(DATE) *dates,
                                  SAFEARRAY(double) *values,
                                  IVolatilityCombo **ppVol)
{
  try
  {
    *ppVol = new VolatilityComboImpl(shared_ptr<ihg::VolatilityCombo>
                 ( new ihg::VolatilityCombo
                       (
                         C2A::COM::Translate<ISpotComponent *>::From(pSpotComponent),
                         C2A::COM::ToVector<Date, DATE>(dates),
                         C2A::COM::ToVector<double, double>(values)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
GlobalIHGImpl::NewVolatilityCallback(ICallback *pCallback, IVolatility **ppVol)
{
  try
  {
    *ppVol = new VolatilityCallbackImpl(pCallback);

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

// ----------------------------------------------------------------------------
// Hazard rate objects creation
// ----------------------------------------------------------------------------

STDMETHODIMP
GlobalIHGImpl::NewHazardRateCallback(ICallback *pCallback, IHazardRate **ppHR)
{
  try
  {
    *ppHR = new HazardRateCallbackImpl(pCallback);

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

#if 0

STDMETHODIMP
GlobalIHGImpl::NewFlatHazardRate(double dValue, IHazardRateFlat **ppHR)
{
  try
  {
    *ppHR = new HazardRateFlatImpl(new ihg::HazardRateFlat(dValue));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

#endif // 0

STDMETHODIMP
GlobalIHGImpl::NewHazardRateTimeOnly(SAFEARRAY(DATE) *dates,
                                     SAFEARRAY(double) *values,
                                     IHazardRateTimeOnly **ppHR)
{
  try
  {
    *ppHR = new HazardRateTimeOnlyImpl(shared_ptr<ihg::HazardRateTimeOnly>(
              new ihg::HazardRateTimeOnly(C2A::COM::ToVector<Date, DATE>(dates),
                                          C2A::COM::ToVector<double, double>(values))
            ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
GlobalIHGImpl::NewHazardRatePower(double dAlpha,
                                  double dBeta,
                                  double dS0,
                                  IHazardRatePower **ppHR)
{
  try
  {
    *ppHR = new HazardRatePowerImpl(shared_ptr<ihg::HazardRatePower>(
              new ihg::HazardRatePower(dAlpha, dBeta, dS0)
            ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}


STDMETHODIMP
GlobalIHGImpl::NewHRSpotComponentPower(double dBeta,
                                       double dS0,
                                       IHRSpotComponentPower **ppHRSpotComponent)
{
  try
  {
    *ppHRSpotComponent
         = new HRSpotComponentPowerImpl(shared_ptr<ihg::HRSpotComponentPower>
               ( new ihg::HRSpotComponentPower(dBeta, dS0)));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}


STDMETHODIMP 
GlobalIHGImpl::NewHazardRateCombo(ISpotComponent *pSpotComponent,
                                  SAFEARRAY(DATE) *dates,
                                  SAFEARRAY(double) *values,
                                  IHazardRateCombo **ppHR)
{
  try
  {
    *ppHR = new HazardRateComboImpl(shared_ptr<ihg::HazardRateCombo>
                ( new ihg::HazardRateCombo
                      (
                        C2A::COM::Translate<ISpotComponent *>::From(pSpotComponent),
                        C2A::COM::ToVector<Date, DATE>(dates),
                        C2A::COM::ToVector<double, double>(values)
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

/*
STDMETHODIMP 
GlobalIHGImpl::NewImpliedCDSSpreads(DATE contractingDate,
                                    DATE firstPaymentDate,
                                    DayCountConvention dcc,
                                    Frequency freq,
                                    double dRecoveryRate,
                                    IImpliedCDSSpreads **ppImpliedCDSSpreads)
{
  try
  {
    *ppImpliedCDSSpreads = new ImpliedCDSSpreadsImpl
                               ( new ihg::ImpliedCDSSpreads
                                     (
                                       C2A::COM::Translate<ito33::Date>::From(contractingDate),
                                       C2A::COM::Translate<ito33::Date>::From(firstPaymentDate),
                                       static_cast<Date::DayCountConvention>(dcc),
                                       static_cast<finance::Frequency>(freq),
                                       dRecoveryRate
                                     )
                               );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
*/

STDMETHODIMP 
GlobalIHGImpl::NewParametrizationHRWithSpotComponentPower(IVolatility *pVol, 
                                                          IParametrizationHRWithSpotComponentPower **ppParametrization)
{
  try
  {
    *ppParametrization = new ParametrizationHRWithSpotComponentPowerImpl
                             ( shared_ptr<ihg::ParametrizationHRWithSpotComponentPower>(
                                new ihg::ParametrizationHRWithSpotComponentPower
                                    (
                                      C2A::COM::Translate<IVolatility *>::From(pVol)
                                    )
                              ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalIHGImpl::NewParametrizationVolPower(IHazardRate *pHR, 
                                          IParametrizationVolPower **ppParametrization)
{
  try
  {
    *ppParametrization = new ParametrizationVolPowerImpl
                             (shared_ptr<ihg::ParametrizationVolPower>(
                                new ihg::ParametrizationVolPower
                                    (
                                      C2A::COM::Translate<IHazardRate *>::From(pHR)
                                    )
                             ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalIHGImpl::NewParametrizationVolTanh(IHazardRate *pHR, 
                                         IParametrizationVolTanh **ppParametrization)
{
  try
  {
    *ppParametrization = new ParametrizationVolTanhImpl
                             (shared_ptr<ihg::ParametrizationVolTanh>(
                                new ihg::ParametrizationVolTanh
                                    (
                                      C2A::COM::Translate<IHazardRate *>::From(pHR)
                                    )
                             ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalIHGImpl::NewParametrizationVolWithTimeComponent(IHazardRate *pHR, 
                                                      IParametrizationVolWithTimeComponent **ppParametrization)
{
  try
  {
    *ppParametrization = new ParametrizationVolWithTimeComponentImpl
                             (
                                shared_ptr<ihg::ParametrizationVolWithTimeComponent>(
                                new ihg::ParametrizationVolWithTimeComponent
                                    (
                                      C2A::COM::Translate<IHazardRate *>::From(pHR)
                                    )
                             ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
