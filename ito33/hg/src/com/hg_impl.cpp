/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/com/hg_impl.cpp
// Purpose:     declaration and implementation of IHG interface
// Created:     2005/01/17
// RCS-ID:      $Id: hg_impl.cpp,v 1.9 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/hg/showversion.h"
#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"
#include "ito33/com/unknown_impl.h"
#include "ito33/com/c2a.h"

#include "com/varianceswap_impl.h"

#include "com/theoreticalmodel_impl.h"
#include "com/parametrization_impl.h"

#include "hg.h"

using namespace ito33;
using std::string;

DEFINE_COM_IFACE(GlobalHG, IUnknown);
 
// ----------------------------------------------------------------------------
// HGImpl declaration
// ----------------------------------------------------------------------------

class GlobalHGImpl : public COM::ImplementCoClass
                                 <
                                   IGlobalHG,
                                   GlobalHGImpl,
                                   TYPE_LIST_2(ISupportErrorInfo, IDispatch)
                                 >
{
public:

  GlobalHGImpl() { }

  STDMETHODIMP get_Version(BSTR *version);

  STDMETHODIMP NewUnderlyingProcess(long nNbRegimes,
                                    SAFEARRAY(double) *vols,
                                    SAFEARRAY(double) *defaultIntensities,
                                    IUnderlyingProcess **ppUnderlyingProcess);
     
  STDMETHODIMP NewTheoreticalModel(IUnderlyingProcess *pUnderlyingProcess, 
                                   ITheoreticalModel **ppTM);

  STDMETHODIMP NewParametrization(IUnderlyingProcess *pUnderlyingProcess, 
                                  IParametrization **ppPM);

  STDMETHODIMP ComputeVarianceSwapByLog(IVarianceSwap *pVarianceSwap,
                                        ITheoreticalModel *pTM,
                                        IModelOutput **ppMO);

  STDMETHODIMP ComputeImpliedVolatilityStrikeByLog(IVarianceSwap *pVarianceSwap,
                                                   ITheoreticalModel *pTM,
                                                   double *pVol);

private:

  GlobalHGImpl(const GlobalHGImpl&);
  GlobalHGImpl& operator=(const GlobalHGImpl&);
};


// ============================================================================
// HG implementation
// ============================================================================

DEFINE_COM_COCLASS_FULL(CLSID_GlobalHG, "HG.GlobalHG", 1, GlobalHGImpl);

STDMETHODIMP 
GlobalHGImpl::get_Version(BSTR *version)
{
  try
  {
    *version = C2A::COM::Translate<string>::To( ito33::hg::ShowVersion() );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalHGImpl::NewUnderlyingProcess(long nNbRegimes,
                                   SAFEARRAY(double) *vols,
                                   SAFEARRAY(double) *defaultIntensities,
                                   IUnderlyingProcess **ppUnderlyingProcess)
{
  try
  {
    *ppUnderlyingProcess = new UnderlyingProcessImpl
                           (
                             shared_ptr<hg::UnderlyingProcess>
                             ( new hg::UnderlyingProcess(nNbRegimes, 
                                                         C2A::COM::ToVector<double, double>(vols),
                                                         C2A::COM::ToVector<double, double>(defaultIntensities)
                                                        )
                             )
                           );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalHGImpl::NewTheoreticalModel(IUnderlyingProcess *pUnderlyingProcess, 
                                  ITheoreticalModel **ppTM)
{
  try
  {
    *ppTM = new TheoreticalModelImpl(
              shared_ptr<hg::TheoreticalModel>
              ( 
                new hg::TheoreticalModel
                (C2A::COM::Translate<IUnderlyingProcess *>::From(pUnderlyingProcess)) 
              ) 
            );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
GlobalHGImpl::NewParametrization(IUnderlyingProcess *pUnderlyingProcess, 
                                 IParametrization **ppPM)
{
  try
  {
    *ppPM = new ParametrizationImpl(
              shared_ptr<hg::Parametrization>
              ( new hg::Parametrization
                    (C2A::COM::Translate<IUnderlyingProcess *>::From(pUnderlyingProcess)) 
              )
            );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
GlobalHGImpl::ComputeVarianceSwapByLog(IVarianceSwap *pVarianceSwap,
                                       ITheoreticalModel *pTM,
                                       IModelOutput **ppMO)
{
  try
  {
    shared_ptr<hg::TheoreticalModel>
      pTMImpl( C2A::COM::Translate<ITheoreticalModel *>::From(pTM) );

    shared_ptr<finance::VarianceSwap>
      pVSImpl( C2A::COM::Translate<IVarianceSwap *>::From(pVarianceSwap) );

    *ppMO = new ModelOutputImpl(
                  pTMImpl->Compute
                  ( *pVSImpl, 
                    reinterpret_cast<hg::TheoreticalModel::ModelPriceFunction>
                    (hg::TheoreticalModel::PriceVarianceSwapByLog)
                ));

     return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

// The following code is mainly copied from ImpliedParameterCalculator
// since it sounds not worth it to try to share the code
// This code is really temporary
STDMETHODIMP
GlobalHGImpl::ComputeImpliedVolatilityStrikeByLog(IVarianceSwap *pVarianceSwap,
                                                  ITheoreticalModel *pTM,
                                                  double *pVol)
{
  try
  {
    shared_ptr<hg::TheoreticalModel>
      pTMImpl( C2A::COM::Translate<ITheoreticalModel *>::From(pTM) );

    shared_ptr<finance::VarianceSwap>
      pVSImpl( C2A::COM::Translate<IVarianceSwap *>::From(pVarianceSwap) );

    double dPrice = 0.;
    if ( pVSImpl->HasMarketPrice() )
      dPrice = pVSImpl->GetMarketPrice();  

    // Clone the model, don't use the flags inside the original model
    shared_ptr<hg::TheoreticalModel> pMyModel( pTMImpl->Clone() );

    // Don't use flags inside Derivative either, use default flags
    pMyModel->SetExternalFlagsToDefaults();

    const double dVolStrikeTmp = 0.1;
    finance::VarianceSwap vsTmp(*pVSImpl);
    vsTmp.SetVolatilityStrike(dVolStrikeTmp);

    double dPriceTmp = pMyModel->Compute
                       ( vsTmp, 
                         reinterpret_cast<hg::TheoreticalModel::ModelPriceFunction>
                         (hg::TheoreticalModel::PriceVarianceSwapByLog)
                       )->GetPrice();
    
    const finance::SessionData& sessionData( *pVSImpl->GetSessionData() );

    double dSlope = - sessionData.GetYieldCurve()->GetForwardDiscountFactor
                      ( GetDoubleFrom( sessionData.GetValuationDate() ), 
                        GetDoubleFrom( pVSImpl->GetMaturityDate() ) );

    double dVolStrike = (dPrice - dPriceTmp) / dSlope;
    dVolStrike += dVolStrikeTmp * dVolStrikeTmp;
 
    if ( dVolStrike > 0 )
      dVolStrike = sqrt(dVolStrike);
    else // don't allow negative value because of round off error
      dVolStrike = 0.;

    *pVol = dVolStrike;

     return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
