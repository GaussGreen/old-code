/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/impliedparametercalculator.cpp
// Purpose:     global functions for computing implied spread/strike etc
// Created:     2006/06/01
// RCS-ID:      $Id: impliedparametercalculator.cpp,v 1.12 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/useexception.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/derivative_visitor.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/gammavarianceswap.h"

#include "ito33/finance/impliedparametercalculator.h"

extern const ito33::finance::Error
  ITO33_VARIANCESWAP_IMPLIEDSTRIKE_NOT_VANILLA,
  ITO33_VARIANCESWAP_IMPLIEDSTRIKE_NOT_IMPLEMENTED,
  ITO33_VARIANCESWAP_IMPLIEDSTRIKE_FOR_STARTED;

namespace ito33
{

namespace finance 
{

/* static */ double   
ImpliedParameterCalculator::ComputeImpliedCDSSpread
(const CDSLike& cds, const TheoreticalModel& model)
{
  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( cds.GetSpreadStream() );

  // Clone the model, don't use the flags inside the original model
  shared_ptr<TheoreticalModel> pMyModel( model.Clone() );

  // Don't use flags inside Derivative either, use default flags
  pMyModel->SetExternalFlagsToDefaults(); 

  const double dSpread1 = 0.1;
  const double dSpread2 = 0.2;

  CDS cds1(cds.GetRecoveryRate(), pSpreadStream->ChangeSpread(dSpread1));
  cds1.SetSessionData( cds.GetSessionData() );

  const double dPrice1 = pMyModel->Compute(cds1)->GetPrice();

  CDS cds2(cds.GetRecoveryRate(), pSpreadStream->ChangeSpread(dSpread2));
  cds2.SetSessionData( cds.GetSessionData() );

  const double dPrice2 = pMyModel->Compute(cds2)->GetPrice();

  const double dSlope = (dPrice2 - dPrice1) / (dSpread2 - dSpread1);
  const double dC = dPrice1 - dSlope * dSpread1;
  
  if ( dSlope >= 0 )
    return 0;

  return - dC / dSlope;
}

/* static */ double   
ImpliedParameterCalculator::ComputeImpliedEDSSpread
(const EDS& eds, const TheoreticalModel& model)
{
  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( eds.GetSpreadStream() );

  // Clone the model, don't use the flags inside the original model
  shared_ptr<TheoreticalModel> pMyModel( model.Clone() );

  // Don't use flags inside Derivative either, use default flags
  pMyModel->SetExternalFlagsToDefaults(); 

  const double dSpread1 = 0.1;
  const double dSpread2 = 0.2;

  EDS eds1( eds.GetRecoveryRate(), pSpreadStream->ChangeSpread(dSpread1),  
            eds.GetBarrier() );
  eds1.SetSessionData( eds.GetSessionData() );

  const double dPrice1 = pMyModel->Compute(eds1)->GetPrice();

  EDS eds2( eds.GetRecoveryRate(), pSpreadStream->ChangeSpread(dSpread2),
            eds.GetBarrier() );
  eds2.SetSessionData( eds.GetSessionData() );

  const double dPrice2 = pMyModel->Compute(eds2)->GetPrice();

  const double dSlope = (dPrice2 - dPrice1) / (dSpread2 - dSpread1);
  const double dC = dPrice1 - dSlope * dSpread1;
  
  if ( dSlope >= 0 )
    return 0;

  return - dC / dSlope;
}

/// Visitor class helps to see if implied vol strike can be computed
class VSVisitor : protected DerivativeVisitor
{
public:
  
  VSVisitor(const VarianceSwapLike& vs) : m_dVolStrike(0.1)
  {
    try
    {
      vs.Visit(*this);
    }
    catch( const DerivativeVisitor::Exception& )
    {
      throw EXCEPTION(ITO33_VARIANCESWAP_IMPLIEDSTRIKE_NOT_IMPLEMENTED);
    }
  }

  shared_ptr<VarianceSwapLike> operator()() const { return m_pVS; }

  virtual void OnVarianceSwap(const VarianceSwap& vs)
  {
    shared_ptr<VarianceSwap> pVarianceSwap( new VarianceSwap(vs) );
    pVarianceSwap->SetVolatilityStrike(m_dVolStrike);

    m_pVS = pVarianceSwap;
  }
 
  virtual void OnGammaVarianceSwap(const GammaVarianceSwap& vs)
  {
    shared_ptr<GammaVarianceSwap> pVarianceSwap( new GammaVarianceSwap(vs) );
    pVarianceSwap->SetVolatilityStrike(m_dVolStrike);

    m_pVS = pVarianceSwap;
  }

private:
  const double m_dVolStrike;
 
  shared_ptr<VarianceSwapLike> m_pVS;

  NO_COPY_CLASS(VSVisitor);
};

/* static */ double 
ImpliedParameterCalculator::ComputeImpliedVolatilityStrike
(const VarianceSwapLike& varianceSwap, const TheoreticalModel& model)
{  
  double dPrice = 0.;
  if ( varianceSwap.HasMarketPrice() )
    dPrice = varianceSwap.GetMarketPrice();  
  
  CHECK_COND( varianceSwap.GetTerms()->GetCapMultiplier() <= 0,
              ITO33_VARIANCESWAP_IMPLIEDSTRIKE_NOT_VANILLA );

  // Clone the model, don't use the flags inside the original model
  shared_ptr<TheoreticalModel> pMyModel( model.Clone() );

  // Don't use flags inside Derivative either, use default flags
  pMyModel->SetExternalFlagsToDefaults();

  VSVisitor visitor(varianceSwap);
  
  shared_ptr<VarianceSwapLike> pVSTmp = visitor();
  double dVolStrikeTmp = pVSTmp->GetVolatilityStrike();

  double dPriceTmp = pMyModel->Compute(*pVSTmp)->GetPrice();
  
  const SessionData& sessionData( *varianceSwap.GetSessionData() );

  double dSlope = - sessionData.GetYieldCurve()->GetForwardDiscountFactor
                    ( GetDoubleFrom( sessionData.GetValuationDate() ), 
                      GetDoubleFrom( varianceSwap.GetMaturityDate() ) );

  double dVolStrike = (dPrice - dPriceTmp) / dSlope;
  if ( varianceSwap.GetTerms()->GetSwapType() == Swap_Variance )
    dVolStrike += dVolStrikeTmp * dVolStrikeTmp;
  else
    dVolStrike += dVolStrikeTmp;
  
  if ( varianceSwap.GetTerms()->GetSwapType() == Swap_Variance )
    if ( dVolStrike > 0 )
      dVolStrike = sqrt(dVolStrike);
    else // don't allow negative value because of round off error
      dVolStrike = 0.;

  return dVolStrike;
}

/* static */ double 
ImpliedParameterCalculator::ComputeImpliedVolatilityStrike
(const shared_ptr<VarianceSwapTerms>& pTerms,
 const shared_ptr<SessionData>& pSessionData,
 const TheoreticalModel& model)
{  
  CHECK_COND(     pTerms->GetStartOfSamplingPeriod()
               >= pSessionData->GetValuationDate(),
               ITO33_VARIANCESWAP_IMPLIEDSTRIKE_FOR_STARTED );

  VarianceSwap vs(pTerms, 0.1);
  vs.SetSessionData(pSessionData);

  return ComputeImpliedVolatilityStrike(vs, model);
}

} // namespace finance

} // namespace ito33
