/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/basket_visitor.cpp
// Purpose:     Visitor for Derivative-derived classes used in calibration
// Created:     2006/08/07
// RCS-ID:      $Id: basket_visitor.cpp,v 1.9 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/vector.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"
#include "ito33/finance/derivatives.h"

#include "ito33/finance/basket_visitor.h"

namespace ito33
{

namespace finance
{

typedef std::vector< shared_ptr<Derivative> > Elements;

/**
    Class prepares for the (eventual) updates of the weights according to the
    specified options.
 */
class MarketValueVisitor : protected DerivativeVisitor
{
public:
  MarketValueVisitor(const Elements& elements);

  /**
      Enables/disables advanced target eventually different than market price.

      @param bUseAdvancedTarget Use or not adavanced target
   */
  void EnableAdvancedTarget(bool bUseAdvancedTarget = true)
  {
    m_bUseAdvancedTarget = bUseAdvancedTarget;
  }

  /**
      Enables/disables relative error on target.

      @param bUseRelativeError Use or not relative error on target
   */
  void EnableRelativeError(bool bUseRelativeError = true)
  {
    m_bUseRelativeError = bUseRelativeError;
  }

  /**
      Runs the visitor, updates if necessary the derivative weight.

      @return The derivatives with default weights
   */
  shared_ptr<Derivatives> Run();

  virtual void OnOption(const Option&);

  virtual void OnOneTouch(const OneTouch&);

  virtual void OnFXOneTouch(const FXOneTouch&);

  virtual void OnCDS(const CDS&);

  virtual void OnReferenceCDS(const ReferenceCDS&);

  virtual void OnParBond(const ParBond&);

  virtual void OnEDS(const EDS&);

  virtual void OnVarianceSwap(const VarianceSwap&);

  virtual void OnVarianceSwaption(const VarianceSwaption&);

private:
  
  /// For generic derivatives, the target is always the market price.
  void OnGeneric();

  /// For cdslike
  void OnCDSLike();

  /// The result
  shared_ptr<Derivatives> m_pDerivatives;

  /// The derivatives
  const Elements& m_elements;

  /// If advanced target that might be different than market price will be used
  bool m_bUseAdvancedTarget;

  /// If relative error of the target will be used
  bool m_bUseRelativeError;

  /// Current derivative
  shared_ptr<Derivative> m_pDerivative;

  /// Current market value
  double m_dMarketValue;

  /// Current coefficient to transform between market price and value
  double m_dCoe;

  NO_COPY_CLASS(MarketValueVisitor);
};

MarketValueVisitor::MarketValueVisitor(const Elements& elements)
                                     : m_elements(elements),
                                       m_bUseAdvancedTarget(false),
                                       m_bUseRelativeError(false)
{
  m_pDerivatives = make_ptr( new Derivatives ); 
}

shared_ptr<Derivatives> MarketValueVisitor::Run()
{
  Elements::const_iterator iter;
  for (iter = m_elements.begin(); iter != m_elements.end(); ++iter)
  {
    m_pDerivative = *iter;

    // Validate at first the derivative since we'll make use of its details
    m_pDerivative->Validate();

    try
    {
      m_pDerivative->Visit(*this);
    }
    catch( const DerivativeVisitor::Exception& )
    {
      OnGeneric();
    }  
    
    double dScale = m_dCoe;
    if ( m_bUseRelativeError )
      dScale *= m_dMarketValue;

    dScale = fabs(dScale);

    if ( dScale < 1.e-5 ) // value taken from old code, need review
      dScale = 1;
   
    m_pDerivatives->AddWithWeight(m_pDerivative, 1. / (dScale * dScale) );
  }

  return m_pDerivatives;
}

void MarketValueVisitor::OnGeneric()
{
  m_dMarketValue = m_pDerivative->GetMarketPrice();
  
  m_dCoe = 1;
}

void MarketValueVisitor::OnOption(const Option&)
{
  // Current derivative is an option.
  shared_ptr<Option> pOption( static_pointer_cast<Option>(m_pDerivative) );

  // Cache the option market price by using a copy if the market price
  // is given by impliedvol
  if ( pOption->IsImpliedVolSet() )
  {
    shared_ptr<Option> pNewOption( new Option(*pOption) );
    pNewOption->SetMarketPrice( pOption->GetMarketPrice() );
  
    m_pDerivative = pNewOption;
  }

  if ( !m_bUseAdvancedTarget )
  {
    OnGeneric();
    return;
  }

  // Otherwise, vega can be used to transform price to implied vol
  m_dMarketValue = pOption->GetImpliedVol();

  m_dCoe = pOption->GetVegaFrom(m_dMarketValue);
}

void MarketValueVisitor::OnFXOneTouch(const FXOneTouch&)
{
  // Current derivative is a FX one touch
  shared_ptr<FXOneTouch>
    pFXOneTouch( static_pointer_cast<FXOneTouch>(m_pDerivative) );

  // Cache the barrier level by creating a new one touch
  shared_ptr<OneTouch> 
    pOneTouch( new OneTouch( pFXOneTouch->GetMaturityDate(),
                             pFXOneTouch->GetBarrier(),
                             pFXOneTouch->GetBarrierType(),
                             pFXOneTouch->GetRebateType() ) );

  pOneTouch->SetSessionData( pFXOneTouch->GetSessionData() );
  pOneTouch->SetMarketPrice( pFXOneTouch->GetMarketPrice() );

  m_pDerivative = pOneTouch;

  OnGeneric();
}

void MarketValueVisitor::OnOneTouch(const OneTouch&)
{
  OnGeneric();
}

void MarketValueVisitor::OnCDSLike()
{
  if ( !m_bUseAdvancedTarget )
  {
    OnGeneric();
    return;
  }

  // Current derivative is a cdslike
  shared_ptr<CDSLike> pCDS( static_pointer_cast<CDSLike>(m_pDerivative) );

  m_dMarketValue = pCDS->GetSpread();
  
  const SessionData& sessionData( *pCDS->GetSessionData() );

  m_dCoe = pCDS->GetSpreadStream()->GetDiscount
           ( *sessionData.GetYieldCurve(), sessionData.GetValuationDate() );

  m_dCoe /= m_dMarketValue;
}

void MarketValueVisitor::OnCDS(const CDS&)
{
  OnCDSLike();
}

void MarketValueVisitor::OnReferenceCDS(const ReferenceCDS&)
{
  OnCDSLike();
}

void MarketValueVisitor::OnParBond(const ParBond&)
{
  OnGeneric();
}

void MarketValueVisitor::OnEDS(const EDS&)
{
  if ( !m_bUseAdvancedTarget )
  {
    OnGeneric();
    return;
  }

  // Current derivative is an EDS
  shared_ptr<EDS> pEDS( static_pointer_cast<EDS>(m_pDerivative) );

  m_dMarketValue = pEDS->GetSpreadStream()->GetAnnualPaymentAmount();
  
  const SessionData& sessionData( *pEDS->GetSessionData() );

  m_dCoe = pEDS->GetSpreadStream()->GetDiscount
           ( *sessionData.GetYieldCurve(), sessionData.GetValuationDate() );

  m_dCoe /= m_dMarketValue;
}

void MarketValueVisitor::OnVarianceSwap(const VarianceSwap&)
{
  if ( !m_bUseAdvancedTarget )
  {
    OnGeneric();
    return;
  }

  // Current derivative is a variance swap
  shared_ptr<VarianceSwap> 
    pVS( static_pointer_cast<VarianceSwap>(m_pDerivative) );
  const SessionData& sessionData( *pVS->GetSessionData()  );

  double dSlope = - sessionData.GetYieldCurve()->GetForwardDiscountFactor
                    ( sessionData.GetValuationDate(), pVS->GetMaturityDate() );

  if ( pVS->GetTerms()->GetSwapType() == Swap_Variance )
    m_dCoe = 2 * dSlope * pVS->GetVolatilityStrike();
  else
    m_dCoe = dSlope;

  m_dMarketValue = pVS->GetVolatilityStrike();
}

void MarketValueVisitor::OnVarianceSwaption(const VarianceSwaption&)
{
  OnGeneric();
}

BasketVisitor::BasketVisitor(const Derivatives& derivatives)
                           : m_derivatives(derivatives),
                             m_bUseForward(true),
                             m_bUseAdvancedTarget(false),
                             m_bUseRelativeError(false)   
{
}

void BasketVisitor::Run()
{
  // First pass, update derivative weights
  Elements elements;
  Derivatives::Elements::const_iterator iterOrig;
  for (iterOrig = m_derivatives.begin(); 
       iterOrig != m_derivatives.end();
       ++iterOrig)
    elements.push_back(iterOrig->first);
  
  MarketValueVisitor visitor(elements);

  visitor.EnableAdvancedTarget(m_bUseAdvancedTarget);
  visitor.EnableRelativeError(m_bUseRelativeError);

  shared_ptr<Derivatives> pDerivatives( visitor.Run() );

  // Second pass, group same type of instrument to speed up
  Derivatives::const_iterator iter;
  for (iter = pDerivatives->begin(), iterOrig = m_derivatives.begin(); 
       iter != pDerivatives->end(); 
       ++iter, ++iterOrig)
  {
    m_pDerivative = iter->first;
    m_dWeight = iterOrig->second * iter->second;

    try
    {
      m_pDerivative->Visit(*this);
    }
    catch( const DerivativeVisitor::Exception& )
    {
      OnGeneric();
    }
  }
}

void BasketVisitor::OnGeneric()
{
  // Current derivative is a generic one
  if ( !m_pGenericDerivatives )
    m_pGenericDerivatives = make_ptr(new Derivatives);

  m_pGenericDerivatives->AddWithWeight(m_pDerivative, m_dWeight);
}

void BasketVisitor::OnOption(const Option&)
{
  // Current derivative is an option.
  shared_ptr<Option> pOption( static_pointer_cast<Option>(m_pDerivative) );

  if ( pOption->GetExerciseType() == ExerciseType_European && m_bUseForward )
  {
    if ( !m_pForwardOption )
      m_pForwardOption = make_ptr(new ForwardOption);

    m_pForwardOption->Add(pOption, m_dWeight);
  }
  else
    OnGeneric();
}

void BasketVisitor::OnOneTouch(const OneTouch&)
{
  // TODO: group one touch if applicable
  OnGeneric();
}

} // namespace finance

} // namespace ito33
