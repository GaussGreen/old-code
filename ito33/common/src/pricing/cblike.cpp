/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cblike.cpp
// Purpose:     Implementation of the contract for CB-like
// Created:     2004/08/19
// RCS-ID:      $Id: cblike.cpp,v 1.25 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/issuer.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/numeraire.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/convertiblelike.h"

#include "ito33/pricing/cashflows.h"
#include "ito33/pricing/cblike.h"

namespace ito33
{

namespace pricing
{


void 
CBLike::GetBondLikeTermsData(const finance::BondLikeTerms& blt, 
                             const finance::SessionData& sessionData,
                             const shared_ptr<finance::Numeraire>& pNumeraire)
{
  m_dRecoveryRate = blt.GetRecoveryRate();

  m_dMaturityTime = GetDoubleFrom( blt.GetMaturityDate() );
  
  m_dNominal = blt.GetNominal();
  
  m_dIssuePrice = blt.GetIssuePrice() * m_dNominal;

  m_dIssueTime = GetDoubleFrom( blt.GetIssueDate() );

  m_fiscalYearStartDate = 
    sessionData.GetEquity()->GetIssuer()->GetFiscalYearStartDate();  

  GetCrossCurrencyData(sessionData, pNumeraire);
}

void 
CBLike::GetCrossCurrencyData(const finance::SessionData& sessionData,
                             const shared_ptr<finance::Numeraire>& pNumeraire)
{
  // Check for cross-currency.  
  m_bIsCrossCurrency  
    = ( pNumeraire &&
      ( *pNumeraire != *(sessionData.GetEquity()->GetNumeraire()) ) );

  if ( m_bIsCrossCurrency )
  {
    m_pDerivativeCurve = 
      sessionData.GetRateData()->GetYieldCurve( pNumeraire );
    
    m_dSpotFXRate = sessionData.GetSpotFXRate( 
      sessionData.GetEquity()->GetNumeraire(), pNumeraire );
  }
  else
    m_pDerivativeCurve = sessionData.GetYieldCurve();

}

void CBLike::GetConvertibleLikeData(const finance::ConvertibleLike& cbLike)
{  
  m_bNewShare = cbLike.GetConvertIntoNewShare();

  m_bIsFixedQuanto = cbLike.IsFixedQuanto();

  GetConversions()->SetTriggerAsPercentageOf
                    ( cbLike.GetConversionTriggerAsPercentageOf() );

  if (m_bIsCrossCurrency)
  {
    GetConversions()->SetTriggerInCurrencyOf(cbLike.GetTriggerInCurrencyOf());    
    
    GetConversions()->SetFixedFXRate(cbLike.GetFixedFXRate());
    
    if ( m_bIsFixedQuanto )
    {
      m_dFXRateVolatility = cbLike.GetFXRateVolatility();

      m_dCorrelation = cbLike.GetCorrelationBetweenUnderlyingAndFXRate();
    }
  }

  if (cbLike.IsExchangeable() )
    SetExchangeable(cbLike.IsExchangeableUponDefault());
}

double CBLike::GetFXRateVolatility() const
{
  ASSERT_MSG( m_bIsFixedQuanto, "FX rate volatility required for a non Fixed "
    "quanto");

  return m_dFXRateVolatility;
}

double CBLike::GetCorrelationBetweenUnderlyingAndFXRate() const
{
  ASSERT_MSG( m_bIsFixedQuanto, "Correlation required for a non Fixed quanto");

  return m_dCorrelation;
}

} //namespace pricing

} //namespace ito33
