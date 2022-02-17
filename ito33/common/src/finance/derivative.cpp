/////////////////////////////////////////////////////////////////////////////
// Name:        common/src//finance/derivative.cpp
// Purpose:     do the necessary for Derivative class
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: derivative.cpp,v 1.25 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/constants.h"

#include "ito33/finance/error.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/numeraire.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"

extern const ito33::finance::Error
  ITO33_INVALID_MARKETPRICE_NEGATIVE,
  ITO33_INVALID_MARKETPRICE_TOOBIGTOBESET_1,
  ITO33_NO_MARKETPRICE,
  ITO33_MATURITYBEFOREVALUATION,
  ITO33_DERIVATIVE_INVALID_SESSIONDATA,
  ITO33_DERIVATIVE_UNDEFINED_SESSIONDATA,
  ITO33_DERIVATIVE_VISITOR_NOT_IMPLEMENTED,
  ITO33_DERIVATIVE_INVALID_NUMERAIRE,
  ITO33_NULL_COMPUTATIONALFLAGS;

namespace ito33
{

namespace finance
{


Derivative::Derivative()           
{
  UnsetMarketPrice(); 
}

void Derivative::UnsetMarketPrice()
{
  m_dPrice = INVALIDPRICE;
}

bool Derivative::MarketPriceIsSet() const
{
  return m_dPrice < INVALIDPRICE;
}

void Derivative::SetSessionData(const shared_ptr<SessionData>& pSessionData)
{
  CHECK_COND( pSessionData, ITO33_DERIVATIVE_INVALID_SESSIONDATA );

  m_pSessionData = pSessionData;
}

void Derivative::SetComputationalFlags
     ( const shared_ptr<ComputationalFlags>& pFlags )
{  
  CHECK_COND(pFlags, ITO33_NULL_COMPUTATIONALFLAGS);
  
  m_pFlags = pFlags;
}

void Derivative::SetNumeraire(const shared_ptr<Numeraire>& pCurrency)
{
  CHECK_COND( pCurrency, ITO33_DERIVATIVE_INVALID_NUMERAIRE );

  m_pNumeraire = pCurrency;
}

void Derivative::CheckSessionData() const
{
  // Check that the session data is set
  CHECK_COND( m_pSessionData, ITO33_DERIVATIVE_UNDEFINED_SESSIONDATA );
}

void Derivative::ValidateWith(const SessionData& sessionData) const
{
  // Validate at first the derivative itself
  Validate();
  
  // cross validation with the session data
  CHECK_COND( GetMaturityDate() > sessionData.GetValuationDate(),
              ITO33_MATURITYBEFOREVALUATION );
}

void Derivative::ValidateAll() const
{
  // Check that the session data is set
  CheckSessionData();

  // Validate with the session data
  ValidateWith(*m_pSessionData);
}

const shared_ptr<SessionData>& Derivative::GetSessionData() const
{
  return m_pSessionData;
}

void Derivative::SetMarketPrice(double dPrice) 
{
  CHECK_COND( dPrice > 0, ITO33_INVALID_MARKETPRICE_NEGATIVE );

  DoSetMarketPrice(dPrice);
}

void Derivative::DoSetMarketPrice(double dPrice)
{
  m_dPrice = dPrice;

  CHECK_COND_1 ( MarketPriceIsSet(),
                 ITO33_INVALID_MARKETPRICE_TOOBIGTOBESET_1,
                 m_dPrice);
  // Actually MarketPriceIsSet() is (m_dPrice < INVALIDPRICE).
  // This check is not that readable but the rational here
  // is that if DoSetMarketPrice() sucesses, MarketPriceIsSet() must
  // be true 
}

bool Derivative::HasMarketPrice() const
{
  return MarketPriceIsSet();
}

void Derivative::CheckMarketPrice() const
{
  CHECK_COND( HasMarketPrice(), ITO33_NO_MARKETPRICE );
}

bool Derivative::IsCrossCurrency(const SessionData& sessionData) const
{
  return m_pNumeraire && *m_pNumeraire != *( sessionData.GetNumeraire() );
}

bool Derivative::IsCrossCurrency() const
{
  CheckSessionData();

  return IsCrossCurrency(*m_pSessionData);
}

void Derivative::PerturbFXRate(double dFXRateShift) const
{
  ASSERT_MSG( IsCrossCurrency(), 
              "Cannot perturb FX rate for a non cross-currency derivative!" );

  double dFXRate;

  dFXRate = GetSessionData()->GetRateData()->GetSpotFXRates()->GetFXRate
            ( GetSessionData()->GetNumeraire(), m_pNumeraire );

  dFXRate += dFXRateShift;

  GetSessionData()->GetRateData()->GetSpotFXRates()->SetFXRate
  (
    GetSessionData()->GetEquity()->GetNumeraire(), 
    m_pNumeraire,
    dFXRate
  );
}

void Derivative::DumpMe(XML::Tag& tagParent) const
{
  if ( m_pIssuer )
    m_pIssuer->Dump(tagParent);

  if ( m_pNumeraire )
    tagParent.Element(XML_TAG_CURRENCY)( m_pNumeraire->GetCode() );
}

void Derivative::DumpMarketPrice(XML::Tag& tagParent) const
{
  if ( MarketPriceIsSet() )
    tagParent.Element(XML_TAG_FINANCE_MARKETPRICE)(m_dPrice);
}

void Derivative::Visit(DerivativeModifyingVisitor&)
{
  throw EXCEPTION(ITO33_DERIVATIVE_VISITOR_NOT_IMPLEMENTED);
} 


} // namespace finance

} // namespace ito33
