/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/convertiblelike.cpp
// Purpose:     financial convertible-like bond class
// Author:      ITO33
// Created:     2004/10/13
// RCS-ID:      $Id: convertiblelike.cpp,v 1.26 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/issuer.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/cashflowstream.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"
#include "ito33/finance/bondlike/triggerincurrencyof.h"
#include "ito33/finance/termstructurecds.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"
#include "ito33/xml/finance/bondlike/trigger.h"

extern const ito33::finance::BondError
  ITO33_BONDLIKE_CROSSCURRENCY_FIXEDFXRATE,
  ITO33_BONDLIKE_NOT_CROSSCURRENCY,
  ITO33_BONDLIKE_NOT_EXCHANGEABLE,
  ITO33_BONDLIKE_TRIGGERASPERCENTAGEOF,
  ITO33_BONDLIKE_TRIGGERINCURRENCYOF,
  ITO33_BONDLIKE_NULL_BONDLIKETERMS,
  ITO33_BONDLIKE_NEWSHARE_AND_EXCHANGEABLE,
  ITO33_BONDLIKE_FIXEDQUANTO_FXVOLATILITY,
  ITO33_BONDLIKE_FIXEDQUANTO_CORRELATION,
  ITO33_BONDLIKE_NOT_FIXEDQUANTO,
  ITO33_BONDLIKE_FIXED_QUANTO_NOT_CROSSCURRENCY,
  ITO33_BONDLIKE_FIXEDFXRATE_NOT_SET;

namespace ito33
{

namespace finance
{

ConvertibleLike::ConvertibleLike
                 (const shared_ptr<BondLikeTerms>& pBondLikeTerms)
  : Derivative(), 
    m_pBondLikeTerms( CHECK_PTR(pBondLikeTerms,
                                ITO33_BONDLIKE_NULL_BONDLIKETERMS) ),
    m_bIsExchangeable(false),
    m_triggerAsPercentageOf(TriggerAsPercentageOf_Principal),
    m_triggerInCurrencyOf(TriggerInCurrencyOf_Derivative),
    m_bConvertIntoNewShare(false),
    m_bIsFixedQuanto(false),
    m_dFixedFXRate(0.0),
    m_dFXRateVolatility(0),
    m_dCorrelation(0) 
{ 
}

void ConvertibleLike::SetExchangeable(bool bExchangeableUponDefault)
{
  CHECK_COND( !m_bConvertIntoNewShare, 
              ITO33_BONDLIKE_NEWSHARE_AND_EXCHANGEABLE );

  m_bIsExchangeable = true;
  m_bExchangeableUponDefault = bExchangeableUponDefault;
}

void ConvertibleLike::SetConvertIntoNewShare(bool bConvertIntoNewShare)
{ 
  CHECK_COND( !( bConvertIntoNewShare && IsExchangeable() ), 
              ITO33_BONDLIKE_NEWSHARE_AND_EXCHANGEABLE );

  m_bConvertIntoNewShare = bConvertIntoNewShare;
}

bool ConvertibleLike::IsExchangeable() const
{
  return m_bIsExchangeable;
}

bool ConvertibleLike::IsExchangeableUponDefault() const
{
  CHECK_COND(IsExchangeable(), ITO33_BONDLIKE_NOT_EXCHANGEABLE);

  return m_bExchangeableUponDefault;
}

void ConvertibleLike::PerturbFXRate(double dFXRateShift) const
{
  if ( IsFixedQuanto() )
  {
    ASSERT_MSG( IsCrossCurrency(), 
                "Cannot perturb FX rate for a non cross-currency derivative!" );

    m_dFixedFXRate += dFXRateShift;
  }
  else // Not fixed quanto case
    Derivative::PerturbFXRate( dFXRateShift );
}

void 
ConvertibleLike::SetConversionTriggerAsPercentageOf
                 ( TriggerAsPercentageOf asPercentageOf )
{
  CHECK_COND
  ( 
    asPercentageOf < TriggerAsPercentageOf_Max,
    ITO33_BONDLIKE_TRIGGERASPERCENTAGEOF
  );

  m_triggerAsPercentageOf = asPercentageOf;
}

void ConvertibleLike::SetTriggerInCurrencyOf(TriggerInCurrencyOf inCurrencyOf) 
{
  CHECK_COND
  ( 
    inCurrencyOf < TriggerInCurrencyOf_Max,
    ITO33_BONDLIKE_TRIGGERINCURRENCYOF
  );

  m_triggerInCurrencyOf = inCurrencyOf; 
}

void ConvertibleLike::SetFixedFXRate(double dFixedFXRate) 
{
  CHECK_COND(dFixedFXRate > 0., ITO33_BONDLIKE_CROSSCURRENCY_FIXEDFXRATE); 

  m_dFixedFXRate = dFixedFXRate;
}

void ConvertibleLike::SetFixedQuanto(double dFXVolatility, double dCorrelation)
{
  CHECK_COND(dFXVolatility >= 0. && dFXVolatility <= 5., 
             ITO33_BONDLIKE_FIXEDQUANTO_FXVOLATILITY);
  
  CHECK_COND(dCorrelation >= -1. && dCorrelation <= 1., 
             ITO33_BONDLIKE_FIXEDQUANTO_CORRELATION);
  
  m_bIsFixedQuanto = true;
  
  m_dFXRateVolatility = dFXVolatility;

  m_dCorrelation = dCorrelation;
}

TriggerInCurrencyOf ConvertibleLike::GetTriggerInCurrencyOf() const
{ 
  return m_triggerInCurrencyOf;
}

double ConvertibleLike::GetFixedFXRate() const 
{ 
  return m_dFixedFXRate; 
}

bool ConvertibleLike::IsFixedQuanto() const 
{
  return m_bIsFixedQuanto;
}

double ConvertibleLike::GetFXRateVolatility() const 
{ 
  CHECK_COND(IsFixedQuanto(), ITO33_BONDLIKE_NOT_FIXEDQUANTO);

  return m_dFXRateVolatility; 
}

double ConvertibleLike::GetCorrelationBetweenUnderlyingAndFXRate() const 
{ 
  CHECK_COND(IsFixedQuanto(), ITO33_BONDLIKE_NOT_FIXEDQUANTO);

  return m_dCorrelation; 
}

Date ConvertibleLike::GetMaturityDate() const
{
  return m_pBondLikeTerms->GetMaturityDate();
}

double ConvertibleLike::GetAccruedInterestValue() const
{
  CheckSessionData();

  if ( GetBondLikeTerms()->GetCashDistribution() )
    return GetBondLikeTerms()->GetCashDistribution()->
                        GetAccrued(GetSessionData()->GetValuationDate())
          * GetBondLikeTerms()->GetNominal();
  else
    return 0;
}

void ConvertibleLike::ValidateWith(const SessionData& sessionData) const
{
  Derivative::ValidateWith(sessionData);

  // Fixed FX is required for fixed quanto or trigger in currency of underlying
  if ( IsCrossCurrency(sessionData) )
    if (   GetTriggerInCurrencyOf() == TriggerInCurrencyOf_Underlying 
        || m_bIsFixedQuanto  )
    {
      CHECK_COND( GetFixedFXRate() > 0., ITO33_BONDLIKE_FIXEDFXRATE_NOT_SET );
    }
}

void ConvertibleLike::DumpMe(XML::Tag& tagParent) const
{
  Derivative::DumpMe(tagParent);
  
  tagParent.Element(*m_pBondLikeTerms);  

  tagParent.Element(XML_TAG_CONVERTIBLELIKE_NEWSHARE)(m_bConvertIntoNewShare);

  if ( m_triggerInCurrencyOf != TriggerInCurrencyOf_Derivative )
    tagParent.Element(XML_TAG_BONDLIKE_TRIGGERINCURRENCYOF)
              (
                GetNameFromEnumValue(
                  m_triggerInCurrencyOf,
                  SIZEOF(g_triggerInCurrencyOfs),
                  g_triggerInCurrencyOfs)
              );

  if ( m_dFixedFXRate > 0 )
    tagParent.Element(XML_TAG_BONDLIKE_FIXEDFXRATE)(m_dFixedFXRate);  
  
  if ( IsFixedQuanto() )
  {
    XML::Tag tagFixedQuanto(XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_ROOT,
                             tagParent);

    tagFixedQuanto.Element(XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_FXVOL)
                           (m_dFXRateVolatility);

    tagFixedQuanto.Element(XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_CORRELATION)
                           (m_dCorrelation);
    
  }

  if ( m_triggerAsPercentageOf != TriggerAsPercentageOf_Principal )
    tagParent.Element (XML_TAG_BONDLIKE_TAPO)
                      (
                        GetNameFromEnumValue(
                          m_triggerAsPercentageOf,
                          SIZEOF(g_triggerAsPercentageOfs),
                          g_triggerAsPercentageOfs)
                      );

  if ( IsExchangeable() )
  {
    XML::Tag tagExchangeable(XML_TAG_CONVERTIBLELIKE_EXCHANGEABLE_ROOT,
                             tagParent);

    tagExchangeable.Element(XML_TAG_CONVERTIBLELIKE_EXCHANGEABLE_COD)
                           (m_bExchangeableUponDefault);
  }

  DumpMarketPrice(tagParent);
}

} // namespace fianance

} // namespace ito33
