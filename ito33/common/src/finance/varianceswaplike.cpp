/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/varianceswaplike.cpp
// Purpose:     Implementation of the base variance swap class
// Created:     2006/08/04
// RCS-ID:      $Id: varianceswaplike.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswaplike.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/varianceswap.h"

extern const ito33::finance::Error 
  ITO33_VARIANCESWAP_NULL_TERMS,
  ITO33_VARIANCESWAP_NEGATIVE_VOLATILITYSTRIKE,
  ITO33_VARIANCESWAP_NEGATIVE_CURRENTVOLATILITY,
  ITO33_VARIANCESWAP_INVALID_NB_SAMPLES_USED,
  ITO33_VARIANCESWAP_ZERO_NB_SAMPLES_USED,
  ITO33_PREVIOUS_SHARE_PRICE_UNDEFINED,
  ITO33_VARIANCESWAP_NOCURRENTVOLATILITY;

namespace ito33
{

namespace finance
{

VarianceSwapLike::VarianceSwapLike(const shared_ptr<VarianceSwapTerms>& pTerms,
                                   double dVolatilityStrike)
  : m_pTerms(pTerms),
    m_dVolatilityStrike(dVolatilityStrike),
    m_dCurrentVolatility(-1.0),
    m_nNbSamplesUsed(0)
{
  // standard input checks
  CHECK_COND( m_pTerms, ITO33_VARIANCESWAP_NULL_TERMS );

  CHECK_COND( m_dVolatilityStrike >= 0.0, 
              ITO33_VARIANCESWAP_NEGATIVE_VOLATILITYSTRIKE );
}

void VarianceSwapLike::SetCurrentValues(double dCurrentVolatility,
                                        size_t nNbSamplesUsed)
{
  CHECK_COND( dCurrentVolatility >= 0.0, 
              ITO33_VARIANCESWAP_NEGATIVE_CURRENTVOLATILITY);

  CHECK_COND( nNbSamplesUsed < m_pTerms->GetNbSamplingReturns(),
              ITO33_VARIANCESWAP_INVALID_NB_SAMPLES_USED);

  CHECK_COND( nNbSamplesUsed > 0,
              ITO33_VARIANCESWAP_ZERO_NB_SAMPLES_USED);

  m_dCurrentVolatility = dCurrentVolatility;

  m_nNbSamplesUsed = nNbSamplesUsed;
}

Date VarianceSwapLike::GetMaturityDate() const
{
  return m_pTerms->GetMaturityDate();
}


void VarianceSwapLike::ValidateWith(const SessionData& sessionData) const
{
  // validation depends on relation of valuation date to sampling start date
  Date valuationDate = sessionData.GetValuationDate();
  Date startOfSamplingPeriod = m_pTerms->GetStartOfSamplingPeriod();

  if ( valuationDate >= startOfSamplingPeriod )
  {
    // The previous share price must be set
    double dPreviousSharePrice = sessionData.GetPreviousSharePrice();
    CHECK_COND(dPreviousSharePrice > 0., ITO33_PREVIOUS_SHARE_PRICE_UNDEFINED);
  }

  if ( valuationDate > startOfSamplingPeriod )
    CHECK_COND( HasCurrentValues(), ITO33_VARIANCESWAP_NOCURRENTVOLATILITY);
}


void VarianceSwapLike::DumpMe(XML::Tag& tagParent) const
{
  // dump values common to all derivatives
  Derivative::DumpMe(tagParent);

  tagParent.Element(XML_TAG_VARIANCESWAP_VOLATILITYSTRIKE)(m_dVolatilityStrike);
  
  m_pTerms->DumpMe(tagParent);

  // dump the market price, if set
  DumpMarketPrice(tagParent);  
}


} // namespace finance

} // namespace ito33
