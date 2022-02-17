/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/varianceswap.cpp
// Purpose:     Implementation of variance swap class
// Created:     2006/02/21
// RCS-ID:      $Id: conditionalvarianceswap.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/conditionalvarianceswap.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/varianceswap.h"

extern const ito33::finance::Error
  ITO33_VARIANCESWAP_INVALID_CURRENT_COUNT;

namespace ito33
{

namespace finance
{

ConditionalVarianceSwap::ConditionalVarianceSwap(
    const shared_ptr<VarianceSwapTerms>& pTerms,
    double dVolatilityStrike)
  : VarianceSwapLike(pTerms, dVolatilityStrike),
    m_iCurrentConditionalCount(-1)
{
}


void ConditionalVarianceSwap::SetCurrentValues(
  double dCurrentConditionalVolatility,
  size_t nNbSamplesUsed,
  int iCurrentConditionalCount)
{
  VarianceSwapLike::SetCurrentValues
                    ( dCurrentConditionalVolatility, nNbSamplesUsed );

  CHECK_COND( iCurrentConditionalCount >= 0
              && (size_t)iCurrentConditionalCount < m_pTerms->GetNbSamplingReturns(),
              ITO33_VARIANCESWAP_INVALID_CURRENT_COUNT);

  m_iCurrentConditionalCount = iCurrentConditionalCount;
}

void ConditionalVarianceSwap::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnConditionalVarianceSwap(*this);
}


XML::Tag ConditionalVarianceSwap::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // conditional variance swap tag
  XML::Tag tagSwap(XML_TAG_CONDITIONALVARIANCESWAP_ROOT, tagParent);

  VarianceSwapLike::DumpMe(tagSwap);

  // the current values are optional
  if ( HasCurrentValues() )
  {
    tagSwap.Element(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY)
                   (m_dCurrentVolatility);
    tagSwap.Element(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED)(m_nNbSamplesUsed);
    tagSwap.Element(XML_TAG_VARIANCESWAP_CURRENT_COUNT)
                   (m_iCurrentConditionalCount);
  }

  return tagSwap;
}


} // namespace finance

} // namespace ito33
