/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/gammavarianceswap.cpp
// Purpose:     Implementation of the gamma variance swap class
// Created:     2006/08/04
// RCS-ID:      $Id: gammavarianceswap.cpp,v 1.4 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/gammavarianceswap.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/varianceswap.h"

extern const ito33::finance::Error 
  ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_START_PRICE;

namespace ito33
{

namespace finance
{

GammaVarianceSwap::GammaVarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
                                     double dVolatilityStrike)
  : VarianceSwapLike(pTerms, dVolatilityStrike),
    m_dStartSharePrice(0)
{
}


void GammaVarianceSwap::SetCurrentValues(double dCurrentGammaVolatility,
                                         size_t nNbSamplesUsed,
                                         double dStartSharePrice)
{
  VarianceSwapLike::SetCurrentValues
                    ( dCurrentGammaVolatility, nNbSamplesUsed );

  CHECK_COND( dStartSharePrice > 0.0,
              ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_START_PRICE);

  m_dStartSharePrice = dStartSharePrice;
}


void GammaVarianceSwap::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnGammaVarianceSwap(*this);
}


XML::Tag GammaVarianceSwap::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // gamma variance swap tag
  XML::Tag tagSwap(XML_TAG_GAMMAVARIANCESWAP_ROOT, tagParent);

  VarianceSwapLike::DumpMe(tagSwap);

  // the current values are all optional
  if ( HasCurrentValues() )
  {
    tagSwap.Element(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY)
                   (m_dCurrentVolatility);
    tagSwap.Element(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED)(m_nNbSamplesUsed);
    tagSwap.Element(XML_TAG_VARIANCESWAP_START_SHARE_PRICE)
                   (m_dStartSharePrice);
  }

  return tagSwap;
}


} // namespace finance

} // namespace ito33
