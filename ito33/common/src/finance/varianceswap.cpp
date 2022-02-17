/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/varianceswap.cpp
// Purpose:     Implementation of variance swap class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswap.cpp,v 1.22 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/varianceswap.h"

namespace ito33
{

namespace finance
{

VarianceSwap::VarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
                           double dVolatilityStrike)
  : VarianceSwapLike(pTerms, dVolatilityStrike)
{
}


void VarianceSwap::SetCurrentValues(double dCurrentVolatility,
                                    size_t nNbSamplesUsed)
{
  VarianceSwapLike::SetCurrentValues(dCurrentVolatility, nNbSamplesUsed);
}

void VarianceSwap::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnVarianceSwap(*this);
}


XML::Tag VarianceSwap::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // variance swap tag
  XML::Tag tagSwap(XML_TAG_VARIANCESWAP_ROOT, tagParent);

  VarianceSwapLike::DumpMe(tagSwap);

  // the current volatility value and number of samples used are optional
  if ( HasCurrentValues() )
  {
    tagSwap.Element(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY)
                   (m_dCurrentVolatility);
    tagSwap.Element(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED)(m_nNbSamplesUsed);
  }

  return tagSwap;
}


} // namespace finance

} // namespace ito33
