/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/optionvarianceswap.cpp
// Purpose:     Implementation of the option variance swap class
// Created:     2006/08/04
// RCS-ID:      $Id: optionvarianceswap.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/optionvarianceswap.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/optiontype.h"
#include "ito33/xml/finance/varianceswap.h"

extern const ito33::finance::Error 
  ITO33_VARIANCESWAP_INVALID_OPTIONTYPE;

namespace ito33
{

namespace finance
{

OptionVarianceSwap::OptionVarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
                                       double dVolatilityStrike,
                                       OptionType optionType)
  : VarianceSwapLike(pTerms, dVolatilityStrike),
    m_optionType(optionType) 
{
  CHECK_COND( IsValidOptionType(m_optionType) , 
              ITO33_VARIANCESWAP_INVALID_OPTIONTYPE );
}


void OptionVarianceSwap::SetCurrentValues(double dCurrentVolatility,
                                          size_t nNbSamplesUsed)
{
  VarianceSwapLike::SetCurrentValues(dCurrentVolatility, nNbSamplesUsed);
}


void OptionVarianceSwap::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnOptionVarianceSwap(*this);
}


XML::Tag OptionVarianceSwap::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // option variance swap tag
  XML::Tag tagSwap(XML_TAG_OPTIONVARIANCESWAP_ROOT, tagParent);

  VarianceSwapLike::DumpMe(tagSwap);

  tagSwap.Element(XML_TAG_OPTION_TYPE)
                     (
                      GetNameFromEnumValue(
                        m_optionType,
                        SIZEOF(g_optionTypes),
                        g_optionTypes)
                     );

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
