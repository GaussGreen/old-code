/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/varianceswaption.cpp
// Purpose:     Implementation of variance swaption class
// Created:     2006/07/21
// RCS-ID:      $Id: varianceswaption.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/optionerror.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswaption.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/optiontype.h"
#include "ito33/xml/finance/varianceswap.h"

extern const ito33::finance::Error 
  ITO33_VARIANCESWAPTION_NULL_VSTERMS,
  ITO33_VARIANCESWAPTION_MATURITY_AFTER_STARTOFVSSAMPLING;

extern const ito33::finance::OptionError
  ITO33_OPTION_INVALID_OPTIONTYPE,
  ITO33_OPTION_NEGATIVE_STRIKE;

namespace ito33
{

namespace finance
{

VarianceSwaption::VarianceSwaption
   (const shared_ptr<VarianceSwapTerms>& pTerms,
    OptionType optionType, double dStrike, Date maturityDate)
  : m_pTerms(pTerms),
    m_optionType(optionType),
    m_dStrike(dStrike),
    m_maturityDate(maturityDate)
{
  // standard input checks
  CHECK_COND( m_pTerms, ITO33_VARIANCESWAPTION_NULL_VSTERMS );

  CHECK_COND( IsValidOptionType(m_optionType) , 
              ITO33_OPTION_INVALID_OPTIONTYPE );

  CHECK_COND( m_dStrike >= 0.0, 
              ITO33_OPTION_NEGATIVE_STRIKE );

  CHECK_COND( m_pTerms->GetStartOfSamplingPeriod() >= m_maturityDate,
              ITO33_VARIANCESWAPTION_MATURITY_AFTER_STARTOFVSSAMPLING );
}

void VarianceSwaption::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnVarianceSwaption(*this);
}

XML::Tag VarianceSwaption::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // variance swap tag
  XML::Tag tagSwaption(XML_TAG_VARIANCESWAPTION_ROOT, tagParent);

  DumpMe(tagSwaption);

  tagSwaption.Element(XML_TAG_OPTION_TYPE)
                     (
                      GetNameFromEnumValue(
                        m_optionType,
                        SIZEOF(g_optionTypes),
                        g_optionTypes)
                     );

  tagSwaption.Element(XML_TAG_FINANCE_STRIKE)(m_dStrike);

  tagSwaption.Element(XML_TAG_FINANCE_MATURITY)(m_maturityDate);

  {
    XML::Tag tagTerms(XML_TAG_VARIANCESWAPTION_VSTERMS, tagSwaption);

    m_pTerms->DumpMe(tagTerms);
  }

  // dump the market price, if set
  DumpMarketPrice(tagSwaption);

  return tagSwaption;
}


} // namespace finance

} // namespace ito33
