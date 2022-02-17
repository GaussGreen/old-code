/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/option.cpp
// Purpose:     Implementation of shared ptr for Option class
// Author:      ZHANG Yunzhi
// Created:     April 28, 2003
// RCS-ID:      $Id: option.cpp,v 1.32 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/option.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"
#include "ito33/finance/optionerror.h"
#include "ito33/finance/impliedvolpriceconverter.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/option.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/optiontype.h"

extern const ito33::finance::OptionError
  ITO33_OPTION_NEGATIVE_STRIKE,
  ITO33_OPTION_NEGATIVE_IMPLIED_VOLATILITY,
  ITO33_OPTION_NO_IMPLIED_VOLATILITY;

extern const ito33::finance::Error ITO33_INVALID_MARKETPRICE_NEGATIVE;

namespace ito33
{

namespace finance
{


Option::Option(double dStrike, Date maturityDate,
               OptionType optionType, ExerciseType exerciseType)
               :OptionLike(maturityDate, optionType, exerciseType),
               m_dStrike(dStrike), m_dImpliedVol(-1.0)

{
  CHECK_COND( m_dStrike > 0, ITO33_OPTION_NEGATIVE_STRIKE);
}

void Option::SetImpliedVol(double dImpliedVol)
{
  CHECK_COND( dImpliedVol >= 0, ITO33_OPTION_NEGATIVE_IMPLIED_VOLATILITY);

  m_dImpliedVol = dImpliedVol;
  UnsetMarketPrice(); // invalidate market price
}

void Option::SetMarketPrice(double dPrice)
{
  Derivative::SetMarketPrice(dPrice);
  m_dImpliedVol = -1.; // invalidate implied volatility
}

double Option::DoGetMarketPrice() const
{
  if ( IsImpliedVolSet() )    
    return GetPriceFromImpliedVol(*this, m_dImpliedVol);
  else
    return Derivative::DoGetMarketPrice();
}

bool Option::HasMarketPrice() const
{
  return IsImpliedVolSet() || Derivative::HasMarketPrice();
}

double Option::GetImpliedVolFrom(double dPrice) const
{
  CHECK_COND ( dPrice > 0.0, ITO33_INVALID_MARKETPRICE_NEGATIVE );

  return GetImpliedVolFromPrice(*this, dPrice);
}

double Option::GetPriceFrom(double dVol) const
{
  CHECK_COND( dVol >= 0, ITO33_OPTION_NEGATIVE_IMPLIED_VOLATILITY);

  return GetPriceFromImpliedVol(*this, dVol);
}

double Option::GetVegaFrom(double dVol) const
{
  CHECK_COND( dVol >= 0, ITO33_OPTION_NEGATIVE_IMPLIED_VOLATILITY);

  return GetVegaFromImpliedVol(*this, dVol);
}

double Option::GetImpliedVol() const
{
  if ( IsImpliedVolSet() )
    return m_dImpliedVol;
  else if ( HasMarketPrice() )
    return GetImpliedVolFrom(m_dPrice);
  else
    throw EXCEPTION(ITO33_OPTION_NO_IMPLIED_VOLATILITY);
}

void Option::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnOption(*this);
}

void Option::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnOption(*this);
}

XML::Tag Option::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // option tag
  XML::Tag tagOption(XML_TAG_OPTION_ROOT, tagParent);
  
  DumpMe(tagOption);

  tagOption.Element(XML_TAG_FINANCE_STRIKE)(m_dStrike);

  if ( IsImpliedVolSet() )
    tagOption.Element(XML_TAG_OPTION_IMPLIED_VOL)(m_dImpliedVol);

  return tagOption;
}

} // namespace finance

} // namespace ito33
