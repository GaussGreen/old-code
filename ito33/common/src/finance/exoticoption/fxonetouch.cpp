/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/exoticoption/fxonetouch.cpp
// Purpose:     Implement financial FX OneTouch class
// Created:     2006/01/31
// RCS-ID:      $Id: fxonetouch.cpp,v 1.6 2006/08/19 22:40:30 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/optionerror.h"
#include "ito33/finance/exoticoption/fxonetouch.h"
#include "ito33/finance/exoticoption/onetouchutils.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/exoticoption/onetouch.h"

#include "ito33/finance/derivative_visitor.h"

extern const ito33::finance::OptionError ITO33_ONETOUCH_NEGATIVE_BSBARRIER,
                                         ITO33_ONETOUCH_BARRIER_NOTFOUND,
                                         ITO33_ONETOUCH_USE_MARKET_QUOTE,
                                         ITO33_ONETOUCH_NO_MARKET_QUOTE,
                                         ITO33_ONETOUCH_QUOTE_TOO_LARGE;

namespace ito33
{

namespace finance
{

FXOneTouch::FXOneTouch
(Date maturityDate, double dBSBarrier, BarrierType barrierType, double dVol)
  : OneTouch(maturityDate), 
    m_dBSBarrier(dBSBarrier), 
    m_dReferenceVol(dVol)
{
  CHECK_COND(m_dBSBarrier >= 0, ITO33_ONETOUCH_NEGATIVE_BSBARRIER);

  // Set base class variables
  m_barrierType = barrierType;
  m_rebateType = Rebate_Immediate;

  m_dBarrier = -1;
}

double FXOneTouch::GetBarrier() const
{
  double dBarrier = GetBarrierFromDelta
                    ( m_pSessionData, m_maturityDate, 
                      m_dBSBarrier, m_barrierType, m_dReferenceVol );

  // Check if the barrier is found
  CHECK_COND( dBarrier != 0, ITO33_ONETOUCH_BARRIER_NOTFOUND );

  return dBarrier;
}

void FXOneTouch::SetMarketQuote(double dMarketQuote)
{
  m_dQuote = dMarketQuote;

  m_dPrice = m_dQuote + m_dBSBarrier;

  CHECK_COND( HasMarketPrice(),
              ITO33_ONETOUCH_QUOTE_TOO_LARGE );
}

double FXOneTouch::GetMarketQuote() const
{
  CHECK_COND( HasMarketPrice(), ITO33_ONETOUCH_NO_MARKET_QUOTE);
  return m_dQuote;
}

void FXOneTouch::SetMarketPrice(double /* dPrice */)
{
  // Must call SetMarketQuote to set the market price
  throw EXCEPTION(ITO33_ONETOUCH_USE_MARKET_QUOTE);
}

void FXOneTouch::SetBarrier(double /* dBarrier */)
{ 
  FAIL("Not implemented.");
}

void FXOneTouch::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnFXOneTouch(*this);
}

XML::Tag FXOneTouch::Dump(ito33::XML::Tag& tagParent) const 
{
  XML::Tag tagFXOneTouch(XML_TAG_FXONETOUCH_ROOT, tagParent);
	
  DumpMe(tagFXOneTouch);

  tagFXOneTouch.Element(XML_TAG_FINANCE_MATURITY)(m_maturityDate);

  tagFXOneTouch.Element(XML_TAG_FXONETOUCH_BSBARRIER)(m_dBSBarrier);

  tagFXOneTouch.Element(XML_TAG_BARRIER_TYPE)
                       (
                         GetNameFromEnumValue
                         (m_barrierType, SIZEOF(g_barrierTypes), g_barrierTypes)
                       );

  tagFXOneTouch.Element(XML_TAG_FXONETOUCH_REFVOL)(m_dReferenceVol);

  if ( HasMarketPrice() )
    tagFXOneTouch.Element(XML_TAG_FXONETOUCH_QUOTE)(m_dQuote);

  return tagFXOneTouch;
}

} // namespace finance

} // namespace ito33
