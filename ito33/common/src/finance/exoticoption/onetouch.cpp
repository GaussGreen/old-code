/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/exoticoption/onetouch.cpp
// Purpose:     Implement financial OneTouch class
// Created:     2005/07/04
// RCS-ID:      $Id: onetouch.cpp,v 1.11 2006/08/19 22:40:30 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/optionerror.h"
#include "ito33/finance/exoticoption/onetouch.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/exoticoption/onetouch.h"

#include "ito33/finance/derivative_visitor.h"

extern const ito33::finance::Error ITO33_INVALID_BARRIER;

namespace ito33
{

namespace finance
{


OneTouch::OneTouch
(Date maturityDate, double dBarrier, BarrierType barrierType, 
 RebateType rebateType)
  : Derivative(), m_maturityDate(maturityDate), 
    m_dBarrier(dBarrier), m_barrierType(barrierType),
    m_rebateType(rebateType)
{
  CHECK_COND(dBarrier > 0, ITO33_INVALID_BARRIER);
}

void OneTouch::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnOneTouch(*this);
}

XML::Tag OneTouch::Dump(ito33::XML::Tag& tagParent) const 
{
  XML::Tag tagOneTouch(XML_TAG_ONETOUCH_ROOT, tagParent);
	
  DumpMe(tagOneTouch);

  tagOneTouch.Element(XML_TAG_FINANCE_MATURITY)(m_maturityDate);

  tagOneTouch.Element(XML_TAG_BARRIER)(m_dBarrier);

  tagOneTouch.Element(XML_TAG_BARRIER_TYPE)
                     (
                       GetNameFromEnumValue
                       (m_barrierType, SIZEOF(g_barrierTypes), g_barrierTypes)
                     );

  tagOneTouch.Element(XML_TAG_REBATE_TYPE)
                     ( 
                       GetNameFromEnumValue
                       (m_rebateType, SIZEOF(g_rebateTypes), g_rebateTypes)
                     );

  DumpMarketPrice(tagOneTouch);

  return tagOneTouch;
}

} // namespace finance

} // namespace ito33
