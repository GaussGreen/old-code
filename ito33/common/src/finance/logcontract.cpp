/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/logcontract.cpp
// Purpose:     Implementation of LogContract
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/logcontract.h"
#include "ito33/finance/derivative_visitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/logcontract.h"

namespace ito33
{

namespace finance
{


void LogContract::SetStartSharePrice(double dS0)
{ 
  ASSERT(dS0 > 0);

  m_dS0 = dS0;
}

LogContract::LogContract(Date maturityDate, Date startOfReturnPeriod)
                       : m_maturityDate(maturityDate),
                         m_startOfReturnPeriod(startOfReturnPeriod),
                         m_dS0(-1)
{
}

void LogContract::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnLogContract(*this);
}

XML::Tag LogContract::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagLogContract(XML_TAG_LOGCONTRACT_ROOT, tagParent);

  DumpMe(tagLogContract);

  tagLogContract.Element(XML_TAG_FINANCE_MATURITY)(m_maturityDate);

  tagLogContract.Element(XML_TAG_LOGCONTRACT_T0)(m_startOfReturnPeriod);

  if ( m_dS0 > 0 )
    tagLogContract.Element(XML_TAG_LOGCONTRACT_S0)(m_dS0);

  return tagLogContract;
}


} // namespace finance

} // namespace ito33
