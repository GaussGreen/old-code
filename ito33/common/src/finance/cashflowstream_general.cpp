/////////////////////////////////////////////////////////////////////////////
// Name:        finance/cashflowstream_general.cpp
// Purpose:     implementation of non inline methods of 
//              finance::CashFlowStreamGeneral
// Author:      Pedro Ferreira
// Created:     Mar 25, 2004
// RCS-ID:      $Id: cashflowstream_general.cpp,v 1.20 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/arraycheckers.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream_general.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/cashflowstream_all.h"

extern const ito33::finance::Error ITO33_CASHFLOWSTREAM_BAD_FIRST_PAYMENT_DATE;

namespace ito33
{

namespace finance
{

CashFlowStreamGeneral::CashFlowStreamGeneral
                       (
                         Date contractingDate,
                         const std::vector<Date>& paymentDates,
                         const std::vector<double>& paymentRates,
                         Date::DayCountConvention dcc,
                         Frequency freq
                       )
                       : CashFlowStream(contractingDate, dcc, freq) 
{
  CheckIncreasingOrder(paymentDates, "Cash flow payment dates");

  CheckNonNegativity(paymentRates, "Cash flow payment rates");
 
  CHECK_COND
  (
    contractingDate < paymentDates[0],
    ITO33_CASHFLOWSTREAM_BAD_FIRST_PAYMENT_DATE
  );

  for (size_t n = 0; n < paymentDates.size(); n++)
    AddCashFlow(paymentDates[n], paymentRates[n]);
}

ito33::XML::Tag 
CashFlowStreamGeneral::Dump(const char* name, ito33::XML::Tag& tagParent) const
{
  XML::Tag tagName(name, tagParent);
  XML::Tag tagMe(XML_TAG_CASHFLOWSTREAMGENERAL_ROOT, tagName);

  CashFlowStream::DumpMe(tagMe);

  Elements::const_iterator iter;

  for (iter = m_elements.begin(); iter != m_elements.end(); ++iter)
  {
    XML::Tag tagElement(XML_TAG_CASHFLOWSTREAMGENERAL_CASHFLOW, tagMe);

    tagElement.Element(XML_TAG_FINANCE_DATE)(iter->first);
    tagElement.Element(XML_TAG_FINANCE_RATE)(iter->second);
  }

  return tagName;
}


} // namespace finance

} // namespace ito33
