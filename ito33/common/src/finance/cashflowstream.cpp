/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/cashflowstream.cpp
// Purpose:     implementation of non inline methods of finance::CashFlowStream
// Author:      Pedro Ferreira
// Created:     Mar 24, 2004
// RCS-ID:      $Id: cashflowstream.cpp,v 1.27 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"

extern const ito33::finance::Error
  ITO33_INVALID_FREQUENCY;

namespace ito33
{

namespace finance
{

CashFlowStream::CashFlowStream
    (
      Date contractingDate, 
      Date::DayCountConvention dcc,
      Frequency freq
    ) 
    : m_contractingDate(contractingDate), m_dcc(dcc), m_freq(freq)
{
  CHECK_COND( IsValid(m_freq), ITO33_INVALID_FREQUENCY );
}

double CashFlowStream::GetAccrued(Date date) const
{
  if ( date <= GetContractingDate() || date >= GetLastPaymentDate() )
    return 0;

  Elements::const_iterator iter = m_elements.upper_bound(date);
  
  Date nextPaymentDate = iter->first;
  double dNextAmount = iter->second;
  Date previousPaymentDate = (iter == m_elements.begin()
                           ? GetContractingDate() : (--iter)->first);

  return dNextAmount 
         * Date::DaysDiffWithDayCount(previousPaymentDate, date, m_dcc)
         / Date::DaysDiffWithDayCount(previousPaymentDate, nextPaymentDate, 
                                      m_dcc);
}

double CashFlowStream::GetDiscount(const YieldCurve& yc, Date date) const
{
  double dDiscount = 0.;
  Elements::const_iterator iter = m_elements.upper_bound(date);
  for ( ; iter != m_elements.end(); ++iter)
    dDiscount += yc.GetForwardDiscountFactor(date, iter->first) * iter->second;

  return dDiscount;
}

const CashFlowStream::CollectionElements& CashFlowStream::GetAll() const
{
  if ( m_collectionElements.empty() )
    DoGetAll();

  return m_collectionElements;
}

void CashFlowStream::DoGetAll() const
{
  Elements::const_iterator iter;
  CashFlowStream *pCF = const_cast<CashFlowStream *>(this);

  for (iter = m_elements.begin(); iter != m_elements.end(); ++iter)
    pCF->m_collectionElements.push_back
                              ( CashFlow(iter->first, iter->second) );
}

void CashFlowStream::DumpMe(XML::Tag& tagRoot) const
{

  tagRoot.Element(XML_TAG_CASHFLOWSTREAM_CONTRACTINGDATE)
                 ( GetContractingDate() );
  
  tagRoot.Element(XML_TAG_DAYCOUNTCONVENTION)
                 (XML::GetNameOfDayCountConvention(GetDayCountConvention()));
  
  tagRoot.Element(XML_TAG_PAYMENTFREQUENCY)
                 (XML::GetNameOfFrequency(GetPaymentFrequency()));
}


} // namespace finance

} // namespace ito33
