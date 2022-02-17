/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/parbond.cpp
// Purpose:     Implement financial ParBond class
// Created:     2005/05/20
// autor:       ZHANG Yunzhi
// RCS-ID:      $Id: parbond.cpp,v 1.6 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/parbond.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/parbond.h"

#include "ito33/finance/derivative_visitor.h"

extern const ito33::finance::Error ITO33_PARBOND_INVALID_CONTRACTING_DATE,
                          ITO33_PARBOND_NULL_MATURITY,
                          ITO33_PARBOND_NEGATIVE_YTM,
                          ITO33_PARBOND_NEGATIVE_SPREAD,
                          ITO33_PARBOND_INCONSISTENT_FREQUENCY_MATURITY,
                          ITO33_INVALID_RECOVERYRATE_1;

namespace ito33
{

namespace finance
{

void ParBond::SelfValidate() const
{
  CHECK_COND(m_fakeContractingDate.IsValid(),
             ITO33_PARBOND_INVALID_CONTRACTING_DATE);

  finance::Validate(m_paymentFrequency);

  ito33::Validate(m_dcc);

  CHECK_COND(m_nMaturity > 0, ITO33_PARBOND_NULL_MATURITY);

  CHECK_COND(m_dYTM >= 0, ITO33_PARBOND_NEGATIVE_YTM);

  CHECK_COND(m_dSpread >= 0, ITO33_PARBOND_NEGATIVE_SPREAD);

  CHECK_COND_1(m_dRecoveryRate >= 0 && m_dRecoveryRate <= 1,
               ITO33_INVALID_RECOVERYRATE_1,
               m_dRecoveryRate);

  CHECK_COND( (m_nMaturity * m_paymentFrequency) % 12 == 0,
              ITO33_PARBOND_INCONSISTENT_FREQUENCY_MATURITY);
}


Date ParBond::GetContractingDate() const
{
  if(m_pSessionData)
    return m_pSessionData->GetValuationDate();
  else
    return m_fakeContractingDate;
}

Date ParBond::GetMaturityDate() const
{
  return AddMonthsAdjustedForEndOfMonth
                    (GetContractingDate(), m_nMaturity);
}

double ParBond::GetCouponRate() const
{
  return GetCashFlowStream()->begin()->second * m_paymentFrequency;  
}

shared_ptr<CashFlowStream> ParBond::GetCashFlowStream() const
{
  // GetSessionData() validate session();
  Date dtContracting = GetContractingDate();

  size_t nNbPayments = m_nMaturity * m_paymentFrequency / 12 ;
  std::vector<Date> pPaymentDates(nNbPayments);

  size_t nFrac = 12 / m_paymentFrequency;
  size_t paymentDuration = nFrac;
  size_t nIdx = 0;

  // dSum = Sum_all_coupon ( 1+ (ytm+spread)/frequency)^-(frequency*DeltaT) )
  double
    dSum = 0,
    dBase = pow( 1 + (m_dYTM + m_dSpread) / m_paymentFrequency,
                 -m_paymentFrequency),
    dLastDiscount = 0; 

  while(paymentDuration <= m_nMaturity)
  {
    pPaymentDates[nIdx] = AddMonthsAdjustedForEndOfMonth
                            (dtContracting, paymentDuration);
    dLastDiscount = 
       pow(dBase, Date::YearsDiff(dtContracting, pPaymentDates[nIdx++], m_dcc));

    dSum += dLastDiscount;
    paymentDuration += nFrac;
  }

  shared_ptr<CashFlowStream> pcfs = 
    shared_ptr<CashFlowStream> 
      (new CashFlowStreamGeneral
            (
              dtContracting,
              pPaymentDates,
              std::vector<double>(nNbPayments, (1 - dLastDiscount) / dSum),
              m_dcc,
              m_paymentFrequency
            )
       );
  return pcfs;
}


void ParBond::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnParBond(*this);
}

XML::Tag ParBond::Dump(ito33::XML::Tag& tagParent) const 
{
  XML::Tag tagParBond(XML_TAG_PARBOND_ROOT,tagParent);
	
  DumpMe(tagParBond);

  tagParBond.Element(XML_TAG_PARBOND_CONTRACTING_DATE)
                    (GetContractingDate());
  tagParBond.Element(XML_TAG_PARBOND_YTM)(m_dYTM);

  tagParBond.Element(XML_TAG_PARBOND_SPREAD)(m_dSpread);
  tagParBond.Element(XML_TAG_PARBOND_MATURITY)(m_nMaturity);
  tagParBond.Element(XML_TAG_PAYMENTFREQUENCY)
            (XML::GetNameOfFrequency(m_paymentFrequency));
  tagParBond.Element(XML_TAG_DAYCOUNTCONVENTION)
            (XML::GetNameOfDayCountConvention(m_dcc));
  tagParBond.Element(XML_TAG_FINANCE_RECOVERYRATE)(m_dRecoveryRate);

  return tagParBond;
}


} // namespace finance

} // namespace ito33
