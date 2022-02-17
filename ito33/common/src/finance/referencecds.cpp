/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/referencecds.cpp
// Purpose:     implement financial reference cds class
// Created:     2006/05/17
// RCS-ID:      $Id: referencecds.cpp,v 1.11 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/constants.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/referencecds.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/referencecds.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"

extern const ito33::finance::Error
  ITO33_REFERENCECDS_MATURITY,
  ITO33_REFERENCECDS_INVALID_SPREAD,
  ITO33_REFERENCECDS_NONZERO_PRICE,
  ITO33_REFERENCECDS_MATURITYDATE_NOT_AVAILABLE,
  ITO33_REFERENCECDS_FIRSTPAYMENTDATE_NOT_AVAILABLE,
  ITO33_REFERENCECDS_SPREAD_NOT_AVAILABLE,
  ITO33_REFERENCECDS_SPREADSTREAM_NOT_AVAILABLE,
  ITO33_REFERENCECDS_SPREADSTREAM_UNDEFSESSIONDATA;

//-----------------------------------------------------------------------------
// helper function
//-----------------------------------------------------------------------------
namespace ito33
{

/**
    Calculates the schedule of the spread payments and returns the first
    payment date.

    @param issueDate the issue date of reference date.
    @param nMonths maturity as number of months
    @param[out] maturityDate the maturity date of the CDS
    @return the first payment date of the CDS.
 */
Date CalculateReferenceCDSSchedule(const Date& issueDate,
                                   size_t nMonths,
                                   Date& maturityDate)
{
  Date firstPaymentDate;

  // Determine the maturity of the reference CDS
  // The year of today (issue date should be the valuation date)
  size_t nYear = issueDate.GetYear();

  // Find the next regular date: March 20, June 20, Sep 20, Dec 20
  // If valuation is a regular date, use the next regular date for the first
  // payment
  if ( issueDate >= Date(nYear, Date::Dec, 20) )
    firstPaymentDate = Date(nYear + 1, Date::Mar, 20);
  else
    for (size_t n = 0; n < 4; n++)
      if ( Date(nYear, Date::Month(3 * ( n + 1)), 20) > issueDate )
      {
        firstPaymentDate = Date(nYear, Date::Month(3 * ( n + 1)), 20);

        break;
      }

  // If there is too few days between the issue date and the current regular
  // date, we may need to use the regular date next to the current one as 
  // first payment date
  // For now, this is not implemented because of missing of specif
  // This brings up also the problem of maturity: Is it from the current 
  // regular date or the first payment date?

  maturityDate = firstPaymentDate;
  maturityDate.AddMonths(nMonths);

  return firstPaymentDate;
}

}

namespace ito33
{

namespace finance
{

ReferenceCDS::ReferenceCDS(size_t maturity,
                           Frequency freq, Date::DayCountConvention dcc,
                           double dRecoveryRate)
  : CDSLike(dRecoveryRate),
    m_nMonthsToMaturity(maturity),
    m_dSpread(-1.0),
    m_freq(freq),
    m_dcc(dcc)
{
  CHECK_COND(maturity > 0, ITO33_REFERENCECDS_MATURITY);

  // The price must be 0 at issue
  SetMarketPrice(0.0);
}

void ReferenceCDS::MakeSpreadStream()
{
  // do nothing if we are not ready
  if ( !(HasSpread() && m_pSessionData) )
    return;

  Date issueDate(m_pSessionData->GetValuationDate());
  Date firstPaymentDate;
  Date maturityDate;

  firstPaymentDate = CalculateReferenceCDSSchedule(issueDate,
                                                   m_nMonthsToMaturity,
                                                   maturityDate);

  m_pSpreadStream =  shared_ptr<CashFlowStreamUniform>
                      ( new CashFlowStreamUniform
                            (
                              issueDate, firstPaymentDate, maturityDate,
                              m_dSpread, m_dcc, m_freq
                            )
                      );
}

double ReferenceCDS::GetSpread() const
{
  // If the spread was set, the spread stream will have been created.
  CHECK_COND(HasSpread(), ITO33_REFERENCECDS_SPREAD_NOT_AVAILABLE);

  return m_dSpread;
}


shared_ptr<CashFlowStreamUniform> ReferenceCDS::GetSpreadStream() const
{ 
  CHECK_COND(HasSpread(), ITO33_REFERENCECDS_SPREADSTREAM_NOT_AVAILABLE);

  CHECK_COND(m_pSessionData, ITO33_REFERENCECDS_SPREADSTREAM_UNDEFSESSIONDATA);

  ASSERT_MSG( m_pSpreadStream,
              "m_pSpreadStream should have already been calculated.");

  return m_pSpreadStream;
}


void ReferenceCDS::SetSpread(double dSpread)
{
  CHECK_COND(dSpread > 0, ITO33_REFERENCECDS_INVALID_SPREAD);

  m_dSpread = dSpread;

  MakeSpreadStream();
}


void ReferenceCDS::SetSessionData(const shared_ptr<SessionData>& pSessionData)
{
  CDSLike::SetSessionData(pSessionData);

  MakeSpreadStream();
}


Date ReferenceCDS::GetMaturityDate() const
{
  CHECK_COND(m_pSessionData, ITO33_REFERENCECDS_MATURITYDATE_NOT_AVAILABLE);

  Date maturityDate;

  CalculateReferenceCDSSchedule(m_pSessionData->GetValuationDate(),
                                m_nMonthsToMaturity,
                                maturityDate);

  return maturityDate;
}


Date ReferenceCDS::GetFirstPaymentDate() const
{
  CHECK_COND(m_pSessionData, ITO33_REFERENCECDS_FIRSTPAYMENTDATE_NOT_AVAILABLE);

  Date maturityDate;

  return CalculateReferenceCDSSchedule(m_pSessionData->GetValuationDate(),
                                       m_nMonthsToMaturity,
                                       maturityDate);
}


bool ReferenceCDS::HasSpread() const
{
  return ( m_dSpread > 0 );
}


void ReferenceCDS::SetMarketPrice(double dPrice)
{
  CHECK_COND(dPrice == 0.0, ITO33_REFERENCECDS_NONZERO_PRICE);

  DoSetMarketPrice(dPrice);
}


void ReferenceCDS::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnReferenceCDS(*this);
}


void ReferenceCDS::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnReferenceCDS(*this);
}


XML::Tag ReferenceCDS::Dump(ito33::XML::Tag& tagParent) const 
{
  XML::Tag tagRefCDS(XML_TAG_REFERENCECDS_ROOT, tagParent);
	
  DumpMe(tagRefCDS);

  // dump the parameters need to construct this class
  tagRefCDS.Element(XML_TAG_REFERENCECDS_MATURITY_MONTHS)(m_nMonthsToMaturity);

  tagRefCDS.Element(XML_TAG_DAYCOUNTCONVENTION)
                 (XML::GetNameOfDayCountConvention(GetDayCountConvention()));

  tagRefCDS.Element(XML_TAG_PAYMENTFREQUENCY)
                 (XML::GetNameOfFrequency(GetPaymentFrequency()));
  
  // base class parameter
  tagRefCDS.Element(XML_TAG_FINANCE_RECOVERYRATE)(m_dRecoveryRate); 

  // dump the spread, if set
  if ( HasSpread() )
    tagRefCDS.Element(XML_TAG_REFERENCECDS_SPREAD)(m_dSpread);

  // dump the market price (defaults to zero in constructor)
  DumpMarketPrice(tagRefCDS);

  return tagRefCDS;
}


} // namespace finance

} // namespace ito33
