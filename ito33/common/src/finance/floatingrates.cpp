/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/floatingrates.cpp
// Purpose:     implementation of non inline methods of finance::FloatingRates
// Author:      Nabil
// Created:     2005/07/15
// RCS-ID:      $Id: floatingrates.cpp,v 1.11 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"

#include "ito33/arraycheckers.h"
#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/floatingrates.h"
#include "ito33/finance/rateutils.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/lastpaymenttype.h"
#include "ito33/xml/finance/floatingrates.h"

extern const ito33::Error 
              ITO33_BAD_PARAM;
extern const ito33::finance::Error
              ITO33_FLOATINGRATES_INCOHERENT_PAYMENTSTREAM,
              ITO33_FLOATINGRATES_INVALID_PAYMENTSTREAM,
              ITO33_FLOATINGRATES_NO_UNKNOWN_PAYMENTS,
              ITO33_CASHFLOWSTREAM_INCONSISTENT_DATES,
              ITO33_FLOATINGRATES_FIRSTKNOWNPAYMENT_BEFORE_STARTOFACCRUED,
              ITO33_FLOATINGRATES_FIRSTUNKNOWN_AFTER_LASTUNKNOWN;

namespace ito33
{

namespace finance
{

FloatingRates::FloatingRates(double dMargin, 
                             Date startOfAccruedDate, 
                             Date firstUnknownPaymentDate,
                             Date lastUnknownPaymentDate,
                             Frequency frequency,
                             LastPaymentType lastPaymentType )
              : m_dMargin(dMargin),
                m_bHasUnknownPayments(true),
                m_startOfAccruedDate(startOfAccruedDate), 
                m_firstUnknownPaymentDate(firstUnknownPaymentDate),
                m_lastUnknownPaymentDate(lastUnknownPaymentDate),
                m_paymentfrequency(frequency),
                m_lastPaymentType(lastPaymentType) 
{
  finance::Validate(lastPaymentType);

  Validate();

  Init();
}

FloatingRates::FloatingRates(double dMargin, 
                             Date startOfAccruedDate, 
                             Date firstUnknownPaymentDate,
                             Date lastUnknownPaymentDate,
                             Frequency frequency)
              : m_dMargin(dMargin),
                m_bHasUnknownPayments(true),
                m_startOfAccruedDate(startOfAccruedDate), 
                m_firstUnknownPaymentDate(firstUnknownPaymentDate),
                m_lastUnknownPaymentDate(lastUnknownPaymentDate),
                m_paymentfrequency(frequency),
                m_lastPaymentType(LastPaymentType_Max)
{
  Validate();

  Init();
}

FloatingRates::FloatingRates(Date startOfAccruedDate, Frequency frequency)
              : m_dMargin(0.),
                m_bHasUnknownPayments(false),
                m_startOfAccruedDate(startOfAccruedDate),
                m_paymentfrequency(frequency),
                m_lastPaymentType(LastPaymentType_Max)
{ 
  Validate();

  Init();
}

void FloatingRates::Validate() const
{
  if ( m_bHasUnknownPayments )
  {
    CHECK_COND_MSG
      ( m_firstUnknownPaymentDate > m_startOfAccruedDate, 
        ITO33_BAD_PARAM,
        "FloatingRates definition: first unknown payment date must be greater "
        "than the start of accrued date.");

    CHECK_COND(m_firstUnknownPaymentDate <= m_lastUnknownPaymentDate,
               ITO33_FLOATINGRATES_FIRSTUNKNOWN_AFTER_LASTUNKNOWN);
  }

  finance::Validate(m_paymentfrequency);
}

void FloatingRates::SetMultiplier(double dMultiplier)
{
  CHECK_COND_MSG( dMultiplier > 0., ITO33_BAD_PARAM,
      "FloatingRates definition: Multiplier value for the floating payment must be"
      " greater than zero.");
  
  m_dMultiplier = dMultiplier;
}

void FloatingRates::SetCap(double dCap)
{
  CHECK_COND_MSG( (dCap > 0. && dCap <= 1.), ITO33_BAD_PARAM,
      "FloatingRates definition: Cap value for the floating payment must be"
      " greater than zero and lower than 1.");
  
  m_dCap = dCap;
}

void FloatingRates::SetFloor(double dFloor)
{  
  CHECK_COND_MSG( dFloor >= 0., ITO33_BAD_PARAM,
      "FloatingRates definition: Floor value for the floating payment must be"
      " positive.");

  m_dFloor = dFloor;
}

void FloatingRates::SetFixingDelay(int iFixingDelay)
{  
  m_iFixingDelay = iFixingDelay;
}

void FloatingRates::SetDayCountConvention(Date::DayCountConvention dcc)
{
  m_dcc = dcc;
}

void FloatingRates::SetKnownPaymentStream
     ( const std::vector<Date>& paymentDates, 
       const std::vector<double>& paymentRates )
{  
  CheckIncreasingOrder
  (
    paymentDates, 
    "paymentDates must be a non empty array of increasing dates."
  );

  CheckNonNegativity
  (
    paymentRates, 
    "paymentRates must be an array of non negative rates."
  );

  CHECK_COND( paymentDates.size() == paymentRates.size(), 
              ITO33_FLOATINGRATES_INCOHERENT_PAYMENTSTREAM ); 
  
  CHECK_COND( paymentDates.front() > m_startOfAccruedDate, 
              ITO33_FLOATINGRATES_FIRSTKNOWNPAYMENT_BEFORE_STARTOFACCRUED );
  
  if ( m_bHasUnknownPayments )
    CHECK_COND( paymentDates.back() < m_firstUnknownPaymentDate, 
                ITO33_FLOATINGRATES_INVALID_PAYMENTSTREAM ); 

  m_pKnownPaymentRatesDates = paymentDates;
  m_pKnownPaymentRates = paymentRates;
}

Date FloatingRates::GetFirstUnknownPaymentDate() const
{
  CHECK_COND( m_bHasUnknownPayments, 
              ITO33_FLOATINGRATES_NO_UNKNOWN_PAYMENTS );
  
  return m_firstUnknownPaymentDate;
}

Date FloatingRates::GetLastUnknownPaymentDate() const
{
  CHECK_COND( m_bHasUnknownPayments, 
              ITO33_FLOATINGRATES_NO_UNKNOWN_PAYMENTS );

  return m_lastUnknownPaymentDate;
}

std::vector<Date> 
FloatingRates::GenerateUnknownPaymentDates() const
{
  if ( m_firstUnknownPaymentDate == m_lastUnknownPaymentDate )
    return std::vector<Date>(1, m_firstUnknownPaymentDate);

  std::list<Date> 
    pUnknownPaymentDates(GenerateRegularDates(m_firstUnknownPaymentDate,
                                              m_lastUnknownPaymentDate,
                                              12 / m_paymentfrequency,
                                              IsEOM(m_dcc)));

  // Check if last coupon is irregular
  if (    pUnknownPaymentDates.empty() 
       || pUnknownPaymentDates.back() != m_lastUnknownPaymentDate )
  {
    CHECK_COND( IsValid(m_lastPaymentType),
                ITO33_CASHFLOWSTREAM_INCONSISTENT_DATES );

    if (  m_lastPaymentType == LastPaymentType_Long )
      if ( !pUnknownPaymentDates.empty() )
        pUnknownPaymentDates.pop_back();

    pUnknownPaymentDates.push_back(m_lastUnknownPaymentDate);
  }

  // Add the first unknown payment date
  pUnknownPaymentDates.insert(pUnknownPaymentDates.begin(), 
                              m_firstUnknownPaymentDate);

  return std::vector<Date>( pUnknownPaymentDates.begin(), 
                            pUnknownPaymentDates.end() );
}

shared_ptr<CashFlowStream> 
FloatingRates::ComputeFloatingPayments
( const shared_ptr<YieldCurve>& pYieldcurve ) const
{
  std::vector<Date> pUnknownPaymentDates;
  std::vector<double> pUnknownPaymentRates;
  
  Date previousDate = m_startOfAccruedDate;
  if ( m_pKnownPaymentRatesDates.size() )
    previousDate = 
      m_pKnownPaymentRatesDates[m_pKnownPaymentRatesDates.size() - 1];

  if ( m_bHasUnknownPayments )
  {
    pUnknownPaymentDates = GenerateUnknownPaymentDates();

    pUnknownPaymentRates = 
      ComputeFloatingRates(previousDate, pUnknownPaymentDates,
                          m_dcc, m_iFixingDelay,
                          m_dMargin, m_dMultiplier, m_dCap, m_dFloor,
                          pYieldcurve);
  }
  
  std::vector<Date> pPaymentDates;
  std::vector<double> pPaymentRates;
    
  size_t
    nIdxCoupon,
    nNbCoupons; 

  // Known payment treatment
  
  nNbCoupons = m_pKnownPaymentRatesDates.size();
  for (nIdxCoupon = 0; nIdxCoupon < nNbCoupons; ++nIdxCoupon)
  {    
    pPaymentDates.push_back( m_pKnownPaymentRatesDates[nIdxCoupon] );
    pPaymentRates.push_back( m_pKnownPaymentRates[nIdxCoupon] );
  }

  // Unknown coupons treatment
  
  nNbCoupons = pUnknownPaymentDates.size();
  for (nIdxCoupon = 0; nIdxCoupon < nNbCoupons; ++nIdxCoupon)
  {    
    pPaymentDates.push_back( pUnknownPaymentDates[nIdxCoupon] );
    pPaymentRates.push_back( pUnknownPaymentRates[nIdxCoupon] );
  }

  shared_ptr<CashFlowStream> 
    pcashflowstream;
  
  if ( pPaymentDates.size() )
    pcashflowstream = shared_ptr<CashFlowStream> 
                      ( new CashFlowStreamGeneral
                            ( 
                              m_startOfAccruedDate, pPaymentDates, 
                              pPaymentRates, m_dcc, m_paymentfrequency 
                            )
                      );

  return pcashflowstream;
}  

ito33::XML::Tag FloatingRates::Dump(ito33::XML::Tag& tagParent) const
{
  XML::Tag tagHere(XML_TAG_FLOATINGRATES_ROOT, tagParent);
  
  tagHere.Element(XML_TAG_FLOATINGRATES_MARGIN)
                 ( GetMargin() );

  tagHere.Element(XML_TAG_FLOATINGRATES_MULTIPLIER)
                 ( GetMultiplier() );
  
  tagHere.Element(XML_TAG_FLOATINGRATES_CAP)
                 ( GetCap() );
  
  tagHere.Element(XML_TAG_FLOATINGRATES_FLOOR)
                 ( GetFloor() );
  
  tagHere.Element(XML_TAG_FLOATINGRATES_FIXINGDELAY)
                 ( GetFixingDelay() ); 
  
  //----------------------
  //day count convention
  //----------------------
  tagHere.Element(XML_TAG_FLOATINGRATES_DAYCOUNTCONVENTION)
                 (XML::GetNameOfDayCountConvention(GetDayCountConvention()));

  tagHere.Element(XML_TAG_FLOATINGRATES_STARTOFACCRUEDDATE)
                 ( GetStartOfAccruedDate() );

  if ( m_pKnownPaymentRatesDates.size() )
  {
    XML::Tag tagKnownCoupons(XML_TAG_FLOATINGRATES_KNOWNPAYMENTS, tagHere);
    for ( size_t nIdx = 0; nIdx < m_pKnownPaymentRatesDates.size(); ++nIdx )
    {
      XML::Tag tagElement(XML_TAG_FLOATINGRATES_KNOWNPAYMENT, 
                          tagKnownCoupons);

      tagElement.Element(XML_TAG_FINANCE_DATE)
                        (m_pKnownPaymentRatesDates[nIdx]);
      
      tagElement.Element(XML_TAG_FINANCE_RATE)
                        (m_pKnownPaymentRates[nIdx]);
    }
  }

  if ( m_bHasUnknownPayments )
  {
    tagHere.Element(XML_TAG_FLOATINGRATES_FIRSTPAYMENTDATE)
                  ( GetFirstUnknownPaymentDate() );

    if (m_lastPaymentType != LastPaymentType_Max)
    {
      tagHere.Element(XML_TAG_LASTPAYMENTTYPE)
                     (
                      GetNameFromEnumValue(
                        m_lastPaymentType,
                        SIZEOF(g_LastPaymentType),
                        g_LastPaymentType)
                     );
    }

    tagHere.Element(XML_TAG_FLOATINGRATES_LASTDATE)
                  ( m_lastUnknownPaymentDate );
  }

  //---------------------------
  // payment frequency
  //---------------------------
  tagHere.Element(XML_TAG_PAYMENTFREQUENCY)
                 (
                   GetNameFromEnumValue(
                   GetPaymentFrequency(),
                   SIZEOF(g_frequencys),
                   g_frequencys)
                 );

  return tagHere;
}


} // namespace finance

} // namespace ito33
