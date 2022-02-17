/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/bondterms.cpp
// Purpose:     financial bondcharacteristics class
// Author:      ZHANG Yunzhi
// Created:     2004 may 6
// RCS-ID:      $Id: bondterms.cpp,v 1.48 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/enum_values_names.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/floatingrates.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/utils.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::Error
  ITO33_NULL_FLOATING_RATES;

extern const ito33::finance::BondError
  ITO33_BONDLIKE_WITH_FLOATING_RATES_AND_CASHFLOWS,
  ITO33_BONDLIKE_WITH_FLOATING_RATES_AFTER_MATURITY,
  ITO33_BONDLIKE_TERMS_BADNOMINAL, 
  ITO33_BONDLIKE_TERMS_INCONSISTENT_FREQUENCIES,
  ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY,
  ITO33_BONDLIKE_TERMS_INCONSISTENT_DCCS,
  ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC,
  ITO33_BONDLIKE_TERMS_EITHER_OID_OR_CASHPAYTOZERO,
  ITO33_BONDLIKE_TERMS_MUSTBE_ACCRETING,
  ITO33_BONDLIKE_TERMS_MUSTBE_CASHPAYTOZERO,
  ITO33_BONDLIKE_TERMS_ISSUEPRICE,
  ITO33_BONDLIKE_ISSUEPRICE_TOO_SMALL,
  ITO33_BONDLIKE_MISSING_CASHFLOWS,
  ITO33_BONDLIKE_FIRST_UNKNOWN_PAYMENT_BEFORE_VALUATION,
  ITO33_BONDLIKE_PREVIOUS_REF_DATE_BEFORE_VALUATION_DATE;

namespace ito33
{

namespace finance
{


BondTerms::BondTerms(Date issueDate, 
                     double dIssuePrice,
                     Date maturityDate,
                     double dNominal,                       
                     double dRedemptionPrice,
                     double dRecoveryRate)
  : BondLikeTerms(issueDate, dIssuePrice, maturityDate, dNominal,
                  dRecoveryRate),
    m_dRedemptionPrice(dRedemptionPrice),
    m_bIsAccretingBond(false),
    m_bIsCashPayToZeroBond(false),
    m_cmpFreq(Frequency_Undefined),
    m_yieldDCC(Date::DayCountConvention_Max)
{
}

void BondTerms::Validate() const
{  
  // checks consistency between yieldcompoundingFrequency and paymentFrequency
  // if both have been set
  if ( IsValid(m_cmpFreq) )
  {
    if ( m_pCashDistribution )
    {
      CHECK_COND
        (
          m_cmpFreq == m_pCashDistribution->GetPaymentFrequency(),
          ITO33_BONDLIKE_TERMS_INCONSISTENT_FREQUENCIES
        );
    }
    else if ( m_pFloatingRates )
      CHECK_COND
        (
          m_cmpFreq == m_pFloatingRates->GetPaymentFrequency(),
          ITO33_BONDLIKE_TERMS_INCONSISTENT_FREQUENCIES
        );
  }

  // checks consistency between yield day count convention and coupon
  // day count convention if both have been set
  if ( IsValid(m_yieldDCC) )
  {
    if ( m_pCashDistribution )
    {
      CHECK_COND
        (
          m_yieldDCC == m_pCashDistribution->GetDayCountConvention(),
          ITO33_BONDLIKE_TERMS_INCONSISTENT_DCCS
        );
    }
    else if ( m_pFloatingRates )
      CHECK_COND
        (
          m_yieldDCC == m_pFloatingRates->GetDayCountConvention(),
          ITO33_BONDLIKE_TERMS_INCONSISTENT_DCCS
        );
  }

  if (    IsAccretingBond()  
       && !IsValid( GetYieldCompoundingFrequencyEvenIfUndefined(*this) ) )
    throw
       EXCEPTION(ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY);

  if (    IsAccretingBond()  
       && !IsValid( GetYieldDayCountConventionEvenIfUndefined(*this) ) )
    throw
       EXCEPTION(ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC);

  // Issue price should be at least one percent smaller than redemption price
  // when the bond accretes
  if ( IsAccretingBond() || IsCashPayToZeroBond() )
  {
    double dYield = - 1;

    if ( IsAccretingBond() )
      dYield = GetYieldToMaturityOfAccretingBond();
    else
      dYield = GetAccretionRateOfCashPayToZero();

    CHECK_COND( dYield <= 0 || GetIssuePrice() <= 0.99 * GetRedemptionPrice(),
                ITO33_BONDLIKE_TERMS_ISSUEPRICE );
  }
  else // non accreting bond, issue price shouldn't be too smaller than 1
  {
    CHECK_COND( GetIssuePrice() >= 0.94, ITO33_BONDLIKE_ISSUEPRICE_TOO_SMALL);  
  }

  // CashPayToZeroBond requires the existence of coupons
  if ( IsCashPayToZeroBond() )
  {
    shared_ptr<CashFlowStream> pcashFlowStream = GetCashDistribution();

    if ( GetFloatingRates() )
    {
      // TODO: This is really a ugly work around because of limited interface of
      // FloatingRates and needs to be fixed.
      shared_ptr<YieldCurveFlat> pYC( new YieldCurveFlat(0.02) );
      
      pcashFlowStream = GetFloatingRates()->ComputeFloatingPayments(pYC);
    }

    CHECK_COND( !pcashFlowStream->GetAll().empty(), 
                ITO33_BONDLIKE_MISSING_CASHFLOWS );
  }
}

void BondTerms::ValidateWith(const SessionData& sessionData) const
{
  if ( m_pFloatingRates && m_pFloatingRates->HasUnknownPayments() )
  {
    CHECK_COND
    ( 
      m_pFloatingRates->GetFirstUnknownPaymentDate() > 
      sessionData.GetValuationDate(), 
      ITO33_BONDLIKE_FIRST_UNKNOWN_PAYMENT_BEFORE_VALUATION
    );

    /* 
      Checks if (PreviousResetDate - FixingDelay > Valuation Date) otherwise,
      first unknown coupon should be known (if other data are correct and 
      complete)
      Previous reset date is the last known payment date (if any) or the start 
      of accrued date.
    */
    Date previousResetDate;

    if( m_pFloatingRates->GetKnownPaymentDates().size() > 0 )
      previousResetDate = m_pFloatingRates->GetKnownPaymentDates().back();
    else
      previousResetDate = m_pFloatingRates->GetStartOfAccruedDate();

    CHECK_COND
    ( 
      previousResetDate.AddDays( - m_pFloatingRates->GetFixingDelay() ) > 
      sessionData.GetValuationDate(),
      ITO33_BONDLIKE_PREVIOUS_REF_DATE_BEFORE_VALUATION_DATE                
    );
  }
}

void BondTerms::SetCashDistribution
     (const shared_ptr<finance::CashFlowStream>& pCashDistribution)
{
  CHECK_COND
        (
          !m_pFloatingRates,
          ITO33_BONDLIKE_WITH_FLOATING_RATES_AND_CASHFLOWS
        );

  BondLikeTerms::SetCashDistribution(pCashDistribution);
}

void BondTerms::SetFloatingRates
     ( const shared_ptr<FloatingRates>& pFloatingRates )
{
  CHECK_COND
        (
          !m_pCashDistribution,
          ITO33_BONDLIKE_WITH_FLOATING_RATES_AND_CASHFLOWS
        );

  m_pFloatingRates = CHECK_PTR( pFloatingRates, ITO33_NULL_FLOATING_RATES );
  
  Date lastPaymentDateOfFloating;

  if ( m_pFloatingRates->HasUnknownPayments() )
    lastPaymentDateOfFloating = m_pFloatingRates->GetLastUnknownPaymentDate();
  else if ( m_pFloatingRates->GetKnownPaymentDates().size() > 0 )
    lastPaymentDateOfFloating = m_pFloatingRates->GetKnownPaymentDates().back();
  
  if ( lastPaymentDateOfFloating.IsValid() )
    CHECK_COND
      ( m_maturityDate >= lastPaymentDateOfFloating, 
        ITO33_BONDLIKE_WITH_FLOATING_RATES_AFTER_MATURITY );                      
}

void BondTerms::SetAccretingBond(double dYield)
{
  // some jp bond has negative yield, so nothing to check with yield value

  CHECK_COND
    (
      !m_bIsCashPayToZeroBond,
      ITO33_BONDLIKE_TERMS_EITHER_OID_OR_CASHPAYTOZERO
    );

  m_bIsAccretingBond = true;
  m_dYield = dYield;
}


double BondTerms::GetYieldToMaturityOfAccretingBond() const
{ 
  CHECK_COND( IsAccretingBond(), ITO33_BONDLIKE_TERMS_MUSTBE_ACCRETING);
  return m_dYield; 
}

void BondTerms::SetCashPayToZero(double dYield)
{
  // as we support negative accretingBondYield,
  // we also support negative yield here
  
  CHECK_COND(
    !m_bIsAccretingBond,
    ITO33_BONDLIKE_TERMS_EITHER_OID_OR_CASHPAYTOZERO);

  m_bIsCashPayToZeroBond = true;
  m_dAccretionRateOfCashPayToZero = dYield;
}

double BondTerms::GetAccretionRateOfCashPayToZero() const
{
  CHECK_COND( IsCashPayToZeroBond(),
              ITO33_BONDLIKE_TERMS_MUSTBE_CASHPAYTOZERO);
  return m_dAccretionRateOfCashPayToZero;
}

void 
BondTerms::SetYieldCompoundingFrequency(Frequency compoundingFrequency)
{ 
  finance::Validate(compoundingFrequency);

  m_cmpFreq = compoundingFrequency;
}

void 
BondTerms::SetYieldDayCountConvention(Date::DayCountConvention dcc)
{ 
  ::ito33::Validate(dcc);

  m_yieldDCC = dcc;
}

Frequency BondTerms::GetYieldCompoundingFrequency() const
{
  CHECK_COND( IsValid(m_cmpFreq),
              ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY );

  return m_cmpFreq;
}

Date::DayCountConvention BondTerms::GetYieldDayCountConvention() const
{
  CHECK_COND( IsValid(m_yieldDCC), ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC );

  return m_yieldDCC;
}

double BondTerms::GetRedemptionPrice() const
{
  return m_dRedemptionPrice;
}

XML::Tag BondTerms::Dump(ito33::XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDTERMS_ROOT, tagParent);

  BondLikeTerms::DumpMe(tagMe);

  if ( m_bIsAccretingBond )
  {
    tagMe.Element(XML_TAG_BONDTERMS_OIDYIELD)(m_dYield);
  }

  if ( m_bIsCashPayToZeroBond )
  {
    tagMe.Element(XML_TAG_BONDTERMS_CASHPAYTOZEROACCRETIONRATE)
                 (m_dAccretionRateOfCashPayToZero);
  }

  tagMe.Element(XML_TAG_BONDTERMS_REDEMPTIONPRICE)(m_dRedemptionPrice);

  // compounding  frequency
  if ( IsValid(m_cmpFreq) )
    tagMe.Element(XML_TAG_BONDTERMS_COMPOUNDING_FREQUENCY)
                     (
                       GetNameFromEnumValue(
                         m_cmpFreq,
                         SIZEOF(g_frequencys),
                         g_frequencys)
                     );

  // yield day count convention
  if ( IsValid(m_yieldDCC) )
    tagMe.Element(XML_TAG_BONDTERMS_YIELD_DCC)
                 (
                   GetNameFromEnumValue(
                     m_yieldDCC,
                     SIZEOF(g_dayCountConventions),
                     g_dayCountConventions)
                 );

  if ( m_pFloatingRates )
    m_pFloatingRates->Dump(tagMe);
  
  return tagMe;
}

} // namespace finance

} // namespace ito33
