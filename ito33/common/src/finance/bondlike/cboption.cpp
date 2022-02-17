/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/cboption.cpp
// Purpose:     financial convertible bond option class
// Author:      Nabil
// Created:     2005/06/07
// RCS-ID:      $Id: cboption.cpp,v 1.18 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/rateutils.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/cboption.h"
#include "ito33/finance/bondlike/putschedule.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/bondlike/cboption.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"

extern const ito33::Error 
  ITO33_INVALID_FREQUENCY;

extern const ito33::finance::BondError
  ITO33_CBOPTION_WITH_FLOATING_CB,
  ITO33_CBOPTION_WITH_INVALID_CB,
  ITO33_CBOPTION_WITH_BAD_MATURITY,
  ITO33_CBOPTION_INVALID_FLOATINGRATES,
  ITO33_CBOPTION_INCOHERENT_ASWNOTIONAL,
  ITO33_CBOPTION_AND_CB_NOT_SAME_CURRENCY,
  ITO33_CBOPTION_AND_CB_NOT_SAME_SESSIONDATA,
  ITO33_CBOPTION_WITH_FLOATING_RATES_AFTER_MATURITY,
  ITO33_CBOPTION_FIRST_UNKNOWN_PAYMENT_BEFORE_VALUATION,
  ITO33_CBOPTION_PREVIOUS_REF_DATE_BEFORE_VALUATION_DATE;

namespace ito33
{

namespace finance
{

CBOption::CBOption(const shared_ptr<ConvertibleBond>& pConvertibleBond,
                   const shared_ptr<FloatingRates>& pFloatingRates,
                   Date maturityDate)
         : m_maturityDate(maturityDate),
           m_ASWNotionalType(ASWNotionalIs_IssuePrice)
{
  m_pConvertibleBond = CHECK_PTR
                       (
                        pConvertibleBond,
                        ITO33_CBOPTION_WITH_INVALID_CB
                       );

  m_pFloatingRates = 
    CHECK_PTR(pFloatingRates, ITO33_CBOPTION_INVALID_FLOATINGRATES);
}

void CBOption::SetSessionData(const shared_ptr<SessionData>& pSessionData)
{
  if ( m_pConvertibleBond->GetSessionData() )
    CHECK_COND( pSessionData == m_pConvertibleBond->GetSessionData(),
                ITO33_CBOPTION_AND_CB_NOT_SAME_SESSIONDATA );

  Derivative::SetSessionData(pSessionData);

  m_pConvertibleBond->SetSessionData(pSessionData);
}

void CBOption::Validate() const
{
  Derivative::Validate();
  
  // Validate the underlying convertible
  m_pConvertibleBond->Validate();

  CHECK_COND( !m_pConvertibleBond->GetBondTerms()->GetFloatingRates(), 
              ITO33_CBOPTION_WITH_FLOATING_CB );
  
  // The currency of the cboption and the underlying cb must be the same if
  // both are set
  if ( m_pConvertibleBond->GetNumeraire() && GetNumeraire() )
  {
    CHECK_COND
    ( 
      *( m_pConvertibleBond->GetNumeraire() ) == *GetNumeraire(), 
      ITO33_CBOPTION_AND_CB_NOT_SAME_CURRENCY
    );
  }
  else
  {
    CHECK_COND( !m_pConvertibleBond->GetNumeraire() && !GetNumeraire(),
                ITO33_CBOPTION_AND_CB_NOT_SAME_CURRENCY );
  }

  Date lastPaymentDateOfFloating;

  if ( m_pFloatingRates->HasUnknownPayments() )
    lastPaymentDateOfFloating = m_pFloatingRates->GetLastUnknownPaymentDate();
  else if ( m_pFloatingRates->GetKnownPaymentDates().size() > 0 )
    lastPaymentDateOfFloating = 
      m_pFloatingRates->GetKnownPaymentDates().back();
  
  if ( lastPaymentDateOfFloating.IsValid() )
    CHECK_COND
      ( m_maturityDate >= lastPaymentDateOfFloating, 
        ITO33_CBOPTION_WITH_FLOATING_RATES_AFTER_MATURITY);

  shared_ptr<PutSchedule> 
    pPutSchedule( m_pConvertibleBond->GetPutSchedule() );

  // If maturity is not a put date
  if ( !pPutSchedule || !pPutSchedule->HasPutAt(m_maturityDate) )
  {
    CHECK_COND( m_maturityDate == m_pConvertibleBond->GetMaturityDate(),
                ITO33_CBOPTION_WITH_BAD_MATURITY );

    CHECK_COND( m_ASWNotionalType != ASWNotionalIs_PutPrice,
                ITO33_CBOPTION_INCOHERENT_ASWNOTIONAL );
  }
}

void CBOption::ValidateWith(const SessionData& sessionData) const
{
  Derivative::ValidateWith(sessionData);

  // cross validation with the convertible and the underlying
  m_pConvertibleBond->ValidateWith(sessionData);

  if ( m_pFloatingRates->HasUnknownPayments() )
    CHECK_COND
    ( 
      m_pFloatingRates->GetFirstUnknownPaymentDate() > 
      sessionData.GetValuationDate(), 
      ITO33_CBOPTION_FIRST_UNKNOWN_PAYMENT_BEFORE_VALUATION
    );

  /* 
    Checks if (PreviousResetDate - FixingDelay > Valuation Date) otherwise,
    first unknown payment should be known (if other data are correct and 
    complete)
    Previous reset date is the last known payment date (if any) or the start 
    of accrued date.
  */
  Date previousResetDate;

  if ( !m_pFloatingRates->GetKnownPaymentDates().empty() )
    previousResetDate = m_pFloatingRates->GetKnownPaymentDates().back();
  else
    previousResetDate = m_pFloatingRates->GetStartOfAccruedDate();

  CHECK_COND
  ( 
    previousResetDate.AddDays( - m_pFloatingRates->GetFixingDelay() ) > 
    sessionData.GetValuationDate(),
    ITO33_CBOPTION_PREVIOUS_REF_DATE_BEFORE_VALUATION_DATE                
  );
}

void CBOption::ValidateAll() const
{
  Derivative::ValidateAll();

  if ( m_pConvertibleBond->GetSessionData() )
    CHECK_COND( GetSessionData() == m_pConvertibleBond->GetSessionData(),
                ITO33_CBOPTION_AND_CB_NOT_SAME_SESSIONDATA );
}

void CBOption::SetASWNotionalIs(ASWNotionalIs ASWNotionalType)
{
  m_ASWNotionalType = ASWNotionalType;
}

const shared_ptr<finance::YieldCurve>& CBOption::GetDerivativeCurve() const 
{ 
  // Check that the session data is set
  CheckSessionData();

  shared_ptr<SessionData> pSessionData = GetSessionData();

  // find the currency of the cboption.
  shared_ptr<Numeraire> pNumeraire = m_pNumeraire;
  if ( !pNumeraire )
    pNumeraire = pSessionData->GetNumeraire();

  return pSessionData->GetRateData()->GetYieldCurve(pNumeraire);
}  

double CBOption::GetCbNominal() const
{
  return m_pConvertibleBond->GetBondLikeTerms()->GetNominal();    
}

double CBOption::GetASWRedemptionRate() const
{
  // If not at maturity, then must be at a put date
  if ( m_maturityDate != m_pConvertibleBond->GetMaturityDate() )
  {
    double dASWRedemptionRate;

    double dCouponAmount;
    dASWRedemptionRate = m_pConvertibleBond->ComputePutPrice
                              ( m_maturityDate, dCouponAmount );
    dASWRedemptionRate -= dCouponAmount;
    dASWRedemptionRate /= GetCbNominal();
    
    return dASWRedemptionRate;
  }

  return m_pConvertibleBond->GetBondTerms()->GetRedemptionPrice();
}

double CBOption::GetASWNotional() const
{
  double dASWNotional = 0.;
  
  switch( m_ASWNotionalType )
  {
    case ASWNotionalIs_IssuePrice:
      dASWNotional = 
        m_pConvertibleBond->GetBondLikeTerms()->GetIssuePrice() *
        m_pConvertibleBond->GetBondLikeTerms()->GetNominal();
      break;
    case ASWNotionalIs_PutPrice:
      dASWNotional = GetASWRedemptionRate() *
        m_pConvertibleBond->GetBondLikeTerms()->GetNominal();
      break;
  }
  
  return dASWNotional;
}

double CBOption::GetBalloonCoupon() const
{
  return ( GetASWRedemptionRate() 
           * m_pConvertibleBond->GetBondLikeTerms()->GetNominal()
           - GetASWNotional() );
}

shared_ptr<CashFlowStream> CBOption::ComputeASWFloatingPayments() const
{
  shared_ptr<YieldCurve> pYieldcurve = GetDerivativeCurve();

  return m_pFloatingRates->ComputeFloatingPayments(pYieldcurve);
}  
  
shared_ptr<CashFlowStream> CBOption::ComputeASWFixedPayments() const
{
  shared_ptr<CashFlowStream>
    pCBCashFlows = m_pConvertibleBond->GetBondTerms()->GetCashDistribution();

  if ( !pCBCashFlows )
    return pCBCashFlows;

  std::vector<Date> pPaymentDates;
  std::vector<double> pPaymentRates;
  
  std::vector<CashFlow> 
    pCashFlows = pCBCashFlows->GetAll();

  size_t
    nIdxCoupon,
    nNbCoupons = pCashFlows.size();

  for ( nIdxCoupon = 0; nIdxCoupon < nNbCoupons; ++nIdxCoupon )
  {
    if ( pCashFlows[nIdxCoupon].GetDate() <= m_maturityDate )
    {
      pPaymentDates.push_back( pCashFlows[nIdxCoupon].GetDate() );
      pPaymentRates.push_back( pCashFlows[nIdxCoupon].GetAmount() );
    }
    else
      break;
  }
  
  shared_ptr<CashFlowStream> 
    pcashflowstream(
                     new CashFlowStreamGeneral
                         ( 
                           pCBCashFlows->GetContractingDate(), pPaymentDates, 
                           pPaymentRates, 
                           pCBCashFlows->GetDayCountConvention(),
                           pCBCashFlows->GetPaymentFrequency()
                         )
                   );

  return pcashflowstream;
}

double ComputeNPV(const CashFlowStream& flows,
                  const Date& valuation,
                  const shared_ptr<YieldCurve>& pYC)
{
  double dValue = 0;
  double dTime = GetDoubleFrom(valuation);

  CashFlowStream::const_iterator iter;
  for (iter = flows.begin(); iter != flows.end(); ++iter)
  {
    if ( iter->first > valuation )
    {
      dValue += iter->second 
              * pYC->GetForwardDiscountFactor(dTime,
                                              GetDoubleFrom(iter->first));
    }
  }

  return dValue;
}

double CBOption::ComputeStrike() const
{
  // Check that the session data is set
  CheckSessionData();

  shared_ptr<SessionData> pSessionData = GetSessionData();

  // Gets the right yield curve.
  shared_ptr<YieldCurve> pYC = GetDerivativeCurve();

  shared_ptr<CashFlowStream>
    pFloating = ComputeASWFloatingPayments();

  double
    dNPVFixedLeg,
    dNPVFloatingLeg,
    dFloatingAccrued;

  double dFactorToMaturity = 
          pYC->GetForwardDiscountFactor( 
            GetDoubleFrom(pSessionData->GetValuationDate()),
            GetDoubleFrom(GetMaturityDate()) );

  shared_ptr<CashFlowStream>
    pFixed = ComputeASWFixedPayments();

  // Fixed leg
  if ( pFixed )
    dNPVFixedLeg = ComputeNPV( *pFixed, pSessionData->GetValuationDate(), pYC);
  else
    dNPVFixedLeg = 0.;

  // Floating leg
  dNPVFloatingLeg = ComputeNPV
      ( *pFloating, pSessionData->GetValuationDate(), pYC );
  
  dFloatingAccrued = pFloating->GetAccrued(pSessionData->GetValuationDate());

  // The strike of the cb option
  return 
    GetASWNotional()
    + (GetConvertibleBond()->GetBondTerms()->GetNominal() * dNPVFixedLeg)
    + GetASWNotional() * (dFloatingAccrued - dNPVFloatingLeg)
    + dFactorToMaturity * GetBalloonCoupon();
}

void CBOption::PerturbFXRate(double dFXRateShift) const
{
  m_pConvertibleBond->PerturbFXRate(dFXRateShift);
}

void CBOption::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnCBOption(*this);
}

XML::Tag CBOption::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagCBOption(XML_TAG_CBOPTION_ROOT, tagParent);

  DumpMe(tagCBOption);

  tagCBOption.Element(XML_TAG_MATURITYDATE)( GetMaturityDate() );
  
  tagCBOption.Element(XML_TAG_ASWNOTIONALTYPE)
                     (
                       GetNameFromEnumValue(
                         m_ASWNotionalType,
                         SIZEOF(g_ASWNotionalIs),
                         g_ASWNotionalIs)
                     );

  // dump data of the underlying cb
  tagCBOption.Element(*m_pConvertibleBond);
  
  // dump data of the floating rates
  m_pFloatingRates->Dump(tagCBOption);

  DumpMarketPrice(tagCBOption);

  return tagCBOption;
}

} // namespace finance

} // namespace ito33
