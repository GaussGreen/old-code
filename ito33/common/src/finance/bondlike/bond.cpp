/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/bond.cpp
// Purpose:     financial bond class
// Author:      ZHANG Yunzhi
// Created:     2004 may 6
// RCS-ID:      $Id: bond.cpp,v 1.63 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/utils.h"

#include "ito33/numeric/bisecnewton.h"
#include "ito33/numeric/exception.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/bondlike/bond.h"

extern const ito33::finance::BondError
  ITO33_BONDLIKE_NULL_BONDTERMS,
  ITO33_BONDLIKE_BOND_SOFTCALL,
  ITO33_BONDLIKE_COMPUTE_YTP_INEXIST,
  ITO33_BONDLIKE_COMPUTE_YTP_NOPUT,
  ITO33_BONDLIKE_COMPUTE_YTP_NOPUT_AFTER_VALUATIONDATE,
  ITO33_BONDLIKE_COMPUTE_YTM_INEXIST,
  ITO33_BONDLIKE_NULL_CALLSCHEDULE,
  ITO33_BONDLIKE_NULL_PUTSCHEDULE,
  ITO33_BONDLIKE_NO_CALL,
  ITO33_BONDLIKE_NO_PUT;

namespace ito33
{

namespace finance
{

/**
   Class for computing the yield-to-maturity or yield-to-put
 */
class YieldCalculator
{

public:

  /**
     Ctor using Bond

     @param dPrice The price used to compute the yield.
     @param maturityDate The maturity date used to compute the yield.
     @param pBond The bond.
     @param cmpFrequency The yield compound frequency
   */
  YieldCalculator(double dPrice, Date maturityDate, double dAmountAtMaturity,
                  const Bond& Bond, Frequency compoundingFrequency, 
                  Date::DayCountConvention dcc) 
                : m_dPrice(dPrice), m_maturityDate(maturityDate),
                  m_dAmountAtMaturity(dAmountAtMaturity), m_Bond(Bond), 
                  m_cmpFrequency(compoundingFrequency), m_dcc(dcc)
  { }

  // Default dtor is ok

  /** 
      Used by the non-linear root finding algorithm.

      The yield is the root of the function 
      \f$ 
        F(x)= P_0 (1 + \frac{x}{f})^{f(T-t_0)} - AmountAtMaturity - 
        \sum_{T > t_i > t_0} C_i (1 + \frac{x}{f})^{f(T-t_i)} 
      \f$

      where:
       1) \f$ P_0 \f$ is the price used to compute the yield
       2) \f$ T \f$ is the maturity date used to compute the yield:
          it is the maturity of the instrument is we compute the YTM
          and the put date if it is the YTP which we compute.
       3) \f$ AmountAtMaturity = redemption + couponAtMaturity \f$ for the YTM
          \f$ AmountAtMaturity = PutValue \f$ for the YTP

      @param dX The current guess of the root (see details below)
      @param dF (output) The value of the non-linear function at dX
      @param dGradF (output) The gradient of the non-linear function at dX
   */
  void operator () (double dX, double& dF, double& dGradF) const
  {     
    shared_ptr<CashFlowStream> 
      pCashDistribution = ComputeCouponRates( *m_Bond.GetBondTerms(), 
                                              *m_Bond.GetSessionData(),
                                              m_Bond.GetNumeraire() );
 
    double dNominal = m_Bond.GetBondTerms()->GetNominal();

    Date valuationDate = m_Bond.GetSessionData()->GetValuationDate();

    double
      dTimeToMaturity = Date::YearsDiff(valuationDate, m_maturityDate, m_dcc);

    double dY = 1. + dX / m_cmpFrequency;

    dF = m_dPrice
       - m_dAmountAtMaturity * pow(dY, - m_cmpFrequency * dTimeToMaturity) ;

    dGradF = dTimeToMaturity * m_dAmountAtMaturity 
           * pow(dY, - m_cmpFrequency * dTimeToMaturity - 1 );

    if ( pCashDistribution )
    {
      CashFlowStream::Elements::const_iterator iter;
      
      for (iter = pCashDistribution->begin(); 
           iter != pCashDistribution->end() && iter->first <= valuationDate;
           ++iter);

      for (; iter != pCashDistribution->end() && iter->first < m_maturityDate; 
             ++iter)
      {
        dTimeToMaturity = Date::YearsDiff( valuationDate, iter->first, m_dcc );

        dF -= iter->second * dNominal 
            * pow(dY, -m_cmpFrequency * dTimeToMaturity);
        dGradF += dTimeToMaturity * iter->second * dNominal 
                * pow(dY, -m_cmpFrequency * dTimeToMaturity - 1);
      }
    }

  }

protected:

  // Price used to compute the yield to maturity.
  double m_dPrice;

  // Amount at maturity used to compute the yield
  double m_dAmountAtMaturity;

  // The maturity date used to compute the yield.
  Date m_maturityDate;

  const Bond& m_Bond;

  Frequency m_cmpFrequency;

  Date::DayCountConvention m_dcc;

  NO_COPY_CLASS(YieldCalculator);

}; // class YieldToMaturity


Bond::Bond(const shared_ptr<BondTerms>& pBondTerms)
         : Derivative(), m_pBondTerms(pBondTerms)
{ 
  CHECK_PTR(pBondTerms, ITO33_BONDLIKE_NULL_BONDTERMS);
}

void Bond::SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule)
{
  m_pCallSchedule = CHECK_PTR(pCallSchedule, ITO33_BONDLIKE_NULL_CALLSCHEDULE);

  // Check that there is no softcall for a bond, otherwise the conversion price
  // is required for the trigger but not well defined since there is no 
  // conversion for a normal bond.
  CHECK_COND( !( m_pCallSchedule->HasSoftCall() ),
              ITO33_BONDLIKE_BOND_SOFTCALL );
}

void Bond::SetPutSchedule(const shared_ptr<PutSchedule>& pPutSchedule)
{
  m_pPutSchedule = CHECK_PTR(pPutSchedule, ITO33_BONDLIKE_NULL_PUTSCHEDULE);
}

double Bond::GetAccruedInterestValue() const
{
  ValidateAll();

  shared_ptr<CashFlowStream> 
    pCoupons = ComputeCouponRates(*m_pBondTerms, *GetSessionData(), m_pNumeraire);

  if ( pCoupons )
    return  pCoupons->GetAccrued( GetSessionData()->GetValuationDate() )
          * m_pBondTerms->GetNominal();
  else
    return 0;
}

Date Bond::GetMaturityDate() const
{
  return m_pBondTerms->GetMaturityDate();
}

double 
Bond::ComputeYieldToMaturity
(double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const
{
  ValidateAll();

  finance::Validate(cmpFrequency);

  ito33::Validate(dcc);

  Date maturityDate = m_pBondTerms->GetMaturityDate();

  double dAmountAtMaturity = m_pBondTerms->GetNominal() 
                           * m_pBondTerms->GetRedemptionPrice();

  shared_ptr<CashFlowStream> 
    pCashDistribution = ComputeCouponRates( *m_pBondTerms, 
                                            *GetSessionData(), 
                                            m_pNumeraire );
  
  if ( pCashDistribution )
  {
    CashFlowStream::Elements::const_iterator iter;
      
    for(iter = pCashDistribution->begin(); 
        iter != pCashDistribution->end() && iter->first < maturityDate;
        ++iter)
    {
    }

    if ( iter != pCashDistribution->end() && iter->first == maturityDate )
      dAmountAtMaturity += iter->second * m_pBondTerms->GetNominal();
  }

  YieldCalculator 
    YTM(dPrice, maturityDate, dAmountAtMaturity, *this, cmpFrequency, dcc);

  numeric::BisecNewton solver(1.e-6);

  double dResult = -10.;
  
  try
  {
    double
      dYTMmin, dYTMmax;
    
    dYTMmin = -0.9999 * cmpFrequency;
    dYTMmax = 2;
    
    solver.SetInitialGuess( pow(dAmountAtMaturity / dPrice,
                1 / (GetDoubleFrom(maturityDate)
                  - GetDoubleFrom(GetSessionData()->GetValuationDate()))) -1  );
    dResult = solver(YTM, dYTMmin, dYTMmax);
  }
  catch (const ito33::numeric::Exception&)
  {
    throw EXCEPTION(ITO33_BONDLIKE_COMPUTE_YTM_INEXIST);
  }

  return dResult; 
}

double Bond::ComputePutPrice(Date putDate) const
{
  double dCouponAmount;
  return ComputePutPrice(putDate, dCouponAmount);
}

double Bond::ComputePutPrice(Date putDate, double& dCouponAmount) const
{ 
  ValidateAll();

  CHECK_COND(m_pPutSchedule, ITO33_BONDLIKE_NO_PUT);

  dCouponAmount = 0.;

  PutSchedule::Elements::const_iterator iterPut;
  PutSchedule::Elements puts = m_pPutSchedule->GetAll();
      
  for (iterPut = puts.begin(); 
       iterPut != puts.end() && iterPut->first != putDate;
       ++iterPut);
  
  CHECK_COND( iterPut != puts.end(), ITO33_BONDLIKE_NO_PUT);
   
  double dNominal, dClaim, dPutValue;

  dNominal = m_pBondTerms->GetNominal();

  // Gets the coupon (if there are one) at the given date
  shared_ptr<CashFlowStream> 
    pCashFlowStream = ComputeCouponRates( *m_pBondTerms, 
                                          *GetSessionData(),
                                          m_pNumeraire );

  if ( pCashFlowStream )
  {
    CashFlowStream::Elements::const_iterator iter;
      
    for (iter = pCashFlowStream->begin(); 
         iter != pCashFlowStream->end() && iter->first != putDate;
         ++iter);
    
    if ( iter != pCashFlowStream->end() )
      dCouponAmount = iter->second * dNominal;
  }
  
  // Compute the put value at the given date
  if ( iterPut->second.bYield )
  {
    dClaim = GetClaimFromYield( putDate, iterPut->second.dYield, *m_pBondTerms, 
                                *GetSessionData(), m_pNumeraire );
    dPutValue = dClaim + dCouponAmount;
  }
  else
  {
    double 
      dAccruedAmount = 0.,
      dBase;

    if ( pCashFlowStream )
      dAccruedAmount = pCashFlowStream->GetAccrued(putDate) * dNominal;
    
    // strike is applied to principal = claim minus accrued
    dClaim = GetClaim( putDate, *m_pBondTerms, *GetSessionData(), m_pNumeraire );
    dBase = dClaim - dAccruedAmount;
    dPutValue = iterPut->second.dStrike * dBase;

    // accrued interest and forfeit coupon only apply for strikes
    if ( m_pPutSchedule->GetKeepAccrued() )
      dPutValue += dAccruedAmount;

    if ( !m_pPutSchedule->GetForfeitCoupon() )
      dPutValue += dCouponAmount;
  }

  return dPutValue;
}

double Bond::ComputeCallPrice(Date callDate) const
{ 
  ValidateAll();

  CHECK_COND(m_pCallSchedule, ITO33_BONDLIKE_NO_CALL);

  return ComputeBondLikeCallPrice( callDate, *m_pBondTerms, *GetSessionData(), 
                                   *m_pCallSchedule, m_pNumeraire );
}

double 
Bond::ComputeYieldToPut
(double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const
{  
  ValidateAll();

  finance::Validate(cmpFrequency);

  ito33::Validate(dcc);

  CHECK_COND
  (
    m_pPutSchedule,
    ITO33_BONDLIKE_COMPUTE_YTP_NOPUT
  );
  
  Date
    valuationDate( GetSessionData()->GetValuationDate() ),
    dateNextPut;

  PutSchedule::Elements::const_iterator iterPut;
  PutSchedule::Elements puts = m_pPutSchedule->GetAll();
      
  for (iterPut = puts.begin(); 
       iterPut != puts.end() && iterPut->first <= valuationDate;
       ++iterPut);
  
  CHECK_COND
    (
      iterPut != puts.end(),
      ITO33_BONDLIKE_COMPUTE_YTP_NOPUT_AFTER_VALUATIONDATE
    );
  
  dateNextPut = iterPut->first;  

  double dAmountAtMaturity = ComputePutPrice(dateNextPut);

  YieldCalculator 
    YTP(dPrice, dateNextPut, dAmountAtMaturity, *this, cmpFrequency, dcc);

  numeric::BisecNewton solver(1.e-6);

  double dResult = -10.;
  try
  {
    solver.SetInitialGuess( pow(dAmountAtMaturity / dPrice,
                1 / (GetDoubleFrom(dateNextPut) 
                  - GetDoubleFrom(GetSessionData()->GetValuationDate()))) -1  );
    dResult = solver(YTP, -.99999 * cmpFrequency, 2);
  }  
  catch (const ito33::numeric::Exception&)
  {
    throw EXCEPTION(ITO33_BONDLIKE_COMPUTE_YTP_INEXIST);
  }

  return dResult;
}

void Bond::Validate() const
{
  Derivative::Validate();

  m_pBondTerms->Validate();

  if ( m_pCallSchedule )
    m_pCallSchedule->Validate();

  if ( m_pPutSchedule )
    m_pPutSchedule->Validate();

  CheckYieldCompoundingFrequency
      (m_pPutSchedule, m_pCallSchedule, *m_pBondTerms);

  CheckYieldDayCountConvention
      (m_pPutSchedule, m_pCallSchedule, *m_pBondTerms);

  // Soft call is not allowed
  if ( m_pCallSchedule )
    CHECK_COND( !( m_pCallSchedule->HasSoftCall() ),
                ITO33_BONDLIKE_BOND_SOFTCALL );
}

void Bond::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnBond(*this);
}

void Bond::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnBond(*this);
}

XML::Tag Bond::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagBond(XML_TAG_BOND_ROOT, tagParent);
  
  DumpMe(tagBond);

  tagBond.Element(*m_pBondTerms);

  if ( m_pCallSchedule && !( m_pCallSchedule->GetAll().empty() ) )
    tagBond.Element(*m_pCallSchedule);

  if ( m_pPutSchedule && !( m_pPutSchedule->GetAll().empty() ) )
    tagBond.Element(*m_pPutSchedule);

  DumpMarketPrice(tagBond);

  return tagBond;
}


} // namespace finance

} // namespace ito33

