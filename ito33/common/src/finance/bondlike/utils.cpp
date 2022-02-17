/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/utils.cpp
// Purpose:     useful functions for using bondlike classes
// Author:      ZHANG Yunzhi
// Created:     2005/02/24
// RCS-ID:      $Id: utils.cpp,v 1.25 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/floatingrates.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/numeraire.h"

#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/utils.h"

extern const ito33::finance::BondError
  ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY,
  ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC,
  ITO33_BONDLIKE_NO_CALL;

namespace ito33
{

namespace finance
{

shared_ptr<CashFlowStream> 
ComputeCouponRates(const BondTerms& terms, const SessionData& sessionData,
                   const shared_ptr<Numeraire>& pNumeraire)
{
  shared_ptr<CashFlowStream> 
    pCashDistribution = terms.GetCashDistribution();
    
  if ( !pCashDistribution && terms.GetFloatingRates() )
  {
    shared_ptr<finance::YieldCurve> pDerivativeCurve;

    if ( pNumeraire && 
         ( *pNumeraire != *(sessionData.GetEquity()->GetNumeraire()) ) ) 
    {
      // if IsCrossCurrency
      pDerivativeCurve = sessionData.GetRateData()->GetYieldCurve(pNumeraire);
    }
    else
      pDerivativeCurve = sessionData.GetYieldCurve();

    pCashDistribution = terms.GetFloatingRates()
                      ->ComputeFloatingPayments( pDerivativeCurve );
  }

  return pCashDistribution;
}

double GetClaimFromYield(Date date, double dYield,
                         const BondTerms& bondTerms,
                         const SessionData& sessionData,
                         const shared_ptr<Numeraire>& pNumeraire)
{  
  Date issueDate = bondTerms.GetIssueDate();
  
  double dNominal = bondTerms.GetNominal();

  Frequency 
    cmpFrequency = GetYieldCompoundingFrequencyEvenIfUndefined(bondTerms);

  Date::DayCountConvention 
    dcc = GetYieldDayCountConventionEvenIfUndefined(bondTerms);

  double dDiscount = 1. + dYield / cmpFrequency;
  double dClaim = bondTerms.GetIssuePrice() 
                * pow( dDiscount, 
                       cmpFrequency * Date::YearsDiff(issueDate, date, dcc) 
                     );
  
  shared_ptr<CashFlowStream> 
    pCashFlowStream = ComputeCouponRates(bondTerms, sessionData, pNumeraire);

  if ( pCashFlowStream )
  {
    CashFlowStream::Elements::const_iterator iter;
      
    for (iter = pCashFlowStream->begin(); 
         iter != pCashFlowStream->end() && iter->first <= date;
         ++iter)
    {
      dClaim -= iter->second
              * pow( dDiscount, 
                     cmpFrequency * Date::YearsDiff(iter->first, date, dcc)
                   );
    }      
  }

  dClaim *= dNominal;

  return dClaim;
}

double GetClaim(Date date, 
                const BondTerms& bondTerms,
                const SessionData& sessionData,
                const shared_ptr<Numeraire>& pNumeraire)
{ 
  double dYield, dClaim;
  
  bool IsAccreting = false;

  dYield = -1.;

  shared_ptr<CashFlowStream> 
    pCashFlowStream = ComputeCouponRates(bondTerms, sessionData, pNumeraire);

  if ( bondTerms.IsAccretingBond() )
  {
    dYield = bondTerms.GetYieldToMaturityOfAccretingBond();
    IsAccreting = true;
  }
  else if ( bondTerms.IsCashPayToZeroBond() )
  {
    dYield = bondTerms.GetAccretionRateOfCashPayToZero();
    
    if ( date >= pCashFlowStream->GetAll().back().GetDate() )
      IsAccreting = true;
  }
  
  if ( !IsAccreting )
  {
    double dAccrued = 0;
    if ( pCashFlowStream )
      dAccrued = pCashFlowStream->GetAccrued(date);

    if ( bondTerms.IsCashPayToZeroBond() )// CashPayToZero
      dClaim = bondTerms.GetIssuePrice() + dAccrued;
    else  // No yield given 
      dClaim = 1 + dAccrued;
  }
  else
  {
    Frequency cmpFrequency = 
      GetYieldCompoundingFrequencyEvenIfUndefined(bondTerms);

    Date::DayCountConvention 
      dcc = GetYieldDayCountConventionEvenIfUndefined(bondTerms);

    double dRedemptionValue = bondTerms.GetRedemptionPrice();

    double dDiscount = 1. + dYield / cmpFrequency;
    
    dClaim = dRedemptionValue 
           * pow( dDiscount, 
                  cmpFrequency * 
                  Date::YearsDiff( bondTerms.GetMaturityDate(), date, dcc )
                );

    if ( pCashFlowStream )
    {
      CashFlowStream::Elements::const_iterator iter;
        
      for (iter = pCashFlowStream->begin(); 
           iter != pCashFlowStream->end();
           ++iter)
      {
        if ( iter->first > date )
          dClaim += iter->second
                  * pow( dDiscount, 
                         cmpFrequency * Date::YearsDiff(iter->first, date, dcc)
                       );
      }      
    }

  }

  return dClaim * bondTerms.GetNominal();
}

double ComputeBondLikeCallPrice(Date callDate, 
                                const BondTerms& bondTerms,
                                const SessionData& sessionData,
                                const CallSchedule& calls,
                                const shared_ptr<Numeraire>& pNumeraire)
{ 
  CallSchedule::Elements::const_iterator iterCall;
  CallSchedule::Elements periods = calls.GetAll();
      
  for (iterCall = periods.begin(); iterCall != periods.end(); ++iterCall)
  {
    if ( (*iterCall)->GetStartDate() == (*iterCall)->GetEndDate() )
    {
      // mono date call
      if ( callDate < (*iterCall)->GetEndDate() )
        break;
    }
    else 
    {
        if ( (*iterCall)->GetStartDate() <= callDate
            && callDate < (*iterCall)->GetEndDate() )
          break;
    }
  }

  if ( iterCall == periods.end() )
  {
    // The date may be end of a call period which is followed by a
    // period without call
    for (iterCall = periods.begin(); iterCall != periods.end(); ++iterCall)
    {
      if ( callDate == (*iterCall)->GetEndDate() )
        break;
    }
  }

  CHECK_COND( iterCall != periods.end(), ITO33_BONDLIKE_NO_CALL );
   
  double
    dNominal,
    dCouponAmount,
    dClaim,
    dCallPrice;

  dNominal = bondTerms.GetNominal();

  // Gets the coupon (if there are one) at the given date
  dCouponAmount = 0.;
  shared_ptr<CashFlowStream> 
    pCashFlowStream = ComputeCouponRates(bondTerms, sessionData, pNumeraire);

  if ( pCashFlowStream )
  {
    CashFlowStream::Elements::const_iterator iter;
      
    for (iter = pCashFlowStream->begin(); 
         iter != pCashFlowStream->end() && iter->first != callDate;
         ++iter);
    
    if ( iter != pCashFlowStream->end() )
      dCouponAmount = iter->second * dNominal;
  }
  
  // Compute the call value at the given date
  if ( (*iterCall)->HasYield() )
  {
    dClaim = GetClaimFromYield(callDate, (*iterCall)->GetGuaranteedYield(), 
                               bondTerms, sessionData, pNumeraire);
    dCallPrice = dClaim + dCouponAmount;
  }
  else
  {
    double 
      dAccruedAmount = 0.,
      dBase;

    if ( pCashFlowStream )
      dAccruedAmount = pCashFlowStream->GetAccrued(callDate) * dNominal;
    
    // strike is applied to principal = claim minus accrued
    dClaim = GetClaim(callDate, bondTerms, sessionData, pNumeraire);
    dBase = dClaim - dAccruedAmount;
    dCallPrice = (*iterCall)->GetStrike() * dBase;

    // accrued interest and forfeit coupon only apply for strikes
    if ( calls.GetKeepAccrued() )
      dCallPrice += dAccruedAmount;

    if ( !calls.GetForfeitCoupon() )
      dCallPrice += dCouponAmount;
  }

  return dCallPrice;
}

Frequency GetYieldCompoundingFrequencyEvenIfUndefined(const BondTerms& terms)
{
  if ( terms.GetCashDistribution() )
    return terms.GetCashDistribution()->GetPaymentFrequency();
  else if ( terms.GetFloatingRates() )
    return terms.GetFloatingRates()->GetPaymentFrequency();
  else
  {
    Frequency freq = Frequency_Undefined;
    try
    {
      freq = terms.GetYieldCompoundingFrequency();
    }
    catch(...)
    {
    }
    return freq;
  }
}

Date::DayCountConvention 
GetYieldDayCountConventionEvenIfUndefined(const BondTerms& terms)
{
  if ( terms.GetCashDistribution() )
    return terms.GetCashDistribution()->GetDayCountConvention();
  else if ( terms.GetFloatingRates() )
    return terms.GetFloatingRates()->GetDayCountConvention();
  else
  {
    Date::DayCountConvention dcc = Date::DayCountConvention_Max;
    try
    {
      dcc = terms.GetYieldDayCountConvention();
    }
    catch(...)
    {
    }
    return dcc;
  }
}

void CheckYieldCompoundingFrequency
    (
      const shared_ptr<PutSchedule>& pPutSchedule, 
      const shared_ptr<CallSchedule>& pCallSchedule,
      const BondTerms& bondTerms
    )
{
  if ( IsValid(GetYieldCompoundingFrequencyEvenIfUndefined(bondTerms)) )
    return;

  CHECK_COND(!pPutSchedule || !pPutSchedule->HasYield(),
             ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY);

  CHECK_COND(!pCallSchedule || !pCallSchedule->HasYield(), 
             ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY);
}

void CheckYieldDayCountConvention
    (
      const shared_ptr<PutSchedule>& pPutSchedule, 
      const shared_ptr<CallSchedule>& pCallSchedule,
      const BondTerms& bondTerms
    )
{
  if ( IsValid(GetYieldDayCountConventionEvenIfUndefined(bondTerms)) )
    return;

  CHECK_COND(!pPutSchedule || !pPutSchedule->HasYield(),
             ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC);

  CHECK_COND(!pCallSchedule || !pCallSchedule->HasYield(), 
             ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC);
}

} // namespace finance

} // namespace ito33
