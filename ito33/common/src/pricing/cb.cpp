/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cb.cpp
// Author:      Nabil
// Created:     2004/03/15
// RCS-ID:      $Id: cb.cpp,v 1.104 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/floatingrates.h"
#include "ito33/finance/payoffconstant.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/cb_base.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/utils.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbcalls.h"
#include "ito33/pricing/cbconversions.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(pricing::CB);

  using namespace numeric;
  using namespace finance;

namespace pricing
{

  
void CB::GetBondData(const SessionData& sessionData, 
                     const BondTerms& bondTerms,
                     const shared_ptr<CallSchedule>& pCalls,
                     const shared_ptr<PutSchedule>& pPuts,
                     const shared_ptr<Numeraire>& pNumeraire)
{
  GetBondLikeTermsData(bondTerms, sessionData, pNumeraire);

  m_dRedemptionValue = m_dNominal * bondTerms.GetRedemptionPrice();

  shared_ptr<finance::CashFlowStream> pcashFlowStream;

  if ( bondTerms.GetFloatingRates() )
  {
    // Compute floating payments.
    pcashFlowStream = bondTerms.GetFloatingRates()
                    ->ComputeFloatingPayments( m_pDerivativeCurve );
  }
  else
    pcashFlowStream = bondTerms.GetCashDistribution();

  m_pCashFlows = make_ptr( new CashFlows
                               ( pcashFlowStream, bondTerms.GetNominal() ) );

  m_dAccretionStartTime = -100.;

  // Check if claim is computed by yield. 
  m_bHasYield = false;
  m_dYield = -1.0;
  if ( bondTerms.IsAccretingBond())
  {
    m_bHasYield = true;
    m_dYield = bondTerms.GetYieldToMaturityOfAccretingBond();
    m_dAccretionStartTime = GetDoubleFrom( bondTerms.GetIssueDate() );
  }

  // Check for partial OID. Only one of OIDYield and AccretionRateOfCashPayToZero
  // should be set in bond terms
  if ( bondTerms.IsCashPayToZeroBond() )
  {
    m_bHasYield = true;
    m_dYield = bondTerms.GetAccretionRateOfCashPayToZero();
    m_dAccretionStartTime = GetDoubleFrom 
                            ( pcashFlowStream->GetLastPaymentDate() );
  }

  Frequency freq = GetYieldCompoundingFrequencyEvenIfUndefined(bondTerms);
  // Compounding frequency treatment
  if( IsValid(freq) )
    m_dCmpFreq = double( freq );
  else
    m_dCmpFreq = 0;

  m_calls = CBCalls( pCalls );

  m_puts = CBPuts( pPuts );
  
  m_pPayoff = make_ptr( new finance::PayoffConstant(GetRedemptionValue()) );
  
}

void CB::GetCBBaseData(const finance::CBBase& cblike)
{  
  // Note: m_conversions must be initialized before calling this function

  GetBondData(*cblike.GetSessionData(), *cblike.GetBondTerms(),
              cblike.GetCallSchedule(), cblike.GetPutSchedule(),
              cblike.GetNumeraire() );

  GetConvertibleLikeData(cblike);
}

CB::CB(const finance::ConvertibleBond& cb) : CBLike()
{
  m_conversions = CBConversions( cb.GetConversionSchedule(),
                                 cb.GetSessionData()->GetValuationDate() );

  // note: must construct m_conversions before calling
  GetCBBaseData(cb);
}


CB::CB(const finance::Bond& bond) : CBLike()
{
  GetBondData(*bond.GetSessionData(), *bond.GetBondTerms(),
              bond.GetCallSchedule(), bond.GetPutSchedule(),
              bond.GetNumeraire() );
}

double CB::GetClaim(double dTime, bool bPlus) const
{ 
  double dClaim = 0.0;
  if( !m_bHasYield )    
    // this is a normal bond without accretion
    dClaim = GetNominal() + m_pCashFlows->GetAccruedInterest(dTime, bPlus);
  else if ( !IsAccreting(dTime, bPlus) ) 
  {
    // this is a Cash pay to zero
    dClaim = m_dIssuePrice + m_pCashFlows->GetAccruedInterest(dTime, bPlus);
  }
  else
  {
    double dDiscount = 1. + GetYield() / m_dCmpFreq;
    
    dClaim = GetRedemptionValue()
           * pow(dDiscount, m_dCmpFreq * (dTime - m_dMaturityTime));

    size_t nNbCashFlows = m_pCashFlows->GetNbCashFlows();

    // Coupons after maturity (call notice for example) shouldn't be counted
    // If there is cash flows, the first will be always before maturity
    // since it's the contracting date of the payments if the contracting date
    // is reasonable
    while (   nNbCashFlows > 1 
           && IsAfter(m_pCashFlows->GetTime(nNbCashFlows - 1), m_dMaturityTime))
    {
      nNbCashFlows--;
    }

    size_t nIdxCashFlow;

    if (bPlus) // don't count the coupon at dTime (if any)
      for (nIdxCashFlow = nNbCashFlows - 1; 
              nIdxCashFlow < nNbCashFlows  
              && numeric::IsAfter(m_pCashFlows->GetTime(nIdxCashFlow), dTime); 
          nIdxCashFlow--)
        dClaim += m_pCashFlows->GetAmount(nIdxCashFlow)* 
          pow( dDiscount, 
                  m_dCmpFreq* ( dTime - m_pCashFlows->GetTime(nIdxCashFlow) ) );
    else
      for (nIdxCashFlow = nNbCashFlows - 1; 
              nIdxCashFlow < nNbCashFlows  
          && numeric::IsEqualOrAfter(m_pCashFlows->GetTime(nIdxCashFlow),dTime); 
          nIdxCashFlow--)
        dClaim += m_pCashFlows->GetAmount(nIdxCashFlow)* 
        pow( dDiscount, 
              m_dCmpFreq*( dTime - m_pCashFlows->GetTime(nIdxCashFlow) ) );
  
  } //end if

  return dClaim;

}//CB::GetClaim(double dTime, bool bPlus)


} //namespace pricing

} //namespace ito33
