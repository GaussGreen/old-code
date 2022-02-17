///////////////////////////////////////////////////////////////////////////
// Name:        pricing/cashflows.cpp
// Purpose:     cashflows class
// Author:      Nabil
// Created:     2004/03/15
// RCS-ID:      $Id: cashflows.cpp,v 1.20 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"
#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/cashflowstream.h"

#include "ito33/pricing/cashflows.h"

namespace ito33
{

//Implementation of the AutoPtrDeleter for CashFlows class.
ITO33_IMPLEMENT_AUTOPTR(pricing::CashFlows);

}

namespace ito33
{

using finance::CashFlowStream;
using namespace numeric;

namespace pricing
{  

CashFlows::CashFlows(const shared_ptr<CashFlowStream>& pCoupons, double dNominal)
                   : m_nNbCashFlows(0), m_pdAmounts(0), m_pdTimes(0)
{
  if ( !pCoupons )
    return;  
  
  CashFlowStream::const_iterator iter = pCoupons->begin();
  
  size_t nNbCashFlows = pCoupons->size();

  if ( !nNbCashFlows )
    return;

  // Contracting date will be added to the schedule of the cash flows 
  // that is why m_nNbCashFlows = nNbCashFlows + 1.
  m_nNbCashFlows = nNbCashFlows + 1;
  
  m_pdTimes = new double [m_nNbCashFlows];
  m_pdAmounts = new double [m_nNbCashFlows];

  m_pdTimes[0]   = GetDoubleFrom( pCoupons->GetContractingDate() );
  m_pdAmounts[0] = 0;

  iter = pCoupons->begin();
  for (size_t nIdx = 0; iter != pCoupons->end(); nIdx++, ++iter)
  {
    m_pdTimes[nIdx + 1]   = GetDoubleFrom(iter->first);
    m_pdAmounts[nIdx + 1] = iter->second * dNominal;
  }

  return;
}

double CashFlows::GetAccruedInterest(double dTime, bool bPlus) const
{
  size_t
    nIdxCashFlow;
  
  if(bPlus)
    for(nIdxCashFlow = 0; nIdxCashFlow < m_nNbCashFlows && 
        IsEqualOrAfter(dTime, m_pdTimes[nIdxCashFlow]); ++nIdxCashFlow)
      ;
  else
    for(nIdxCashFlow = 0; nIdxCashFlow < m_nNbCashFlows && 
        IsAfter(dTime, m_pdTimes[nIdxCashFlow]); ++nIdxCashFlow)
      ;
  
  if(nIdxCashFlow == 0 || nIdxCashFlow == m_nNbCashFlows)
    return 0.;
  else
    return m_pdAmounts[nIdxCashFlow] * 
      (dTime - m_pdTimes[nIdxCashFlow - 1]) /
      ( m_pdTimes[nIdxCashFlow] - m_pdTimes[nIdxCashFlow - 1] );
}


} // namespace pricing

} // namespace ito33
