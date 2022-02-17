/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cashflows.h
// Purpose:     cashflows class
// Author:      Nabil
// Created:     2004/03/15
// RCS-ID:      $Id: cashflows.h,v 1.20 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_CASHFLOWS_H_
#define _ITO33_PRICING_CASHFLOWS_H_

#include "ito33/sharedptr.h"

namespace ito33
{
 
namespace finance
{
  class ITO33_DLLDECL CashFlowStream;
}

namespace pricing
{

/// coupon cashflow paying accrued
class CashFlows
{

public:

  CashFlows() { DefaultInit(); }

  /**
    Creates by financial data.

    @param pCoupons cash flow stream, values are expressed as precentage
           of unpaid nominal.
    @param dNominal Initial nominal
    @param pRepaymentDates at each repayment date, a fraction of initial
           nominal is paid. Note that most of bond is redeemed only once
           at maturity.
    */
  CashFlows(const shared_ptr<finance::CashFlowStream>& pCoupons, 
            double dNominal);

  ~CashFlows()          
  { 
    delete [] m_pdAmounts;
    delete [] m_pdTimes;
    m_pdAmounts = 0;
    m_pdTimes = 0;    
  }
   
  double GetTime(size_t nIdx) const 
  {
    ASSERT_MSG( m_nNbCashFlows > nIdx, "Cash flow index out of range!" );

    return m_pdTimes[nIdx]; 
  }    

  double GetAmount(size_t nIdx) const 
  { 
    ASSERT_MSG( m_nNbCashFlows > nIdx, "Cash flow index out of range!" );

    return m_pdAmounts[nIdx]; 
  }
  
  size_t GetNbCashFlows() const { return m_nNbCashFlows; }

  double GetAccruedInterest(double dTime, bool bPlus = true) const;


protected:
  
  void DefaultInit()
  {
    m_nNbCashFlows = 0;
    m_pdAmounts    = 0;
    m_pdTimes      = 0;
  }
  
  size_t m_nNbCashFlows;

  double *m_pdAmounts;
  
  double *m_pdTimes;
};


} // namespace pricing

} // namespace ito33

#endif  //#ifndef _ITO33_PRICING_CASHFLOWS_H_
