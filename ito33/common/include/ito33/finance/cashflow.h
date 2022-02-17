///////////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cashflow.h
// Purpose:     Class for description of one payment in a cash flow streams
// Author:      Pedro Ferreira
// Created:     Mar 24, 2004
// RCS-ID:      $Id: cashflow.h,v 1.32 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cashflow.h
    @brief Class for description of one payment in a cash flow streams
 */

#ifndef _ITO33_FINANCE_CASHFLOW_H_
#define _ITO33_FINANCE_CASHFLOW_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{


/**
    This class models one element of the cash flow stream that is used by
    derivatives like credit default swaps or convertible bonds.

    @nocreate
 */
class ITO33_DLLDECL CashFlow
{
public:

  // default ctors, assignment operator are ok

  /**
     Constructor initializes all the members.
   
     @param dt the date for this cash flow.
     @param amt the amount to be payed at the specified date.
   */
  CashFlow(Date dt, double amt) : m_date(dt), m_dAmount(amt) 
  {
  }

  // Default dtor is ok

  /** 
     Gets the date for this cash flow.
   
     @return the date for the payment.
   */
  Date GetDate() const;

  /** 
     Gets the amount for this cash flow.
   
     @return the amount for the payment.
   */
  double GetAmount() const;

private:

  /// Date for the payment
  Date m_date;      
  
  /// Amount for the payment
  double m_dAmount;  
};

// Inline functions

/// Get Date function
inline Date CashFlow::GetDate() const
{
  return m_date;
}

/// Get amount function
inline double CashFlow::GetAmount() const
{
  return m_dAmount;
}

} // namespace finance

} // namespace ito33

namespace std
{

/**
   Comparison operator for two CashFlows.
   
   @param x first element to be compared.
   @param y second element to be compared.
   
   @return true if the date of x is less than the date of y.
 */
inline bool operator<(const ito33::shared_ptr<ito33::finance::CashFlow>& x,
                      const ito33::shared_ptr<ito33::finance::CashFlow>& y)
{
  return x->GetDate() < y->GetDate();
}

/**
   Comparison operator for two CashFlows.
   
   @param x first element to be compared.
   @param y second element to be compared.
   
   @return true if the date of x is equal to the date of y.
 */
inline bool operator==(const ito33::shared_ptr<ito33::finance::CashFlow>& x,
                       const ito33::shared_ptr<ito33::finance::CashFlow>& y)
{
  return x->GetDate() == y->GetDate();
}

} // namespace std
 
#endif // #ifndef _ITO33_FINANCE_CASHFLOW_H_
