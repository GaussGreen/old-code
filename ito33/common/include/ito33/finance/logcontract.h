/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/logcontract.h
// Purpose:     financial log contract class
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract.h,v 1.1 2006/07/19 17:38:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/logcontract.h
    @brief declaration of the log contract
 */

#ifndef _ITO33_FINANCE_LOGCONTRACT_H_
#define _ITO33_FINANCE_LOGCONTRACT_H_

#include "ito33/date.h"

#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{

/**
    Log contract is a virtual contract that pays log(S_T/S_0) at maturity,
    where S_T is the spot at the maturity, S_0 is the spot at the start of
    return period.
 */
class ITO33_DLLDECL LogContract : public Derivative
{
public:
  /**
      Creates a LogContract object.

      @param maturityDate the maturity date of the log contract
      @param startOfReturnPeriod the date at which the spot will be used for
                                 payoff
   */
  LogContract(Date maturityDate, Date startOfReturnPeriod);
  
  /**
      The share price at the start of return period if the latter is 
      already passed.

      @param dS0 The share price at the start of return period, required if
                 the latter is already passed.
   */
  void SetStartSharePrice(double dS0);

  /**
      Gets the maturity date of the LogContract.

      @return Maturity date of the LogContract
   */
  Date GetMaturityDate() const
  {
    return m_maturityDate;
  }

  /**
      Gets the start of return period of the LogContract.

      @return The start of return period of the LogContract
   */
  Date GetStartOfReturnPeriod() const
  {
    return m_startOfReturnPeriod;
  }

  /**
      The share price at the start of return period if the latter is 
      already passed.

      @return The share price at the start of return period, required if
              the latter is already passed.
   */
  double GetStartSharePrice() const { return m_dS0; }

  // Implement pure virtual functions in base class
  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual XML::Tag Dump(XML::Tag& tagParent) const;
 

private:

  /// Maturity date of the contract
  Date m_maturityDate;

  /// The date at which the spot will be used for payoff
  Date m_startOfReturnPeriod;

  /**
      The share price at the start of return period if the latter is 
      already passed.
   */
  double m_dS0;

}; // class LogContract


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_LOGCONTRACT_H_
