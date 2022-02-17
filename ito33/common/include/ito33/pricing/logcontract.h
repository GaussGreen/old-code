/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/logcontract.h
// Purpose:     contracts class for LogContract
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract.h,v 1.1 2006/07/19 17:38:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/logcontract.h
    @brief The declaration of the pricing log contract class.
 */

#ifndef _ITO33_PRICING_LOGCONTRACT_H_
#define _ITO33_PRICING_LOGCONTRACT_H_

#include "ito33/common.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL LogContract;
}

namespace pricing
{


/// The declaration of the log contract class.
class LogContract : public Contract 
{
public:

  /**
      The ctor.
    
      @param logContract The financial LogContract object
   */
  LogContract(const finance::LogContract& logContract);
  
  /// The start of return period
  double GetT0() const { return m_dT0; }

  /// The spot at the start of return period
  double GetS0() const { return m_dS0; }

private:

  /// The start of return period
  double m_dT0;

  /// The spot at the start of return period
  double m_dS0;

}; // class LogContract;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_LOGCONTRACT_H_
