/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cboptiondata.h
// Purpose:     cb option data class
// Author:      Nabil
// Created:     2005/06/10
// RCS-ID:      $Id: cboptiondata.h,v 1.4 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cboptiondata.h
    @brief The cb option data class.

    The class for the data of the CB option contracts.  
 */

#ifndef _ITO33_PRICING_CBOPTIONDATA_H_
#define _ITO33_PRICING_CBOPTIONDATA_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/finance/bondlike/cboption.h"

#include "ito33/pricing/cashflows.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL CBOption;
  class ITO33_DLLDECL FloatingRates;
}

namespace pricing
{

/// The declaration of the cb option data class.
class CBOptionData
{

public:

  /**
     Creates a CBOptionData by financial CBOption object.

     @param cboption financial CBOption
   */
  CBOptionData(const finance::CBOption& cboption);

  // default dtor is ok
  
  /// @name Some "sets" functions
  //@{
      //TODO: for make whole...
  //@}
  
  /// @name Some "gets" functions
  //@{

  double GetRedemptionRate() const { return m_dRedemptionRate; }

  double GetCbOptionFaceValue() const { return m_dCbOptionFaceValue; }
  
  double GetMaturityTime() const { return m_dMaturityTime; }
  
  double GetBalloonCoupon() const { return m_dBalloonCoupon; }

  /**
     Gets the fixed cash flows in amount
    
     @return the fixed cash flows in amount
   */
  const shared_ptr<CashFlows>& GetFixedCashFlows() const 
  {
    return m_pFixedCashFlows;
  }
  
  /**
     Gets the floating cash flows in amount
    
     @return the floating cash flows in amount
   */
  const shared_ptr<CashFlows>& GetFloatingCashFlows() const 
  {
    return m_pFloatingCashFlows;
  }

  //@}

protected:

  /// The maturity time of the cb option
  double m_dMaturityTime;

  /// redemption amount of the cb option (used in the fixed leg)
  double m_dRedemptionRate;  
  
  /// The face value of the cb (used in the fixed leg)
  double m_dCbFaceValue;
  
  /// The face value of the cb option (used in the floating leg)
  double m_dCbOptionFaceValue;

  /// The balloon coupon
  double m_dBalloonCoupon;

  /// fixed cash flow IN AMOUNT
  shared_ptr<CashFlows> m_pFixedCashFlows;
  
  /// floating rates
  shared_ptr<finance::FloatingRates> m_pFloatingRates;

  /// floating cash flow IN AMOUNT
  shared_ptr<CashFlows> m_pFloatingCashFlows;

}; // class CBOptionData;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBOPTIONDATA_H_
