/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/reset.h
// Purpose:     contract class for reset (backward)
// Author:      David and Yann
// Created:     2004/11/03
// RCS-ID:      $Id: reset.h,v 1.7 2006/02/27 19:26:03 yann Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/reset.h
    @brief The declaration of the reset contract class.

    The base class for reset contracts.  
 */

#ifndef _ITO33_PRICING_RESET_H_
#define _ITO33_PRICING_RESET_H_

#include "ito33/pricing/cb.h"
#include "ito33/finance/bondlike/resetflooredby.h"

namespace ito33
{

namespace finance
{
  class Reset;
}

namespace pricing
{

/// The declaration of the (backward) cb contract class.
class Reset : public CB 
{

public:

  /**
     Creates a CB by financial Reset object.

     @param reset financial Reset
   */
  Reset(const finance::Reset& reset);

  // default ctor
  Reset() : CB() {}

  /// virtual dtor for base class
  virtual ~Reset() { }
    
  /// @name reset functions
  //@{
 
  /**
     Gets the conversion cap rates.
  
     @return a vector of conversion cap rates
   */
  const std::vector<double>& GetCapRates() const
  { 
    return m_pdCapRates; 
  }

  /**
     Get the conversion floor rates.

     @return a vector of conversion floor rates
   */
  const std::vector<double>& GetFloorRates() const
  { 
    return m_pdFloorRates; 
  }

  /**
     Get the multipliers.

     @return a vector of reset multipliers
   */
  const std::vector<double>& GetMultipliers() const 
  { 
    return m_pdMultipliers; 
  }

  /**
     Gets the method for computing the conversion price floor

     @return the method for computing the conversion price floor
   */
  finance::ResetFlooredBy GetFlooredBy() const 
  { 
    return m_FlooredBy; 
  }

  /**
     Gets the reset times
  
     @return a vector of reset times
   */
  const std::vector<double>& GetResetTimes() const 
  { 
    return m_pdResetTimes; 
  } 

  /**
     Gets current conversion price.

     @return current conversion price
   */
  double GetCurrentConversionPrice() const
  { 
    return m_dCurrentConvPrice; 
  }

  /**
     Gets initial conversion price.

     @return initial conversion price
   */
  double GetInitialConversionPrice() const
  { 
    return m_dInitialConvPrice; 
  }

  /**
     For 1D pricing, the variable is a function of the conversion ratio k.
     For normal path dep pricing, we are capping and flooring the
     conversion price N/k.  Thus, for 1D, need to flip the caps and floors.
   */
  void FlipFloorAndCap();

  /**
     Sets the current conversion price.

     For path dependent pricing, the conversion price along each path is 
     different.  When these ratios are set, the current conversion price
     must also be reset to be consistent.
     
     @param dCurrentConvPrice the current conversion price
   */
  void SetCurrentConversionPrice(double dCurrentConvPrice)
  {
    m_dCurrentConvPrice = dCurrentConvPrice;
  }

  //@}


protected:

  /// the reset times
  std::vector<double> m_pdResetTimes; 

  /// reset cap rate
  std::vector<double> m_pdCapRates;

  /// reset floor rate
  std::vector<double> m_pdFloorRates;

  /// reset multiplier value
  std::vector<double> m_pdMultipliers;

  /// how the floor for resets is computed
  finance::ResetFlooredBy m_FlooredBy;

  /// the current conversion price.
  double m_dCurrentConvPrice;

  /// the initial conversion price.
  double m_dInitialConvPrice;

}; // class Reset;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_RESET_H_
