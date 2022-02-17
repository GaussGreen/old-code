/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/varianceswap.h
// Purpose:     contracts class for variance swap pricing (backward)
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswap.h,v 1.17 2006/08/16 19:09:32 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/varianceswap.h
    @brief The declaration of the contracts class for variance swaps.  
 */

#ifndef _ITO33_PRICING_VARIANCESWAP_H_
#define _ITO33_PRICING_VARIANCESWAP_H_

#include "ito33/common.h"

#include "ito33/finance/swaptype.h"
#include "ito33/finance/returntype.h"
#include "ito33/finance/swappayofftype.h"

#include "ito33/pricing/contract.h"


namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL VarianceSwapLike;
  class ITO33_DLLDECL VarianceSwap;
  class ITO33_DLLDECL GammaVarianceSwap;
  class ITO33_DLLDECL ConditionalVarianceSwap;
  class ITO33_DLLDECL OptionVarianceSwap;
}

namespace pricing
{


/// The declaration of the (backward) variance swap contract class.
class VarianceSwap : public Contract
{
public:

  /**
      The ctor for a standard variance swap.
    
      @param varianceSwap reference to an object of type finance::VarianceSwap
   */
  VarianceSwap(const finance::VarianceSwap& varianceSwap);

  /**
      The ctor for a gamma variance swap.
    
      @param varianceSwap reference to an object of type 
                          finance::GammaVarianceSwap
   */
  VarianceSwap(const finance::GammaVarianceSwap& varianceSwap);

  /**
      The ctor for a conditional variance swap.
    
      @param varianceSwap reference to an object of type 
                          finance::ConditionalVarianceSwap
   */
  VarianceSwap(const finance::ConditionalVarianceSwap& varianceSwap);

  /**
      The ctor for an option variance swap.
    
      @param varianceSwap reference to an object of type 
                          finance::OptionVarianceSwap
   */
  VarianceSwap(const finance::OptionVarianceSwap& varianceSwap);

  /**
      Empty virtual destructor.
   */
  virtual ~VarianceSwap() { }

  /**
      @name Accessors variance swap.
   */
  //@{
 
  /**
      Checks if the variance swap is vanilla, so no cap, no barrier and is
      of type variance.

      @return true if the variance swap is vanilla, false otherwise
   */
  bool IsVanilla() const
  {
    return    m_dCapMultiplier < 0 
           && m_dUpCorridorBarrier == 0.0 && m_dDownCorridorBarrier == 0.0 
           && m_swapType == finance::Swap_Variance;
  }

  /**
      Gets the volatility strike.

      @return the volatility strike
   */
  double GetVolatilityStrike() const
  {
    return m_dVolatilityStrike;
  }

  /**
      Gets the start time of the sampling period.

      @return the start time of the sampling period
   */  
  double GetStartTimeOfSamplingPeriod() const
  {  
    return m_dStartTimeOfSamplingPeriod;   
  }

  /** 
      Gets the cap multiplier
       
      @return the cap multiplier
  */
  double GetCapMultiplier() const
  {
    return m_dCapMultiplier;
  }

  /**
      Gets the number of sampling days/returns.

      @return the number of sampling days/returns
   */
  size_t GetNbSamplingReturns() const
  {
    return m_nNbSamplingReturns;
  }

  /**
      Gets the swap type (variance or volatility).

      @return the swap type
   */
  finance::SwapType GetSwapType() const
  {
    return m_swapType;
  }

  /**
      Gets the return calculation method (actual or log).

      @return the return calculation mathod
   */
  finance::ReturnType GetReturnType() const
  {
    return m_returnType;
  }

  /**
      Gets the current volatility.     

      @return the current realized volatility
   */  
  double GetCurrentVolatility() const
  {  
    return m_dCurrentVolatility;
  }

  /**
      Gets the current average of the squared returns.

      @return the current average of the squared returns
   */
  double GetCurrentAvgSqrReturn() const
  {
    return m_dCurrentAvgSqrReturn;
  }

  /**
      Gets the number of samples used to compute the current volatility.

      @return the number of samples used to compute the current volatility
   */
  size_t GetNbSamplesUsed() const
  {
    return m_nNbSamplesUsed;
  }

  /**
      Gets the current number of days within the corridor for the conditional
      payoff.     

      @return the current number of days spent within the corridor
   */  
  int GetCurrentConditionalCount() const
  {  
    return m_iCurrentConditionalCount;
  }

  /**
      Compute payoff value for the specified average of the squared returns.

      Z is the pricing state variable for the average of the squared returns.
      It is equal to Y X / N, where X is defined in the var swap documentation,
      N is the total number of samples, and Y is possibly some factor to 
      account for exotic variance swap types.

      @param dZ value at which we want the payoff

      @return value of the payoff
   */
  double GetPayoffValue(double dZ) const;

  /**
      Gets the previous share price.

      @return the previous share price
   */
  double GetPreviousSharePrice() const
  {
    return m_dPreviousSharePrice;
  }
    
  /**  
      Gets the number of sampling returns per year. 
      
      By default it is set to 252.

      @return the number of sampling returns per year.
   */
  size_t GetAnnualReturnFrequency() const
  {
    return m_nAnnualReturnFrequency;
  }

  /**  
      Gets the up corridor barrier.

      @return the up corridor barrier
   */
  double GetUpCorridorBarrier() const
  {
    return m_dUpCorridorBarrier;
  }

  /**  
      Gets the down corridor barrier.

      @return the down corridor barrier
   */
  double GetDownCorridorBarrier() const
  {
    return m_dDownCorridorBarrier;
  }

  /**  
      Gets the payoff type.

      @return the payoff type
   */
  finance::SwapPayoffType GetSwapPayoffType() const
  {
    return m_swapPayoffType;
  }

  /**
      Gets the share price at the start of the sampling period.

      @return the share price at the sampling start date
   */
  double GetStartSharePrice() const
  {
    return m_dStartSharePrice;
  }
  
  //@}

protected:

  /**
      Initializes terms common to all variance swap types.

      @param varianceSwap reference to an object of type finance::VarianceSwapLike
   */
  void InitCommonTerms(const finance::VarianceSwapLike& varianceSwap);

  /// volatility strike
  double m_dVolatilityStrike;

  /// the start time of the sampling period (converted from date)
  double m_dStartTimeOfSamplingPeriod;

  /// number of sampling days/returns
  size_t m_nNbSamplingReturns;

  /// swap type (variance or volatility)
  finance::SwapType m_swapType;

  /// return calculation method (actual or log)
  finance::ReturnType m_returnType;

  /// swap payoff type (gamma, conditional, etc)
  finance::SwapPayoffType m_swapPayoffType;

  /// current volatility
  double m_dCurrentVolatility;

  /// current average of the squared returns
  double m_dCurrentAvgSqrReturn;

  /// the number of samples used to compute the current volatility
  size_t m_nNbSamplesUsed;

  /// the previous share price
  double m_dPreviousSharePrice;

  /// cap multiplier
  double m_dCapMultiplier;
  
  /// number of sampling returns in a year
  size_t m_nAnnualReturnFrequency;

  /// up corridor barrier
  double m_dUpCorridorBarrier;

  /// down corridor barrier
  double m_dDownCorridorBarrier;

  /// share price at the start of the sampling period
  double m_dStartSharePrice;

  /// number of days within the corridor from start date to valuation date
  size_t m_iCurrentConditionalCount;

}; // class VarianceSwap;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_VARIANCESWAP_H_
