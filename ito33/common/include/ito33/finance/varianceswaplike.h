/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/varianceswaplike.h
// Purpose:     financial base class for variance swap contracts
// Created:     2006/08/04
// RCS-ID:      $Id: varianceswaplike.h,v 1.6 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/varianceswaplike.h
    @brief Declaration of the base financial variance/volatility swap class.
 */

#ifndef _ITO33_FINANCE_VARIANCESWAPLIKE_H_
#define _ITO33_FINANCE_VARIANCESWAPLIKE_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL VarianceSwapTerms;
class ITO33_DLLDECL SessionData;

/**
    VarianceSwapLike represents the common aspects of all variance swaps.

    @nocreate
 */
class ITO33_DLLDECL VarianceSwapLike : public Derivative
{
public:
  
  // protected constructor

  /// Empty virtual destructor.
  virtual ~VarianceSwapLike() { }


  /// @name Modifiers for VarianceSwapLike.
  //@{

  /**
      @internal
      @brief Sets the volatility strike of the variance swap.

      @param dVolatilityStrike the new volatility strike

      @noexport
   */
  void SetVolatilityStrike(double dVolatilityStrike)
  {
    m_dVolatilityStrike = dVolatilityStrike;
  }

  //@}

  /// @name Accessors for VarianceSwapLike.
  //@{

  /**
      Gets the terms of the swap.

      @return the terms of the swap
   */
  const shared_ptr<VarianceSwapTerms>& GetTerms() const { return m_pTerms; }

  /**
      Gets the volatility strike.

      @return the volatility strike of the swap
   */
  double GetVolatilityStrike() const
  {
    return m_dVolatilityStrike;
  }

  /**
      Gets the maturity date of the swap.

      @return maturity date of the swap
   */
  Date GetMaturityDate() const;

  /**
      The current realized volatility.

      Must be defined in relation to the variance swap payoff. For example,
      a gamma payoff must return the volatility with the extra S_i / S_0 
      factor.

      @return the current realized volatility
   */  
  double GetCurrentVolatility() const { return m_dCurrentVolatility; }

  /**
      The number of samples used to compute the current realized volatility.

      @return the number of samples used to compute the current volatility
   */
  size_t GetNbSamplesUsed() const { return m_nNbSamplesUsed; }

  //@}
  
  /**
      Sets the market price for the VS.
 
      Note: Re-implementation of virtual base class since VS prices can be
      positive or negative.

      @param dPrice market price
   */
  virtual void SetMarketPrice(double dPrice) 
  {
    DoSetMarketPrice(dPrice);
  }

  // base class virtual functions
  virtual void ValidateWith(const SessionData& sessionData) const;


protected:

  /**
      Constructor (only called by derived classes).
      
      All variance swaps must specify the variance swap terms and the 
      volatility strike.  If the payoff is based on variance, the volatility
      strike is squared when computing the payoff value. If the start of the 
      sampling period is T_0, return calculations begin at T_1.  The return 
      at T_1 will use the spot values S_1 and S_0.  A nominal of 1 is assumed.
   
      @param pTerms The terms of the swap
      @param dVolatilityStrike The volatility strike of the swap
   */
  VarianceSwapLike(const shared_ptr<VarianceSwapTerms>& pTerms,
                   double dVolatilityStrike);

#ifndef __CPP2ANY__

  /**
      Sets user computed values at the valuation date.
     
      These values are only relevant if the valuation date is past the start 
      of the sampling period of the contract.  
      
      The current volatility defaults to zero. Volatility calculations 
      must not include the drift term, and should use the return method 
      specified in the variance swap terms.

      The number of samples used should be less than the number of
      sampling returns specified in the variance swap terms.

      If current volatility is set, all current values will be set together
      by the derived class.

      @param dCurrentVolatility current realized volatility
      @param nNbSamplesUsed number of samples used in the current volatility
   */
  void SetCurrentValues(double dCurrentVolatility, size_t nNbSamplesUsed);

#endif

  /**
      Checks if current values are set.

      If current volatility is set, all current values will be set together
      by the derived class.

      @return true if current values are set, false otherwise
   */
  bool HasCurrentValues() const { return m_dCurrentVolatility >= 0.; }

  /**
      Dumps entries common to all variance swaps.

      @param tagParent tag of the derived class
   */
  void DumpMe(XML::Tag& tagParent) const;


  /// common terms of the variance swap
  shared_ptr<VarianceSwapTerms> m_pTerms;

  /// volatility strike
  double m_dVolatilityStrike;

  /// current volatility (if valuation is past the sampling start)
  double m_dCurrentVolatility;

  /// number of samples used to compute the current volatility
  size_t m_nNbSamplesUsed;

}; // class VarianceSwapLike


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_VARIANCESWAPLIKE_H_
