/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/gammavarianceswap.h
// Purpose:     financial class for gamma variance swaps
// Created:     2006/08/04
// RCS-ID:      $Id: gammavarianceswap.h,v 1.4 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/gammavarianceswap.h
    @brief Declaration of the financial gamma variance swap class.
 */

#ifndef _ITO33_FINANCE_GAMMAVARIANCESWAP_H_
#define _ITO33_FINANCE_GAMMAVARIANCESWAP_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/varianceswaplike.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL VarianceSwapTerms;
class ITO33_DLLDECL SessionData;

/**
  GammaVarianceSwap represents the financial aspects of a gamma variance swap.
*/
class ITO33_DLLDECL GammaVarianceSwap : public VarianceSwapLike
{
public:
  /**
      Creates a gamma variance swap by taking the terms and the volatility 
      strike.
       
      @param pTerms The terms of the swap
      @param dVolatilityStrike The volatility strike of the swap
   */
  GammaVarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
                    double dVolatilityStrike);

  /// Empty virtual destructor.
  virtual ~GammaVarianceSwap() { }


  /// @name Modifiers for gamma variance swap.
  //@{

  /**
      Sets user computed values at the valuation date.
     
      These values are only relevant if the valuation date is past the start 
      of the sampling period of the contract.  
      
      The current gamma volatility defaults to zero.  Volatility calculations
      must not include the drift term, and should use the return method 
      specified in the variance swap terms.  It must also include the S_i / S_0
      gamma factor.

      The number of samples used should be less than the number of
      sampling returns specified in the variance swap terms.

      The start share price is the share price at the start of the
      sampling period (S_0).  If the valuation date equals the sampling
      start date, the start share price is set to the previous share price.

      @param dCurrentGammaVolatility current realized gamma volatility
      @param nNbSamplesUsed number of samples used in the current volatility
      @param dStartSharePrice the share price at the sampling start date
   */
  void SetCurrentValues(double dCurrentGammaVolatility, 
                        size_t nNbSamplesUsed,
                        double dStartSharePrice);

  //@}

  /// @name Accessors for gamma variance swap.
  //@{

  /**
      The share price at the start of the sampling period.

      @return the share price at the sampling start date
   */
  double GetStartSharePrice() const
  {
    return m_dStartSharePrice;
  }

  //@}

  virtual void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


protected:

  /// share price at the start of the sampling period
  double m_dStartSharePrice;

}; // class GammaVarianceSwap


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_GAMMAVARIANCESWAP_H_
