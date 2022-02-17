/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/conditionalvarianceswap.h
// Purpose:     financial class for conditional variance swaps
// Created:     2006/08/04
// RCS-ID:      $Id: conditionalvarianceswap.h,v 1.3 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/conditionalvarianceswap.h
    @brief Declaration of the financial conditional variance swap class.
 */

#ifndef _ITO33_FINANCE_CONDITIONALVARIANCESWAP_H_
#define _ITO33_FINANCE_CONDITIOANALVARIANCESWAP_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/varianceswaplike.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL VarianceSwapTerms;
class ITO33_DLLDECL SessionData;

/**
   ConditionalVarianceSwap represents the financial aspects of a conditional
   variance swap.
*/
class ITO33_DLLDECL ConditionalVarianceSwap : public VarianceSwapLike
{
public:

  /**
      Creates a conditional variance swap by taking the terms and the 
      volatility strike.
   
      @param pTerms The terms of the swap
      @param dVolatilityStrike The volatility strike of the swap
   */
  ConditionalVarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
                          double dVolatilityStrike);

  /// Empty virtual destructor.
  virtual ~ConditionalVarianceSwap() { }


  /// @name Modifiers for conditional variance swap.
  //@{

  /**
      Sets user computed values at the valuation date.
     
      These values are only relevant if the valuation date is past the start 
      of the sampling period of the contract.  
      
      The current volatility defaults to zero.  Volatility calculations 
      must not include the drift term, must include the corridor, and 
      should use the return method specified in the variance swap terms.

      The number of samples used should be less than the number of
      sampling returns specified in the variance swap terms.

      The number of days within the corridor should be less than the number 
      of sampling returns specified in the variance swap terms.

      @param dCurrentVolatility current realized volatility
      @param nNbSamplesUsed number of samples used in the current volatility
      @param iCurrentConditionalCount the number of days within the corridor
   */
  void SetCurrentValues(double dCurrentConditionalVolatility, 
                        size_t nNbSamplesUsed,
                        int iCurrentConditionalCount);

  //@}

  /// @name Accessors for variance swap.
  //@{

  /**
      The number of days the share price has been within the corridor
      between the sampling start date and the valuation date.
     
      @return the number of days within the corridor
   */
  int GetCurrentConditionalCount() const
  {
    return m_iCurrentConditionalCount;
  }

  //@}

  virtual void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


protected:

  /// number of days within the corridor from start date to valuation date
  int m_iCurrentConditionalCount;

}; // class ConditionalVarianceSwap


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_CONDITIONALVARIANCESWAP_H_
