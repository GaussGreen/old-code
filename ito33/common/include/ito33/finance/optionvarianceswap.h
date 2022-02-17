/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/optionvarianceswap.h
// Purpose:     financial class for option variance swaps
// Created:     2006/08/04
// RCS-ID:      $Id: optionvarianceswap.h,v 1.3 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/optionvarianceswap.h
    @brief Declaration of the financial option variance swap class.
 */

#ifndef _ITO33_FINANCE_OPTIONVARIANCESWAP_H_
#define _ITO33_FINANCE_OPTIONVARIANCESWAP_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/varianceswaplike.h"
#include "ito33/finance/optiontype.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL VarianceSwapTerms;
class ITO33_DLLDECL SessionData;

/**
   OptionVarianceSwap represents the financial aspects of an option 
   variance swap.
*/
class ITO33_DLLDECL OptionVarianceSwap : public VarianceSwapLike
{
public:
  /**
      Creates an option variance swap by taking the terms, volatility 
      strike, and option type (put or call).
        
      @param pTerms The terms of the swap
      @param dVolatilityStrike The volatility strike of the swap
      @param optionType The option payoff type (call or put)
   */
  OptionVarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
                     double dVolatilityStrike,
                     OptionType optionType);

  /// Empty virtual destructor.
  virtual ~OptionVarianceSwap() { }


  /// @name Modifiers for option variance swap.
  //@{

  /**
      Sets user computed values at the valuation date.
     
      These values are only relevant if the valuation date is past the start 
      of the sampling period of the contract.  
      
      The current volatility defaults to zero.  Volatility calculations 
      must not include the drift term, and should use the return method 
      specified in the variance swap terms.

      The number of samples used should be less than the number of
      sampling returns specified in the variance swap terms.

      @param dCurrentVolatility current realized volatility
      @param nNbSamplesUsed number of samples used in the current volatility
   */
  void SetCurrentValues(double dCurrentVolatility, 
                        size_t nNbSamplesUsed);

  //@}

  /// @name Accessors for option variance swap.
  //@{

  /**
      The option payoff type.

      @return the option payoff type
   */
  OptionType GetOptionType() const
  {
    return m_optionType;
  }

  //@}

  virtual void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


protected:

  /// the option payoff type (put or cal)
  OptionType m_optionType;

}; // class OptionVarianceSwap


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_OPTIONVARIANCESWAP_H_
