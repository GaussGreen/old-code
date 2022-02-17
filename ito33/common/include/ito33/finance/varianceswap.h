/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/varianceswap.h
// Purpose:     financial class for standard variance swaps
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswap.h,v 1.22 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/varianceswap.h
    @brief Declaration of the standard financial variance swap class.
 */

#ifndef _ITO33_FINANCE_VARIANCESWAP_H_
#define _ITO33_FINANCE_VARIANCESWAP_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/varianceswaplike.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL VarianceSwapTerms;
class ITO33_DLLDECL SessionData;

/// VarianceSwap represents the financial aspects of a standard variance swap.
class ITO33_DLLDECL VarianceSwap : public VarianceSwapLike
{

public:
  /**
      Creates a variance swap by taking the terms and the volatility strike.
      
      @param pTerms The terms of the swap
      @param dVolatilityStrike The volatility strike of the swap
   */
  VarianceSwap(const shared_ptr<VarianceSwapTerms>& pTerms,
               double dVolatilityStrike);

  /// Empty virtual destructor.
  virtual ~VarianceSwap() { }


  /// @name Modifiers for variance swap.
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

  virtual void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;

}; // class VarianceSwap


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_VARIANCESWAP_H_
