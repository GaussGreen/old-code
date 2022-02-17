/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/swappayofftype.h
// Purpose:     enum for different swap payoff types
// Created:     2006/07/12
// RCS-ID:      $Id: swappayofftype.h,v 1.2 2006/08/03 20:05:12 dave Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/swappayofftype.h
    @brief Enumeration of the variance swap payoff type (standard, gamma, 
           conditional, etc.)
 */

#ifndef _ITO33_FINANCE_SWAP_PAYOFF_TYPE_H_
#define _ITO33_FINANCE_SWAP_PAYOFF_TYPE_H_

namespace ito33
{

namespace finance
{

/// type of the swap payoff: standard, gamma, conditional
enum SwapPayoffType
{
  SwapPayoff_Standard,

  SwapPayoff_Gamma,

  SwapPayoff_Conditional,

  SwapPayoff_Call,

  SwapPayoff_Put

  #ifndef __CPP2ANY__
  , 

  /// @noexport
  SwapPayoff_Max

  #endif

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_SWAP_PAYOFF_TYPE_H_
