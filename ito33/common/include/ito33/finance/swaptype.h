/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/swaptype.h
// Purpose:     enum for different swap types
// Created:     2006/02/21
// RCS-ID:      $Id: swaptype.h,v 1.2 2006/04/10 10:55:10 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/swaptype.h
    @brief Enumeration of the type of swap (variance or volatility)
 */

#ifndef _ITO33_FINANCE_SWAPTYPE_H_
#define _ITO33_FINANCE_SWAPTYPE_H_

namespace ito33
{

namespace finance
{

/// type of the swap: variance or volatility
enum SwapType
{
  Swap_Variance,

  Swap_Volatility

  #ifndef __CPP2ANY__
  , 

  /// @noexport
  Swap_Max

  #endif

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_SWAPTYPE_H_
