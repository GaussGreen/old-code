/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/aswnotionalis.h
// Purpose:     financial convertible bond option class
// Created:     2005/12/27
// RCS-ID:      $Id: aswnotionalis.h,v 1.3 2006/05/05 14:25:31 nabil Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/aswnotionalis.h
    @brief declaration of the enum for asset swap notional type
 */

#ifndef _ITO33_FINANCE_BONDLIKE_ASWNOTIONALIS_H_
#define _ITO33_FINANCE_BONDLIKE_ASWNOTIONALIS_H_

namespace ito33
{

namespace finance
{

/// Enum for asset swap notional type
enum ASWNotionalIs
{
  /**
     The notional of the asset swap is the issue price of the cb.

     This treats the case of:
       - Classic Par/Par asset swap
       - Premium redemption (bond issued at par (100%) but redeemed at a 
         premium (say 125%) )
       - The OID (if the issue price is used as notional)
    
     For the two previous cases, the net present value of a <b>non null</b> 
     balloon coupon is used in the computation of the strike of the option.
   */
  ASWNotionalIs_IssuePrice,

  /** 
     The notional of the asset swap is the put price.
   */
  ASWNotionalIs_PutPrice

  #ifndef __CPP2ANY__
  ,
  /// @noexport
  ASWNotionalIs_Max
  #endif

}; // enum ASWNotionalIs 

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_ASWNOTIONALIS_H_
