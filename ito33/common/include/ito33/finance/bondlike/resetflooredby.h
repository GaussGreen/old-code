/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/resetflooredby.h
// Purpose:     enum to indicate how the conversion price is floored at reset
//              dates
// Author:      David
// Created:     2004/10/05 
// RCS-ID:      $Id: resetflooredby.h,v 1.3 2005/04/15 12:26:45 zhang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/resetflooredby.h
    @brief enum to indicate how the conversion price is floored at reset dates
 */

#ifndef _ITO33_FINANCE_BONDLIKE_RESETFLOOREDBY_H_
#define _ITO33_FINANCE_BONDLIKE_RESETFLOOREDBY_H_

namespace ito33
{

namespace finance
{


/// How the conversion price is floored
enum ResetFlooredBy
{
  /// floored by the prevailing conversion price
  ResetFlooredBy_PrevailingConversionPrice,

  /// floored by the initial conversion price at issue date
  ResetFlooredBy_InitialConversionPrice
  
  #ifndef __CPP2ANY__
  ,
  /// @noexport
  ResetFlooredBy_Max
  #endif

}; // enum ResetFlooredBy 


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_RESETFLOOREDBY_H_
