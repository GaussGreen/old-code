/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/triggerincurrencyof.h
// Purpose:     enum to define the currency in which the trigger is expressed
// Author:      Nabil
// Created:     2004/09/09 
// RCS-ID:      $Id: triggerincurrencyof.h,v 1.6 2006/03/08 13:07:40 nabil Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/triggerincurrencyof.h
    @brief enum to define the currency in which the trigger is expressed  
 */

#ifndef _ITO33_FINANCE_BONDLIKE_TRIGGERINCURRENCYOF_H_
#define _ITO33_FINANCE_BONDLIKE_TRIGGERINCURRENCYOF_H_

namespace ito33
{

namespace finance
{


/// Trigger currency type.
enum TriggerInCurrencyOf
{
  /// Trigger on the underlying stock currency (trigger on conversion price).
  TriggerInCurrencyOf_Underlying,

  /// Trigger on the currency of derivative(Trigger on parity).
  TriggerInCurrencyOf_Derivative
  
  #ifndef __CPP2ANY__
  ,
  /// @noexport
  TriggerInCurrencyOf_Max
  #endif

}; // enum TriggerInCurrencyOf 


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_TRIGGERINCURRENCYOF_H_
