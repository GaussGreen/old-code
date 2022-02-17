/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/triggeraspercentageof.h
// Purpose:     enum to indicate the type of the reference value used to 
//              express the trigger
// Author:      Nabil
// Created:     2004/09/01 
// RCS-ID:      $Id: triggeraspercentageof.h,v 1.13 2006/03/15 14:45:30 nabil Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/triggeraspercentageof.h
    @brief enum to indicate the type of the reference value used to express 
           the trigger
 */

#ifndef _ITO33_FINANCE_BONDLIKE_TRIGGERASPERCENTAGEOF_H_
#define _ITO33_FINANCE_BONDLIKE_TRIGGERASPERCENTAGEOF_H_

namespace ito33
{

namespace finance
{


/// Type of the reference value used to express the trigger.
enum TriggerAsPercentageOf
{
  /// trigger as a percentage of the accreted conversion price
  TriggerAsPercentageOf_Principal,

  /**
     trigger as a percentage of the fixed conversion price: 
     IssuePrice/ConvertionRatio.
   */
  TriggerAsPercentageOf_IssuePrice,
  
  /** 
     Trigger as a percentage of the accreted conversion price
     including Accrued Interest.
   */
  TriggerAsPercentageOf_Claim

  #ifndef __CPP2ANY__
  ,
  /// @noexport
  TriggerAsPercentageOf_Max
  #endif

}; // enum TriggerAsPercentageOf 


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_TRIGGERASPERCENTAGEOF_H_
