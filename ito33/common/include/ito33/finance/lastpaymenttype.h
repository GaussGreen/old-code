/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/lastpaymenttype.h
// Purpose:     enum for last payment type
// Created:     2006/02/13
// RCS-ID:      $Id: lastpaymenttype.h,v 1.2 2006/06/07 21:43:21 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/lastpaymenttype.h
    @brief Enumeration of the type of the last payment (short/long)
 */

#ifndef _ITO33_FINANCE_LASTPAYMENTTYPE_H_
#define _ITO33_FINANCE_LASTPAYMENTTYPE_H_

namespace ito33
{

namespace finance
{

/// type of the last payment: short or long
enum LastPaymentType
{
  LastPaymentType_Short,
  
  LastPaymentType_Long

  #ifndef __CPP2ANY__
  , 

  /// noexport
  LastPaymentType_Max

  #endif

}; // enum LastPaymentType

/**
    Checks if given last payment type is valid.

    @param lastPaymentType given value
    @return true if last payment type is valid, false if not 
    @noexport
 */
inline bool IsValid(LastPaymentType lastPaymentType)
{
  return lastPaymentType == LastPaymentType_Short ||
         lastPaymentType == LastPaymentType_Long;
}

/**
    @internal
    @brief throw exception because of invalid last payment type.
   
    @noexport
 */
void ThrowInvalidLastPaymentType();

/**
    @internal
    @brief Validates given last payment type and throw exception if it is
           invalid.

    @param lastPaymentType given value

    @noexport
 */
inline void Validate(LastPaymentType lastPaymentType)
{
  if ( !IsValid(lastPaymentType) )
    ThrowInvalidLastPaymentType();
}

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_LASTPAYMENTTYPE_H_
