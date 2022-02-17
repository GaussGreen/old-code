/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/frequency.h
// Purpose:     Definition of the enumeration of the frequency
// Created:     2005/02/22 
// RCS-ID:      $Id: frequency.h,v 1.9 2006/05/26 13:34:31 nabil Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/frequency.h
   @brief Enumeration of the frequency of the payments (coupon, cds) or 
          the frequency of yield compounding
 */

#ifndef _ITO33_FINANCE_FREQUENCY_H_
#define _ITO33_FINANCE_FREQUENCY_H_

namespace ito33
{

namespace finance
{


/// Frequency of the payments or yield compounding
enum Frequency
{
#ifndef __CPP2ANY__
  /// @noexport
  Frequency_Undefined = 0,
#endif

  Frequency_Annual = 1,

  Frequency_SemiAnnual = 2,

  Frequency_Quarterly = 4,

  Frequency_BiMonthly = 6,

  Frequency_Monthly = 12
}; // enum Frequency

/**
   Checks if given frequency is valid.

   @param frequency given value
   @return true if frequency is valid, false if not 
   @noexport
 */
inline bool IsValid(Frequency frequency)
{
  return frequency == Frequency_Annual ||
         frequency == Frequency_SemiAnnual ||
         frequency == Frequency_Quarterly ||
         frequency == Frequency_BiMonthly ||
         frequency == Frequency_Monthly;
  
}

/**
   @internal
   @brief throw exception because of invalid frequency.
   
   @noexport
 */
void ThrowInvalidFrequency();

/**
   @internal
   @brief Validates given frequency and throw exception if it is invalid.

   @param frequency given value

   @noexport
 */
inline void Validate(Frequency frequency)
{
  if(!IsValid(frequency))
    ThrowInvalidFrequency();
}


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_FREQUENCY_H_
