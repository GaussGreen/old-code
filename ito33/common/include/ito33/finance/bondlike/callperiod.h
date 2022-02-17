/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/callperiod.h
// Purpose:     call period
// Author:      Wang (extracted from bondlike/callschedule.h by VZ)
// Created:     2004-09-08
// RCS-ID:      $Id: callperiod.h,v 1.30 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/bondlike/callperiod.h
   @brief declaration of the CallPeriod class for CB-like    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CALLPERIOD_H_
#define _ITO33_FINANCE_BONDLIKE_CALLPERIOD_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
   Defines a call period.

   @nocreate
 */
class ITO33_DLLDECL CallPeriod
{
public:
  
  /**
      Creates a call period object initialized with the strike.

      @param startDate Start date of the call
      @param endDate End date of the call
      @param dStrike Strike as fraction of the principal

      @noexport COM
   */
  static shared_ptr<CallPeriod>
  CreateWithStrike(Date startDate, Date endDate, double dStrike);

  /**
      Creates a call period object initialized with the yield.

      @param startDate Start date of the call
      @param endDate End date of the call
      @param dYield guaranted yield upon call

      @noexport COM
   */
  static shared_ptr<CallPeriod>
  CreateWithYield(Date startDate, Date endDate, double dYield);

  // copy constructor is ok

  /**
     The trigger of the call period, expressed as a percentage of the initial
     or accreted conversion price.

     @param dTriggerRate The trigger rate
   */
  void SetTrigger(double dTriggerRate);

  /**
     Gets the call period start date.

     @return The call period start date
   */
  Date GetStartDate() const { return m_startDate; }

  /**
     Gets the call period end date.

     @return the call period end date
   */
  Date GetEndDate() const { return m_endDate; }

  /**
     Whether the call value is defined by a guaranteed yield.

     @return true if a yield is guaranteed upon call, false if it
     consists of a fraction of principal
   */
  bool HasYield() const
  {
    return m_bHasYield;
  }

  /**
     Gets the call strike expressed as a percentage of the principal.

     @return The call strike expressed as a percentage of the principal.
  */
  double GetStrike() const;

  /**
     The trigger of the call period, expressed as a percentage of the initial
     or accreted conversion price.

     @return The trigger level expressed as a percentage of the initial
             or accreted conversion price
   */
  double GetTrigger() const { return m_dTriggerRate; }

  /**
     Get the guaranteed yield upon call.

     @return the guaranteed yield upon call
   */
  double GetGuaranteedYield() const;

  /**
     Checks if the bond may be called conditional on the underlying share 
     trading above a certain trigger.

     returns true if it is a soft call.
   */
  bool IsSoft() const { return m_dTriggerRate > 0; }

  /**
     @internal

     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// privatve ctor helps the creaters
  CallPeriod(Date startDate, Date endDate);

  /// Start date of the period
  Date m_startDate;

  /// End date of the period
  Date m_endDate;

  /// Strike of the period, in percentage
  double m_dStrike;

  /// Trigger of the period, in percentage
  double m_dTriggerRate;

  /// Guaranteed yield to call
  double m_dYield;

  /// has yield?
  bool m_bHasYield;

}; // class CallPeriod


} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_BONDLIKE_CALLPERIOD_H_
