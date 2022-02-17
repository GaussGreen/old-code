/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/sharedependentconversion.h
// Purpose:     Share dependent conversion for convertible bonds
// Author:      Ito 33
// Created:     2004/12/31
// RCS-ID:      $Id: sharedependentconversion.h,v 1.8 2006/06/15 11:26:39 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/sharedependentconversion.h
    @brief declaration of the share dependent conversion class.
 */

#ifndef _ITO33_FINANCE_BONDLIKE_SHAREDEPENDENTCONVERSION_H_
#define _ITO33_FINANCE_BONDLIKE_SHAREDEPENDENTCONVERSION_H_

#include "ito33/date.h"

#include "ito33/finance/bondlike/conversion.h"
#include "ito33/finance/bondlike/cocotype.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
   Class defining the share dependent conversion for convertible
   bonds with attached warrant
 */
class ITO33_DLLDECL ShareDependentConversion : public Conversion
{
public:

  /**
     Constructor. 

     @param startDate the beginning of the conversion period
     @param endDate the end of the conversion period
     @param dBaseRatio the base conversion ratio
     @param dIncrementalShareFactor the incremental share factor
   */
  ShareDependentConversion(
    Date startDate, Date endDate,
    double dBaseRatio, double dIncrementalShareFactor);

  ///@name modifiers ShareDependentConversion
  //@{

  /**
     The reset date, date at which the conversion stops 
     depending on the stock price.

     @param resetDate the date at which the conversion price becomes fixed
   */
  void SetResetDate(Date resetDate);

  /**
     The fixed strike.

     If not set, the fixed strike is equal to the principal 
     divided by the base conversion ratio, it accretes 
     for an OID or CashToZero bond

     @param dFixedStrike the fixed strike
   */
  void SetFixedStrike(double dFixedStrike); 
  
  /**
     The maximum convertion ratio allowed.

     @param dCapRatio the maximum conversion ratio allowed
   */
  void SetCapRatio(double dCapRatio);
 
  /**
     The current conversion ratio. 
     
     When the valuation date is after the reset date
     the current conversion ratio must be specified by the user.

     @param dCurrentRatio current conversion ratio
   */
  void SetCurrentRatio(double dCurrentRatio);

  //@}  // name methods for setting ShareDependentConversion

  ///@name accessors ShareDependentConversion
  //@{

  /**
     Gets the start of the conversion period.

     @return the start of the conversion period
   */
  Date GetStartDate() const 
  { 
    return m_StartDate; 
  }

  /**
     Gets the end of the conversion period.

     @return the end of the conversion period
   */
  Date GetEndDate() const 
  { 
    return m_EndDate; 
  }

  /**
     Gets the base conversion ratio.

     @return The base conversion ratio
   */
  double GetBaseRatio() const 
  { 
    return m_dBaseRatio; 
  }

  /**
     Gets the incremental share factor.

     @return the incremental share factor
   */
  double GetIncrementalShareFactor() const 
  { 
    return m_dIncrementalShareFactor; 
  }

  /**
     The fixed strike, value above which the conversion ratio increases.

     @return the incremental strike
   */
  double GetFixedStrike() const 
  { 
    return m_dFixedStrike; 
  }

  /**
     Gets the maximum conversion ratio.

     @return the maximum conversion ratio
   */
  double GetCapRatio() const 
  { 
    return m_dCapRatio; 
  }

  /**
     Indicates whether or not a reset date has been specified.

     @return true/false if a reset date has been set/not set
   */
  bool HasResetDate() const
  {
    return m_bHasResetDate;
  }

  /**
     Gets the reset date if it has been set.
 
     @return reset date
   */
  Date GetResetDate() const;
  
  /**
     The current conversion ratio.

     @return the current conversion ratio
   */
  double GetCurrentRatio() const
  {
    return m_dCurrentRatio;
  }


  //@}  // name methods for accessing ShareDependentConversion

  /// @name coco related
  //@{

  /**
     Trigger expressed as a percentage of the conversion price 
     (initial or accreted).
 
     @return trigger
   */
  double GetTrigger() const;

  /**
     The contingent conversion feature type.

     @return The CoCo type of this conversion period
   */
  CoCoType GetCoCoType() const;

  /**
     The annual rate at which the trigger level changes through time.

     @return Annual rate at which the trigger level changes through time
   */
  double GetChangeRate() const;

  /**
     The maximum or minimum trigger level expressed as percentage of
     the conversion price (initial or prevailing)

     @return The maximum or minimum trigger level.
   */
  double GetExtremeTrigger() const;

  /**
     Whether or not the last trigger condition has been met.
     
     Typically used for quarterly periods when the valuation date is within
     a quarter.  The code must be told if the trigger was reached at the
     end of the previous quarter.

     @return true if the conversion right is currently active, false otherwise
   */
  bool GetIsLastTriggerConditionMet() const ;

  /**
     The contingent conversion feature for this conversion period.
     The condition (stock price is very high) is checked at the end
     of each quarter in this period.

     When the valuation date occurs in the middle of a quarterly period,
     bIsLastTriggerConditionMet should be set to true if the stock was above the
     trigger at the end of the previous quarter.

     @param dTrigger Level to be reached by the underlying share price in 
                         order to trigger conversion. Expressed as a percentage 
                         of the initial or accreted conversion price.
     @param coCoType The type of the contingent conversion.
     @param dChangeRate Annual rate at which the trigger level 
                changes through time. It may be a positive or negative number.
     @param dExtremeTrigger Maximum or minimum trigger level. If the change 
                rate is positive, then dExtremeTrigger > dTrigger.
     @param bIsLastTriggerConditionMet Is the trigger condtion met.
   */
  void SetCoCo
       ( 
         double dTrigger, CoCoType coCoType,
         double dChangeRate, double dExtremeTrigger,
         bool bIsLastTriggerConditionMet = false
       );

  /**
     Checks if there is contingent conversion feature for this period, use it 
     to query CoCo related information for this period.

     @return true if there is contingent conversion, false otherwise 
   */
  bool HasCoCo() const { return m_coCoType < CoCoType_Max; }

  //@}  // coco related
  

  /**
     @internal 
     @brief Writes myself to parent tag.

     @param tagParent parent tag
     @return the tag created for the reset terms

     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// Start of the conversion period
  Date m_StartDate;

  /// End of the conversion period
  Date m_EndDate;

  /// The date at which the conversion ratio becomes fixed
  Date m_ResetDate;

  /// indicate if there is a reset date
  bool m_bHasResetDate;

  /// Base conversion ratio
  double m_dBaseRatio;

  /// The incremental share factor
  double m_dIncrementalShareFactor;

  /// The fixed asset value above which the conversion ratio increases
  double m_dFixedStrike;

  /// The maximum conversion ratio
  double m_dCapRatio;

  /// For CoCo, the trigger rate
  double m_dTriggerRate;
  
  /// For CoCo, the type of conversion allowed
  CoCoType m_coCoType;

  /// For CoCo, the rate at which the trigger changes over time
  double m_dChangeInTriggerRate;

  /// For CoCo, the maximum (pos. change) or minimum (neg. change) trigger rate
  double m_dExtremeTriggerRate;

  /// For CoCo, is the trigger active at the start date (eg. hit last quarter)
  bool m_bIsLastTriggerConditionMet;

  /// current conversion ratio
  double m_dCurrentRatio;

}; // class ShareDependentConversion


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_SHAREDEPENDENTCONVERSION_H_
