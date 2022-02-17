/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/ConversionSchedule.h
// Purpose:     standard conversion provision for a bond
// Author:      ZHANG Yunzhi
// Created:     2004 may 3
// RCS-ID:      $Id: conversionschedule.h,v 1.27 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/ConversionSchedule.h
    @brief declaration of the ConversionSchedule class
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_
#define _ITO33_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_

#include "ito33/date.h"
#include "ito33/list.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/conversionperiod.h"
#include "ito33/finance/bondlike/conversion.h"

#include "ito33/dlldecl.h"


namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{


/**
   A ConversionSchedule class for a convertible bond.

   @iterator GetAll
 */
class ITO33_DLLDECL ConversionSchedule : public Conversion
{

public:

  /**
     Constructs an empty conversion schedule. 
   */
  ConversionSchedule() : Conversion(), m_bIsLastTriggerConditionMet(false)
  { }

  /// Type of the conversion data structure
  typedef std::list< shared_ptr<ConversionPeriod> > Elements;

  ///@name methods for initializing ConversionSchedule
  //@{

  /**
     Adds a conversion period to the conversion schedule.

     @param conversionPeriod Conversion period. 
   */
  void 
  AddConversionPeriod(const shared_ptr<ConversionPeriod>& conversionPeriod);

  /**
     Convenient shortcut for AddConversionPeriod.

     Allows to add a normal, i.e. non contingent, conversion more easily.

     @param startDate the beginning of the conversion period
     @param endDate the end of the conversion period
     @param dRatio the conversion ratio
   */
  void AddConversionPeriod(const Date& startDate,
                           const Date& endDate,
                           double dRatio)
  {
    AddConversionPeriod(shared_ptr<ConversionPeriod>(
        new ConversionPeriod(startDate, endDate, dRatio)));
  }

  /**
     Whether or not the CoCo trigger condition had been previously met 
     during the period containing the valuation date.
     
     Default value is false.

     @param bIsActive Is conversion active
   */
  void SetIsLastTriggerConditionMet(bool bIsActive)
  {
    m_bIsLastTriggerConditionMet = bIsActive;
  }

  //@}  // name methods for initializing ConversionSchedule

  ///@name methods for accessing ConversionSchedule
  //@{
  /**
     Whether or not the CoCo trigger condition had been previously met 
     during the period containing the valuation date.
     
     @return Whether or not CoCo is active
   */
  bool GetIsLastTriggerConditionMet() const
  {
    return m_bIsLastTriggerConditionMet;
  }

  //@}  // name methods for accessing ConversionSchedule


  /**
     Gets the list of the conversion periods.

     @return The list of the conversion periods.

     @noexport COM
   */
  const Elements& GetAll() const { return m_pConversions; }

  /**
     Validates the schedule with the given valuation date.
  	
     @param valuationDate the valuation date to validate against
   */
  void ValidateWith(Date valuationDate) const;

  /**
     @internal
     
     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// Conversion data structure
  Elements m_pConversions;

  /**
     Whether or not conversion is allowed during a CoCo quarterly period 
     that contains the valuation date.
   */
  bool m_bIsLastTriggerConditionMet;

}; // class ConversionSchedule


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_
