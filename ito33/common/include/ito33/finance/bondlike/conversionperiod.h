/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/conversionperiod.h
// Purpose:     standard ConversionPeriod provision for a bond
// Author:      ZHANG Yunzhi
// Created:     2004 may 3
// RCS-ID:      $Id: conversionperiod.h,v 1.33 2006/06/15 11:08:41 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/conversionperiod.h
    @brief declaration of the ConversionPeriod class
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CONVERSIONPERIOD_H_
#define _ITO33_FINANCE_BONDLIKE_CONVERSIONPERIOD_H_

#include "ito33/date.h"

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
    ConversionPeriod class.
 */
class ITO33_DLLDECL ConversionPeriod
{
public:
  /**
     Constructs a Conversion Period without trigger.

     @param startDate Start date of the ConversionPeriod.
     @param endDate End date of the ConversionPeriod.
     @param dRatio Conversion ratio applicable during the ConversionPeriod.
   */
  ConversionPeriod(Date startDate, Date endDate, double dRatio);

  // copy constructor is ok

  /// @name general functions
  //@{

  /**
     Gets the conversion period start date.

     @return The conversion period start date.
   */
  Date GetStartDate() const { return m_startDate; }

  /**
     Gets the conversion period end date.

     @return The conversion period end date.
   */
  Date GetEndDate() const { return m_endDate; }

  /**
     Gets the conversion ratio applicable during the conversion period.

     @return The conversion ratio applicable during the conversion period.
   */
  double GetRatio() const { return m_dRatio; }

  /**
     Cash amount paid by the bondholder upon conversion. If 
     negative, the bondholder receives money upon conversion.

     @param dCash Cash amount paid or received by the bondholder upon 
                  conversion.
   */
  void SetCash(double dCash) { m_dCash = dCash; }

  /**
     Cash amount paid by the bondholder upon conversion. If 
     negative, the bondholder receives money upon conversion.

     @return Cash amount paid or received by the bondholder upon conversion.
   */
  double GetCash() const { return m_dCash; }

  ///@} // general functions

  /// @name coco related
  //@{

  /**
     Gets trigger rate expressed as a percentage of the conversion price 
     (initial or accreted).
 
     @return trigger
   */
  double GetTrigger() const;

  /**
     Gets the contingent conversion feature type.

     @return The CoCo type of this conversion period
   */
  CoCoType GetCoCoType() const;

  /**
     Gets the annual rate at which the trigger level changes through time.

     @return Annual rate at which the trigger level changes through time
   */
  double GetChangeRate() const;

  /**
     Gets the maximum or minimum trigger level expressed as percentage of
     the conversion price (initial or prevailing)

     @return The maximum or minimum trigger level.
   */
  double GetExtremeTrigger() const;

  /**
     Sets the contingent conversion feature for this conversion period.
     The condition (stock price is very high) is checked at the end
     of each quarter in this period.

     @param dTriggerRate Level to be reached by the underlying share price in 
                         order to trigger conversion. Expressed as a percentage 
                         of the initial or accreted conversion price.
     @param coCoType The type of the contingent conversion.
     @param dChangeRate Annual rate at which the trigger level 
                changes through time. It may be a positive or negative number.
     @param dExtremeTriggerRate Maximum or minimum trigger level. If the change 
                rate is positive, the dExtremTriggerRate > dTriggerRate.
   */
  void SetCoCo( double dTriggerRate, CoCoType coCoType,
                double dChangeRate, double dExtremeTriggerRate );

  /**
     Checks if there is contingent conversion feature for this period, use it 
     to query CoCo related information for this period.

     @return true if there is contingent conversion, false otherwise 
   */
  bool HasCoCo() const { return m_coCoType < CoCoType_Max; }


  //@}

  /**
     @internal
     
     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  Date m_startDate;

  Date m_endDate;

  double m_dRatio;

  double m_dCash;

  double m_dTriggerRate;
  
  CoCoType m_coCoType;

  double m_dChangeInTriggerRate;

  double m_dExtremeTriggerRate;
  
}; // class ConversionPeriod


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CONVERSIONPERIOD_H_
