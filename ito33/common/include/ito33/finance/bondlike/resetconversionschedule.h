/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/resetconversionschedule.h
// Purpose:     Reset schedule for a resettable convertible bond
// Author:      Yann and David
// Created:     2004/10/19
// RCS-ID:      $Id: resetconversionschedule.h,v 1.14 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/resetconversionschedule.h
    @brief declaration of the reset conversion schedule class.   
 */

#ifndef _ITO33_FINANCE_BONDLIKE_RESETCONVERSIONSCHEDULE_H_
#define _ITO33_FINANCE_BONDLIKE_RESETCONVERSIONSCHEDULE_H_

#include "ito33/list.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/bondlike/resetflooredby.h"
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

  class ITO33_DLLDECL ConversionPriceReset;

/**
   Class defining the reset schedule for resttables.

   @iterator GetAll
 */
class ITO33_DLLDECL ResetConversionSchedule : public Conversion
{
public:

  /**
     Constructor specifies the properties that apply to all reset dates.

     @param startDate the beginning of the conversion period
     @param endDate the end of the conversion period
     @param dInitialConvPrice the conversion price at the issue date.
            the conversion price at the issue date is in the currency
            of the underlying.
     @param dCurrentConvPrice the conversion price at the valuation date.
            the conversion price at the valuation date is in the currency
            of the underlying.
     @param flooredBy the method used to floor the conversion prices
   */
  ResetConversionSchedule(Date startDate, Date endDate, double dInitialConvPrice, 
                          double dCurrentConvPrice, ResetFlooredBy flooredBy);
    

  /// Type of the reset data structure
  typedef std::list< shared_ptr<ConversionPriceReset> > Elements;

  ///@name methods for initializing ResetConversionSchedule
  //@{
  
  /**
     Cash amount paid by the bondholder upon conversion. If 
     negative, the bondholder receives money upon conversion.

     @param dCash Cash amount paid or received by the bondholder upon 
                  conversion.
   */
  void SetCash(double dCash) { m_dCash = dCash; }
  
  /**
     Adds a reset date (with accompanying terms) to the reset schedule.

     @param pConversionPriceReset the reset terms at a given date. 
   */
  void AddConversionPriceReset
       (const shared_ptr<ConversionPriceReset>& pConversionPriceReset);

   /**
      @internal 

      @breif Sets the initial conversion price (at issue date).
      The initial conversion price is in the currency of
      the underlying.

      @param dInitialConvPrice initial conversion price (at issue date)

      @noexport
   */
  void SetInitialConversionPrice(double dInitialConvPrice);

  /**
     @internal

     @brief Sets the current conversion price (at valuation date).
     The current conversion price is in the currency of the
     underlying.

     @param dCurrentConvPrice current conversion price (at valuation date)

     @noexport
   */
  void SetCurrentConversionPrice(double dCurrentConvPrice);

  //@}  // name methods for initializing ResetConversionSchedule


  ///@name methods for accessing ResetConversionSchedule
  //@{

  /**
     Gets the start of the conversion period.

     @return the start of the conversion period
   */
  Date GetStartDate() const { return m_StartDate; }

  /**
     Gets the end of the conversion period.

     @return the end of the conversion period
   */
  Date GetEndDate() const { return m_EndDate; }

  /**
     Gets the initial conversion price (at issue date).

     @return initial conversion price (at issue date)
   */
  double GetInitialConversionPrice() const 
  { 
    return m_dInitialConvPrice; 
  }

  /**
     Gets the current conversion price (at valuation date).

     @return current conversion price (at valuation date)
   */
  double GetCurrentConversionPrice() const 
  { 
    return m_dCurrentConvPrice; 
  }

  /**
     Gets the method for computing the conversion price floor.

     @return the method for computing the conversion price floor
   */
  ResetFlooredBy GetResetFlooredBy() const 
  { 
    return m_ResetFlooredBy; 
  }

  /**
     Cash amount paid by the bondholder upon conversion. If 
     negative, the bondholder receives money upon conversion.

     @return Cash amount paid or received by the bondholder upon conversion.
   */
  double GetCash() const { return m_dCash; }


  //@}  // name methods for accessing ConversionSchedule

  /**
     Gets the list of the reset terms.

     @return the list of the reset terms.

     @noexport COM
   */
  const Elements& GetAll() const { return m_ResetTermList; }

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

  /// initial conversion price (at issue date)
  double m_dInitialConvPrice;

  /// current conversion price (at valuation date)
  double m_dCurrentConvPrice;

  /// How the floor is computed
  ResetFlooredBy m_ResetFlooredBy;

  /// reset term data structure
  Elements m_ResetTermList;

  /// cash amount paid by bondholder upon conversion
  double m_dCash;

}; // class ResetConversionSchedule


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_RESETCONVERSIONSCHEDULE_H_
