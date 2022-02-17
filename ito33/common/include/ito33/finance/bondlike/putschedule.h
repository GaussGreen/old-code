/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/putschedule.h
// Purpose:     standard put provision for a bond
// Author:      ZHANG Yunzhi
// Created:     2004 may 3
// RCS-ID:      $Id: putschedule.h,v 1.29 2006/07/27 18:47:33 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/putschedule.h
    @brief declaration of the put schedule class.   
 */

#ifndef _ITO33_FINANCE_BONDLIKE_PUTSCHEDULE_H_
#define _ITO33_FINANCE_BONDLIKE_PUTSCHEDULE_H_

#include "ito33/beforestd.h"
#include <map>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/date.h"

#include "ito33/dlldecl.h"

#include "ito33/finance/bondlike/putperiod.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
   Simple data structure to store the put data, which is either a strike
   or a yield to maturity.  At most one value should be non-zero.

   @noexport
 */
class PutData
{
public:
  /**
     Default constructor only required by std::map<Date, PutData>.
   */
  PutData() : dYield(0.), dStrike(-1.), bYield(false) {};

  /**
     Creates PutData object. The input must have already been validated.
   */
  PutData(double dValue, bool bIsYield)
  {
    bYield = bIsYield;
    if(bYield)
    {
      dYield = dValue;
      dStrike = -1.;
    }
    else
    {
      dStrike = dValue;
      dYield = 0.;
    }
  }

  /// yield
  double dYield;

  /// strike
  double dStrike;

  // boolean flag to indicate is strike or yield to maturity
  bool bYield;
};

/**
   A put schedule class.
   
   @iterator GetPutPeriods
 */
class ITO33_DLLDECL PutSchedule
{
public:

  /**
     Constructs an empty put schedule. 

     KeepAccrued is set to true and ForfeitCoupon is set to false
   */
  PutSchedule() : m_bKeepAccrued(true), m_bForfeitCoupon(false) {}

  /// type of the put data structure
  typedef std::map<Date, PutData > Elements;

  /// type of the PutPeriod structure
  typedef std::list<PutPeriod> CollectionElements;

  ///@name methods for initializing PutSchedule
  //@{
  
  /**
     The keep accrued flag.

     @param bKeepAccrued keep accrued flag.
   */
  void SetKeepAccrued(bool bKeepAccrued) 
  { 
    m_bKeepAccrued = bKeepAccrued; 
  }
  
  /**
     The forfeit coupon flag.

     @param bForfeitCoupon forfeit coupon flag.
   */
  void SetForfeitCoupon(bool bForfeitCoupon)
  {
    m_bForfeitCoupon = bForfeitCoupon;
  }
  
  /**
     Adds a put defined by its date and strike.

     @param putDate date of the put
     @param dPutStrike strike of put expressed as a percentage of principal.
   */
  void AddPutWithStrike(Date putDate, double dPutStrike);

  /**
     Adds a put defined by its date and yield.

     @param putDate date of the put
     @param dYield guaranteed yield upon put.
   */
  void AddPutWithYield(Date putDate, double dYield);

  //@}  // name methods for initializing PutSchedule


  ///@name methods for accessing PutSchedule
  //@{
  
  /**
     The keep Accrued flag.

     @return keep accrued flag.
   */
  bool GetKeepAccrued() const { return m_bKeepAccrued; }
  
  /**
     The forfeit coupon flag.

     @return forfeit coupon flag.
   */
  bool GetForfeitCoupon() const { return m_bForfeitCoupon; }
  
  /**
     Gets all the put periods.

     @noexport COM
   */
  const CollectionElements& GetPutPeriods() const;

  //@}  // name methods for accessing PutSchedule

  /**
     @internal
     @brief Checks if there is a put at the given date

     @param date The date at which we are looking for a put
   */
  bool HasPutAt(Date date) const
  {
    return m_puts.find(date) != m_puts.end();
  }

  /// Validates the put schedule and throws exception if validation fails
  void Validate() const;

  /**
     @internal
     @brief Indicates if the putschedule contains a put with guaranteed yield.

     @return true if the putschedule contains a put with guaranteed yield, 
             false otherwise

     @noexport
   */
  bool HasYield() const;

  /**
     @internal
     @brief Gets the list of puts.

     @return the list of puts.

     @noexport
   */
  const Elements& GetAll() const { return m_puts; }

  /**
     @internal
     
     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;

private:

  /// initialize m_pPutPeriods
  void CreateCollectionElements() const;

private:

  /// keep accrued flag
  bool m_bKeepAccrued;

  /// forfeit coupon flag
  bool m_bForfeitCoupon;

  /// PutSchedule data structure
  Elements m_puts;
  
  /// container of PutPeriods
  CollectionElements m_pPutPeriods;

}; // class PutSchedule


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_PUTSCHEDULE_H_
