/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/call.h
// Purpose:     base call class for CB-like instrument
// Author:      Wang
// Created:     2004/08/17
// RCS-ID:      $Id: call.h,v 1.46 2006/06/16 10:35:26 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/call.h
    @brief declaration of the base call class for CB-like instrument  
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CALL_H_
#define _ITO33_FINANCE_BONDLIKE_CALL_H_

#include <cstddef> // for size_t

#include "ito33/common.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"
#include "ito33/finance/bondlike/makewholetype.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{


/**
   Base call class for CB-like instrument. The following features are defined 
   in this class: CallNotice, MakeWhole etc.

   @noexport COM
   @nocreate
 */
class ITO33_DLLDECL Call
{

public:

  /// virtual dtor for base class
  virtual ~Call() { } 

  ///@name methods for initializing Calls
  //@{

  /**
     The keep Accrued flag.

     @param bKeepAccrued The keep accrued flag
   */
  void SetKeepAccrued(bool bKeepAccrued) { m_bKeepAccrued = bKeepAccrued; }
  
  /**
     The forfeit coupon flag.

     @param bForfeitCoupon The forfeit coupon flag
   */
  void SetForfeitCoupon(bool bForfeitCoupon)
  {
    m_bForfeitCoupon = bForfeitCoupon;
  }
  
  /**
     The call notice period.

     @param nNoticePeriod The call notice period as a number of days
   */
  void SetNoticePeriod(size_t nNoticePeriod);

  /**
     Specifies the premium make-whole by setting the premium paid to the 
     security holder upon the call of the security.

     @param dInitialPremium Make-whole premium value

     @method
   */
  void SetPremiumMakeWhole(double dInitialPremium);

  /**
     Specifies the coupon make-whole by setting the PVCouponMakeWhole flag
     upon the call of the instrument.

     @param bPVCouponMakeWhole true if the security holder gets the risk-free
                               PVd sum of the unpaid coupons during the soft 
                               call period
     
     @method
   */
  void SetCouponMakeWhole(bool bPVCouponMakeWhole)
  {
    m_makeWholeType = MakeWholeType_Coupon;
    m_bPVCouponMakeWhole = bPVCouponMakeWhole;
  }

  /**
     Specifies how the trigger (if any) will be checked.

     @param nTriggerPeriod number of consecutive days during which the call 
                            trigger condition must be met.
     @param nTriggerHistory number of consecutive days during which the call 
                            trigger condition has been met immediately prior 
                            to the valuation date.
   */
  void SetTriggerCheckPeriod(size_t nTriggerPeriod, size_t nTriggerHistory);
   

  /**
     Specifies whether the call triggers are expressed as a percentage of the fixed 
     or accreted conversion price.

     @param asPercentageOf The call triggers rates are a percentage of
   */
  void SetTriggerAsPercentageOf(TriggerAsPercentageOf asPercentageOf);

  //@}  // name methods for initializing Call

  ///@name methods for accessing Calls
  //@{

   
  /**
     Specifies whether the call triggers are expressed as a percentage of the fixed 
     or accreted call price.

     @return the call triggers as percentage of type.
   */
  TriggerAsPercentageOf GetTriggerAsPercentageOf() const 
  { 
    return m_triggerAsPercentageOf; 
  }

  /**
     The keep Accrued flag

     @return keep accrued flag
   */
  bool GetKeepAccrued() const { return m_bKeepAccrued; }
  
  /**
     The forfeit coupon flag

     @return forfeit coupon flag
   */
  bool GetForfeitCoupon() const { return m_bForfeitCoupon; }
  
  /**
     Checks if a call notice period has been specified.

     @return true if the call has a notice period, false otherwise.
   */
  bool HasNoticePeriod() const { return m_nNoticePeriod != 0; }
  
  /**
     The call notice period.

     @return the call notice period.
   */
  size_t GetNoticePeriod() const
  {
    return m_nNoticePeriod;
  }

  /**
     Checks if the security has a call make-whole feature.

     Note that HasMakeWhole feature is only relevant for soft call periods.

     @return true if the call has a make-whole provision.
   */
  virtual bool HasMakeWhole() const = 0;

  /**
     Gets the make-whole type (coupon or premium)

     @return the type of the make whole provision.
   */
  MakeWholeType GetMakeWholeType() const;

  /**
     Gets the premium if the security has a premium make-whole feature.

     @return the premium if the security has a premium make-whole feature.
   */
  double GetMakeWholePremium() const;

  /**
     Relevant only if the instrument has a coupon make-whole feature.

     When the instrument has a coupon make-whole, this function indicates 
     whether the CB-like instrument holder receives the sum of the risk-free
     present value coupons during the soft call period or the sum of absolute
     coupons.

     @return true if the CB-like instrument holder receives the sum of the 
             risk-free present value coupons during the soft call period, 
             false otherwise
   */
  bool IsPVCouponMakeWhole() const;

  /**
     Gets the number of consecutive days during which the call trigger 
     condition has to be met immediately before notifying the call.

     @return The number of days that the trigger will be checked
   */
  size_t GetTriggerPeriod() const { return m_nTriggerPeriod; }

  /**
     Gets the number of consecutive days during which the call trigger 
     condition has been met immediately prior to the valuation date.
   */
  size_t GetTriggerHistory() const { return m_nTriggerHistory; }

  //@}  // name methods for accessing Call


protected:

  /**
     Ctor sets members to appropriate values. 

     KeepAccrued is set to true and ForfeitCoupon is set to false, 
     no notice period, nor trigger period.
   */
  Call() : m_bKeepAccrued(true), m_bForfeitCoupon(false),
           m_nNoticePeriod(0),
           m_makeWholeType(MakeWholeType_Max),
           m_nTriggerPeriod(0), m_nTriggerHistory(0),
           m_triggerAsPercentageOf(TriggerAsPercentageOf_Principal)
  {
  }

  /**
     Writes myself to tag root of specific Call

     @param tagRoot CallSchedule etc.
   */
  void DumpMe(XML::Tag& tagRoot) const;

  /// keep accrued flag.
  bool m_bKeepAccrued;

  /// forfeit coupon flag.
  bool m_bForfeitCoupon;

  /// notice period expressed as a number of days.
  size_t m_nNoticePeriod;

  /// make-whole type.
  MakeWholeType m_makeWholeType;

  /// make-whole premium.
  double m_dMakeWholePremium;

  /**
     If the security holder gets the risk-free PVd sum of the unpaid coupons during
     the soft call period.
   */
  bool m_bPVCouponMakeWhole;
   
  /** 
     Number of consecutive days during which the call trigger condition must be
     met.
   */
  size_t m_nTriggerPeriod;
 
  /** 
     Number of consecutive days during which the call trigger condition has 
     been met immediately prior to the valuation date.
   */
  size_t m_nTriggerHistory;

  /**
     Of percentage of which value (fixed or accreted conversion price)
     is the trigger rate. Default to accreted conversion price.
   */
  TriggerAsPercentageOf m_triggerAsPercentageOf;

}; // class Call


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CALL_H_
