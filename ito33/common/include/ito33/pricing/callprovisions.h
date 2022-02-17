/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/callprovisions.h
// Purpose:     Abstract base class for call provisions
// Author:      Nabil
// Created:     2004/03/11
// RCS-ID:      $Id: callprovisions.h,v 1.67 2006/07/03 16:08:15 nabil Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/callprovisions.h
   @brief base class at pricing level for call.
 */

#ifndef _ITO33_PRICING_CALLPROVISIONS_H_
#define _ITO33_PRICING_CALLPROVISIONS_H_

#include "ito33/vector.h"

#include "ito33/finance/bondlike/makewholetype.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Call;
}

namespace pricing
{

class CBLikeParams;


/**
   base class at pricing level for call.

   IMPORTANT:
   Instance of this class can't be shared, at least we don't recommend it.
   The reason is that some state variables such as Params, index etc. are
   members of this class.
 */
class CallProvisions
{

public:

  /**
     Default ctor sets approriate default value for members.
   */
  CallProvisions()
    : m_nNbCalls(0), m_bKeepAccrued(true), m_bForfeitCoupon(false),
      m_dNoticePeriod(0.), m_bHasNoticePeriod(false),
      m_iMakeWholeType(finance::MakeWholeType_Max),
      m_bDiscountCouponForMakeWhole(false),
      m_dMakeWholePremium(0),
      m_nTriggerPeriod(0), m_nTriggerHistory(0),
      m_triggerAsPercentageOf(finance::TriggerAsPercentageOf_Max)
  { }
  
  /// virtual dtor for base class
  virtual ~CallProvisions() { }
   
  /**
     Computes the call strike for the given spots without taking into account
     of coupons at a given time.
 
     @param dTime the time at which the strikes will be evaluated.
     @param pdS given spot array
     @param nNbS @sa number given spots
     @param pdNewSharePrices new share prices (or spot array)
     @param pdValues (output) call strikes
     @param bPlus at which side the call strikes will be evaluated at dTime.
   */
  virtual void GetCallStrikesWithoutCoupon
  (double dTime, const double* pdS, size_t nNbS, 
   const double* pdNewSharePrices, double* pdValues, 
   bool bPlus = true) const = 0;

  /**
     Computes the call strike for the given spots.
 
     @param pdS given spot array
     @param nNbS @sa number given spots
     @param pdNewSharePrices new share prices (or spot array) 
     @param pdValues (output) call strikes
   */
  virtual void GetCallStrikes
  (const double* pdS, size_t nNbS, 
   const double* pdNewSharePrices, double* pdValues) const = 0;

  /**
     Evaluates call constraint values, if any, for given spots. If there is 
     call notice, pdValues is also an input, otherwise, use the call strike.

     The array pdValues shouldn't be used when the function returns false.

     @TODO : should the following functions return bool or void? That is,
             should we check index inside or outside the function

     @param pdS given spot array
     @param nNbS @sa number given spots
     @param pdValues (output) call constraint values                  
     @param nIdxStartConversion (output) index of the point where the 
                               instrument is called to be converted.               
     @param pdNewSharePrices the new share prices
     @param bNoticed if there is notcie period
     @return true when we have call clause here. otherwise false,
             and does nothing for pdValues
   */
  bool GetCallConstraintValues
       (const double* pdS, size_t nNbS, double* pdValues, 
        size_t& nIdxStartConversion, 
        const double* pdNewSharePrices, 
        bool bNoticed = false) const;
  
  bool IsKeepAccrued() const { return m_bKeepAccrued; }
  
  bool IsForfeitCoupon() const { return m_bForfeitCoupon; }
  
  size_t GetNbCalls() const { return m_nNbCalls; }

  /**
     Gets start time for given call period
  
     @param nIdx index of given call period
   */
  double GetStartTime(size_t nIdx) const
  {
    return m_pdStartTimes[nIdx];
  }
  
  /**
    Get end time for given call period
  
    @param nIdx index of given call period
    */
  double GetEndTime(size_t nIdx) const
  {
    return m_pdEndTimes[nIdx];
  }  

  /**
     Helper function to check if the call is continous, so not monodate.
   */
  bool IsContinuous(size_t nIdx) const 
  {
    return m_pdStartTimes[nIdx] != m_pdEndTimes[nIdx];
  }

  /**
     Gets trigger for given call period.
  
     @param nIdx index of given call period
   */
  double GetTriggerRate(size_t nIdx) const
  {
    return m_pdTriggerRates[nIdx];
  }  

  /**
     Gets the yield to call the for given call period.
  
     @param nIdx index of given call period
   */
  double GetYieldToCall(size_t nIdx) const
  {
    return m_pdYieldToCall[nIdx];
  }  

  /**
     Sets pointer to related params object

     @param pParams pointer to params
   */
  void SetParams(const CBLikeParams* pParams)
  {
    m_pParams = pParams;
  }

  /**
     Checks if the calls has notice period

     @return wheter the bond Has notice period
   */
  bool HasNoticePeriod() const
  {
    return m_bHasNoticePeriod;
  }

  /**
     Gets the notice period as a fraction of year.

     @return Notice period in a fraction of year;
   */
  double GetNoticePeriod() const
  {
    return m_dNoticePeriod;
  }

  /**
     Gets the Index of continuous call at given time and the given side,
     returns INVALIDINDEX if there is no continuous call.

     @param dTime given time
     @param bPlus at which side we are looking for index
     @return index of continous call at time dTime
   */
  size_t GetIndexContinousCallAt(double dTime, bool bPlus);

  /**
     Gets make-whole premium value

     @return make-whole premium value
   */
  double GetMakeWholePremium() const
  {
    return m_dMakeWholePremium;
  }

  /**
     Discount coupon for Gets make-whole or not

     @return true if discount coupon for Gets make-whole, otherwise false
   */
  bool ShouldDiscountCouponForMakeWhole() const
  {
    return m_bDiscountCouponForMakeWhole;
  }

  /**
     Gets make-whole type

     @return make-whole type
   */
  finance::MakeWholeType GetMakeWholeType() const
  {
    return m_iMakeWholeType;
  }

  /**
     Whether make-whole is paid.

     @return whether make-whole is paid.
   */
  bool HasMakeWhole() const
  {
    return m_iMakeWholeType != finance::MakeWholeType_Max;
  }

  /**
     Gets the number of days that the trigger will be checked

     @return Gets the number of days that the trigger will be checked
  */
  size_t GetTriggerPeriod() const 
  { 
    return m_nTriggerPeriod; 
  }

  /**
     Gets the number of days that the trigger will be checked on the valuation 
     date.
   */
  size_t GetTriggerHistory() const 
  { 
    return m_nTriggerHistory; 
  }

  /// Removes the call periods.  Used when pricing an N-day call window
  virtual void Clear() = 0;

  /**
     Gets the trigger level at current time by using the conversion price.
   */
  double GetTriggerLevel() const;

  /**
     @name helper functions for "virtual"calls.
     
     In some case ( n days contingent call for example), even the calls should
     be turned off, we still need to add the call trigger to the space mesh to
     help interpolation for path dependent code. 

     The current implementation is to keep all call datas, but use a flag to 
     turn off the call constraints and to add the call trigger to the space 
     mesh. The call trigger computation (and claim etc) is duplicated, but it 
     hopefully won't take much time compared to system solving.
   */
  //@{

  /**
     Deactivate all softcall but not monodate.
     This is done to ensure that we can combine contingent
     call with hardcall.
   */
  void DeactivateSoftCall(); 

  /**
     Check if the call is active.
   */
  bool IsActive() const;

  /** 
     Turn off call notice.  Used for the intermediate paths (not paths 0 or the top
     path) in N-day call window pricing, for example. Prevents the inactive paths
     from solving mini-cb problems when constraints are updated.
   */
  void DeactivateNoticePeriod(); 

  /**
     Deactivates the trigger period, the trigger period is
     simply set to zero.
   */
  void DeactivateTriggerPeriod();
 
  /**
     Deactivates all calls, this includes notice period and trigger period.
   */
  void DeactivateAllCalls();

  /**
     Gets the trigger as percentage of.
    
     @param the trigger as percentage of
   */
  finance::TriggerAsPercentageOf GetTriggerAsPercentageOf() const
  {
    return m_triggerAsPercentageOf; 
  }

  //@}

  /**
     Determine if there is the active call periods with a trigger period which 
     requires a path dependent pricing. 

     @return true if path depdedent call, false otherwise
   */
  bool HasPathDepCall(double dValuationTime, double dStoppingTime) const;


protected:
  
  /**
     Helper fucntion to get the common call data from a financial Call objet.

     @param call a financial call objet

     @todo makewhole property should also be included
   */
  void GetCommonCallData(const finance::Call& call);

  /// keep accrued flag
  bool m_bKeepAccrued;

  /// forfeit coupon
  bool m_bForfeitCoupon;

  /// has notice period?
  bool m_bHasNoticePeriod;

  /// notice period
  double m_dNoticePeriod;

  /// makewhole type
  finance::MakeWholeType m_iMakeWholeType;

  /// make-whole premium
  double m_dMakeWholePremium;

  /// make-whole coupon discounted property
  bool m_bDiscountCouponForMakeWhole;

  /// number of call windows
  size_t m_nNbCalls;
  
  /// start time array
  std::vector<double> m_pdStartTimes;

  /// end time array
  std::vector<double> m_pdEndTimes;

  /// trigger rate, as a percentage of conversion price
  std::vector<double> m_pdTriggerRates;

  /// yield to call array
  std::vector<double> m_pdYieldToCall;

  // we save the pointer here, otherwise almost all
  // member function should take params as argument
  /// pointer to params 
  const CBLikeParams* m_pParams;

  /** 
     number of consecutive days during which the call trigger condition must be
     met.
   */
  size_t m_nTriggerPeriod;

  /** 
     number of consecutive days during which the call 
     trigger condition has been met immediately prior 
     to the valuation date.
   */
  size_t m_nTriggerHistory;

  /// boolean indicating if we should take into account of the call constraint
  std::vector<bool> m_pbIsActives;
  
  /**
     Of percentage of which value (fixed or accreted call price)
     is the trigger rate. Default to accreted conversion price.
   */
  finance::TriggerAsPercentageOf m_triggerAsPercentageOf;

}; // class CallProvisions


} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_CALLPROVISIONS_H_
