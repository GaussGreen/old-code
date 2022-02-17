/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/cb_base.h
// Purpose:     base class for convertible bond like classes
// Author:      ITO33
// Created:     2004/10/13
// RCS-ID:      $Id: cb_base.h,v 1.23 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/cb_base.h
    @brief declaration of  base class for convertible bond like classes
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CB_BASE_H_
#define _ITO33_FINANCE_BONDLIKE_CB_BASE_H_

#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/bondterms.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL CallSchedule;
class ITO33_DLLDECL PutSchedule;
class ITO33_DLLDECL Bond;

/**
    Base class for convertible bond like classes.
    
    @nocreate
 */
class ITO33_DLLDECL CBBase : public ConvertibleLike
{
public:


  /// virtual dtor for base class
  virtual ~CBBase() { }

  /**
     The call feature of the bond.

     @param pCallSchedule Shared pointer to the CallSchedule.
   */
  virtual void SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule);
  
  /**
     The put feature of the bond.

     @param pPutSchedule Shared pointer to the PutSchedule.
   */
  void SetPutSchedule(const shared_ptr<PutSchedule>& pPutSchedule);

  /**
      @name Methods for accessing convertible-like properties.
   */
  //@{

  /**
     Gets the terms of the bond.
    
     @return The bond terms
   */
  const shared_ptr<BondTerms>& GetBondTerms() const
  {
    return m_pBondTerms; 
  }

  /**
     The call feature.
    
     @return Call data.
   */
  const shared_ptr<CallSchedule>& GetCallSchedule() const 
  { 
    return m_pCallSchedule; 
  }

  /**
     The put feature.
    
     @return The put data.
   */
  const shared_ptr<PutSchedule>& GetPutSchedule() const 
  { 
    return m_pPutSchedule; 
  }

  //@}

  /**
     Computes the yield-to-maturity, throws exception if it is impossible
     to calculate the yield.

     @param dPrice The price used to compute the yield-to-maturity
     @param cmpFrequency The yield compound frequency
     @param dcc The yield day count convention

     @return The yield-to-maturity
   */
  double 
  ComputeYieldToMaturity
  (double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const;
  
  /**
     Computes the yield-to-put, throws exception if it is impossible to 
     calculate the yield.

     @param dPrice The price used to compute the yield-to-put
     @param cmpFrequency The yield compound frequency
     @param dcc The yield day count convention

     @return The yield-to-put
   */
  double 
  ComputeYieldToPut
  (double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const;
  
  /**
     Computes the put price at the given date.

     Throw an exception if there are no put at the given date.

     @param putDate the date at which the put price needs to be computed.

     @return the put price at the given date if there are one at this date.
   */
  double ComputePutPrice(Date putDate) const;
  
  /**
     Computes the put price at the given date. If the put price contains
     a coupon, that is the put date is also a coupon date and the forfeit
     coupon flag is false, the function returns also the amount of coupon 
     value.

     Throw an exception if there is no put at the given date.

     @param putDate the date at which the put price needs to be computed.
     @param dCouponAmount (output) the coupon value in put price. So 0 if
            the put date is not a coupon date or forfeit coupon is true.

     @return the put price at the given date if there are one at this date.

     @noexport
   */
  double ComputePutPrice(Date putDate, double& dCouponAmount) const;

  /**
     Computes the call price at the given date.

     Throw an exception if there are no call at the given date.

     @param callDate the date at which the call price needs to be computed.

     @return the call price at the given date if there are one at this date.
   */
  double ComputeCallPrice(Date callDate) const;

  /**
     Gets the accrued interest value at valuation date

     @return accrued interest value at valuation date
   */
  virtual double GetAccruedInterestValue() const;

  /**
     Gets the straight bond corresponding to this convertible.

     @return the straight bond of this convertible.
   */
  shared_ptr<Bond> GetStraightBond() const;

  virtual void Validate() const;

protected:
  
  /**
     ctor.

     @param pBondTerms The main characteristics of the bond.
   */
  CBBase(const shared_ptr<BondTerms>& pBondTerms);

  /**
     Writes myself to tag parent which can be ConvertibleBond
     or Reset or some other derived class tag

     @param tagParent tag of the derived class
   */
  void DumpMe(XML::Tag& tagParent) const;

  /// bond like terms
  shared_ptr<BondTerms> m_pBondTerms;

  /// call data.
  shared_ptr<CallSchedule> m_pCallSchedule;

  /// put data.
  shared_ptr<PutSchedule> m_pPutSchedule;

}; // class CBBase


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CB_BASE_H_
