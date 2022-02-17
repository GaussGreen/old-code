/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/bond.h
// Purpose:     financial bond class
// Author:      ZHANG Yunzhi
// Created:     2004 may 3
// RCS-ID:      $Id: bond.h,v 1.54 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/bond.h
    @brief declaration of the financial bond class
 */

#ifndef _ITO33_FINANCE_BONDLIKE_BOND_H_
#define _ITO33_FINANCE_BONDLIKE_BOND_H_

#include "ito33/finance/derivative.h"
#include "ito33/finance/frequency.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/** 
    @name Forward Declaration
*/
//@{
class ITO33_DLLDECL BondTerms;
class ITO33_DLLDECL CallSchedule;
class ITO33_DLLDECL CashFlowStream;
class ITO33_DLLDECL PutSchedule;
//@}

/**
   A bond with calls and puts.
 */
class ITO33_DLLDECL Bond : public Derivative
{
public:
  
  /**
     Constructor using BondTerms.

     @param pBondTerms the main characteristics of the bond.
   */
  Bond(const shared_ptr<BondTerms>& pBondTerms);

  /**
     Empty destructor.
   */
  virtual ~Bond() { }

  /**
      @name Methods for initializing a Bond.
   */
  //@{

  /**
     The call feature of the bond.

     @param pCallSchedule shared pointer to the CallSchedule.
   */
  virtual void SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule);
  
  /**
     The put feature of the bond.

     @param pPutSchedule shared pointer to the PutSchedule.
   */
  void SetPutSchedule(const shared_ptr<PutSchedule>& pPutSchedule);

  //@}

  /**
      @name Methods for accessing to Bond.
   */
  //@{
  
  /**
     Gets the terms of the bond.
    
     @return the bond terms
   */
  const shared_ptr<BondTerms>& GetBondTerms() const
  {
    return m_pBondTerms; 
  }

  /**
     The call feature.
    
     @return call data.
   */
  const shared_ptr<CallSchedule>& GetCallSchedule() const 
  { 
    return m_pCallSchedule; 
  }

  /**
     The put feature.
    
     @return the put data.
   */
  const shared_ptr<PutSchedule>& GetPutSchedule() const 
  { 
    return m_pPutSchedule; 
  }

  //@}  

  /**
     Computes the yield-to-maturity, throws exception if it is impossible
     to calculate the yield.

     The frequency/day count convention might be different than 
     the frequency/day count convention of the coupon/yield of this bond,
     to facilitate comparison with another bond. 

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

     The frequency/day count convention might be different than 
     the frequency/day count convention of the coupon/yield of this bond,
     to facilitate comparison with another bond.

     @param dPrice The price used to compute the yield-to-put
     @param cmpFrequency The yield compound frequency
     @param dcc The yield day count convention

     @return The yield-to-put
   */
  double 
  ComputeYieldToPut
  (double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const;

  /**
     Computes the put price at the given date with consideration of KeepAccrued
     and ForfeitCoupon flags.
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
     Computes the call price at the given date with consideration of 
     KeepAccrued and ForfeitCoupon flags.
     Throw an exception if there are no call at the given date.

     @param callDate the date at which the call price needs to be computed.

     @return the call price at the given date if there are one at this date.
   */
  double ComputeCallPrice(Date callDate) const;
   
  /**
     Gets the accrued interest value at valuation date.

     @return accrued interest value at valuation date
   */
  double GetAccruedInterestValue() const;

  // implements base class pure virtuals
  virtual Date GetMaturityDate() const;
  virtual void Validate() const;

  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual void Visit(DerivativeModifyingVisitor& visitor);
  virtual XML::Tag Dump(XML::Tag& tagParent) const;

  

private:

  /// Bond terms
  shared_ptr<BondTerms> m_pBondTerms;

  /// call data.
  shared_ptr<CallSchedule> m_pCallSchedule;

  /// put data.
  shared_ptr<PutSchedule> m_pPutSchedule;

}; // class Bond 


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_BOND_H_
