/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/cboption.h
// Purpose:     financial convertible bond option class
// Author:      Nabil
// Created:     2005/06/06
// RCS-ID:      $Id: cboption.h,v 1.13 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/cboption.h
    @brief declaration of the financial convertible bond option class 
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CBOPTION_H_
#define _ITO33_FINANCE_BONDLIKE_CBOPTION_H_

#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/finance/floatingrates.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/finance/bondlike/aswnotionalis.h"
#include "ito33/finance/bondlike/convertiblebond.h"

namespace ito33
{

namespace finance
{

/**
   A Convertible bond option.
 */
class ITO33_DLLDECL CBOption : public Derivative
{

public:

  /**
     Creates a convertible bond option.

     @param pConvertibleBond the convertible bond of the option.
     @param pASWFloatingRates the floating rates of the asset swap.
     @param maturityDate the maturity date of the option
   */
  CBOption(const shared_ptr<ConvertibleBond>& pConvertibleBond,
           const shared_ptr<FloatingRates>& pASWFloatingRates, 
           Date maturityDate);

  /// virtual dtor
  virtual ~CBOption() {}

  /** 
     The type of the notional of the asset swap.
   
     @param ASWNotionalType the type of the notional of the asset swap.
   */
  void SetASWNotionalIs(ASWNotionalIs ASWNotionalType);

  /** 
     The type of the notional of the asset swap.
   
     @return the type of the notional of the asset swap.
   */
  ASWNotionalIs GetASWNotionalIs() const
  {
    return m_ASWNotionalType;
  }

  /**
     @internal
     Sets the floating rates of the asset swap. 
     
     Added temporarily for testing stuff.

     @param pFloatingRates the floating rates

     @noexport
   */
  void SetASWFloatingRates(const shared_ptr<FloatingRates>& pFloatingRates)
  {
    m_pFloatingRates = pFloatingRates;
  }
  
  /**
     Gets the maturity date of the cb option.
   
     @return the maturity of the cb option.
   */
  Date GetMaturityDate() const;
  
  /**
     @internal

     @noexport
   */
  double GetCbNominal() const;
   
  /**
     @internal

     @noexport
   */
  double GetASWNotional() const;
   
  /**
     @internal

     @noexport
   */
  double GetASWRedemptionRate() const;

  /**
     @internal
     @brief Gets the balloon coupon if the notional of the asset swap is the
            issue price, 0 otherwise.

     @noexport
   */
  double GetBalloonCoupon() const;


  /**
     Gets the convertible bond.

     @return the convertible bond
   */
  const shared_ptr<ConvertibleBond>& GetConvertibleBond() const
  {
    return m_pConvertibleBond;
  }
  
  /**
     Gets the floating rates of the asset swap.

     @return the floating rates
   */
  const shared_ptr<FloatingRates>& GetASWFloatingRates() const
  {
    return m_pFloatingRates;
  }

  /**
     Computes the floating payments expressed in percentage of notional of 
     the asset swap.
     
     @return the floating payment stream of the asset swap
   */
  shared_ptr<CashFlowStream> ComputeASWFloatingPayments() const; 
  
  /**
     Computes the strike of the option, that is, the value of asset swap.

     @return the option strike
   */
  double ComputeStrike() const;

  /**
     Computes the fixed payments expressed in percentage of the nominal of
     the cb.
     
     @return the fixed payment stream of the asset swap
   */
  shared_ptr<CashFlowStream> ComputeASWFixedPayments() const;  

  // see base class
  virtual void PerturbFXRate(double dFXRateShift) const;

  virtual void SetSessionData(const shared_ptr<SessionData>& pSessionData);
  virtual void Validate() const;
  virtual void ValidateWith(const SessionData& sessionData) const;
  virtual void ValidateAll() const;

  virtual void Visit(DerivativeVisitor& visitor) const;
  
  virtual XML::Tag Dump(XML::Tag& tagParent) const;

protected:

  /**
     Gets shared ptr to derivative curve object 
     
     The derivative curve is the one used to discount the cash flow of the 
     derivative intrument.

     @return derivative curve
   */
  const shared_ptr<finance::YieldCurve>& GetDerivativeCurve() const;

  /// convertible bond
  shared_ptr<ConvertibleBond> m_pConvertibleBond;
  
  /// floating rates
  shared_ptr<FloatingRates> m_pFloatingRates;
  
  /// maturity date of the option
  Date m_maturityDate;

  /// cb option face value type
  ASWNotionalIs m_ASWNotionalType;

}; // class CBOption

// Inline functions

inline Date CBOption::GetMaturityDate() const
{
  return m_maturityDate;
}

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CBOPTION_H_
