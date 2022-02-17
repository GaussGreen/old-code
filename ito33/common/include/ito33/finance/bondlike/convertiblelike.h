/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/convertiblelike.h
// Purpose:     convertible-like bond class
// Author:      ITO33
// Created:     2004/10/13
// RCS-ID:      $Id: convertiblelike.h,v 1.30 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/convertiblelike.h
    @brief declaration of the financial convertible-like bond class.
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CONVERTIBLELIKE_H_
#define _ITO33_FINANCE_BONDLIKE_CONVERTIBLELIKE_H_

#include "ito33/finance/derivative.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"
#include "ito33/finance/bondlike/triggerincurrencyof.h"

namespace ito33
{

namespace finance
{

/// @name Forward declaration
//@{

class ITO33_DLLDECL BondLikeTerms;

//@}

/**
   A Convertible-like security.

   @nocreate
 */
class ITO33_DLLDECL ConvertibleLike : public Derivative
{
public:

  /// virtual dtor for base class
  virtual ~ConvertibleLike() { }
  
  /**
      @name Methods for setting convertible-like properties.
   */
  //@{

  /**
     Defines an exchangeable bond with the ExchangeableUponDefault property.

     The user must define the default of bond's issuer.

     @param bExchangeableUponDefault Specifies whether the bondholder is 
                entitled to convert upon default of the issuer.
   */
  void SetExchangeable(bool bExchangeableUponDefault);

  /**
     Specifies whether the conversion triggers are expressed as a percentage
     of the fixed or accreted conversion price.

     @param asPercentageOf The conversion triggers rates are a percentage of
   */
  void SetConversionTriggerAsPercentageOf(TriggerAsPercentageOf asPercentageOf);
  
  /**
     The trigger currency type (for call and conversion) in case of a
     cross currency: whether the trigger is in the security currency
     or in the underlying currency.

     @param inCurrencyOf The triggers are in currency of.
   */
  void SetTriggerInCurrencyOf(TriggerInCurrencyOf inCurrencyOf);
  
  /**
     The Fixed FX Rate (for cross currency). The FX rate is the number of 
     base currency (i.e the derivative currency) needed to purchase one unit 
     of the foreign currency (i.e the underlying currency).

     @param dFixedFXRate The fixed FX rate.
   */
  void SetFixedFXRate(double dFixedFXRate);

  /**
     Defines the fixed quanto with the FX volatility and the correlation
     between the underlying share and the FX rate.

     @param dFXVolatility The volatility of the FX rate.
     @param dCorrelation The correlation between the undelying share and
                         the FX rate. 
   */
  void SetFixedQuanto(double dFXVolatility, double dCorrelation);
  
  /**
     Whether the derivative is converted into new share.
     True if the conversion gives new shares, false otherwise.

     @param bConvertIntoNewShare new share flag.
   */
  void SetConvertIntoNewShare(bool bConvertIntoNewShare);

  //@}

  /**
      @name Methods for accessing convertible-like properties.
   */
  //@{

  // see base class
  virtual void PerturbFXRate(double dFXRateShift) const;
  
  /**
     Checks if the convertible is a fixed quanto or not.

     @return true if the bond is a fixed quanto.
   */
  bool IsFixedQuanto() const;

  /**
     Checks if the convertible is exchangeable or not.

     @return true if the bond is exchangeable.
   */
  bool IsExchangeable() const;

  /**
     In the case of an exchangeable, this function specifies whether the 
     exchangeable is exchangeable upon default. Otherwise, an exception is
     thrown.

     @return Exchangeable upon default property
   */
  bool IsExchangeableUponDefault() const;

  /**
     Specifies whether the conversion triggers are expressed as a percentage
     of the fixed or accreted conversion price.

     @return the conversion triggers as percentage of type.
   */
  TriggerAsPercentageOf GetConversionTriggerAsPercentageOf() const 
  { 
    return m_triggerAsPercentageOf; 
  }

  /**
     The trigger currency type (for call and conversion) in case of a
     cross currency: whether the trigger is in the security currency
     or in the underlying currency.

     @return The trigger currency type.
   */
  TriggerInCurrencyOf GetTriggerInCurrencyOf() const;
  
  /**
     Whether the derivative is converted into new share.
     True if the conversion gives new shares, false otherwise.

     @return new share flag.
   */
  bool GetConvertIntoNewShare() const { return m_bConvertIntoNewShare; }
  
  /**
     The Fixed FX Rate (for cross currency). It is the number of base currency 
     (i.e the derivative currency) needed to purchase one unit of the foreign
     currency (i.e the underlying currency).

     @return The fixed FX rate
   */
  double GetFixedFXRate() const;
  
  /**
     Gets the volatility of the FX Rate (for a fixed quanto).

     @return The volatility of the FX rate.
   */
  double GetFXRateVolatility() const;
  
  /**
     Gets the correlation between the underlying share and the FX rate 
     (for a fixed quanto).

     @return The correlation between the underlying share and the FX rate. 
   */
  double GetCorrelationBetweenUnderlyingAndFXRate() const;

  /**
     The terms of the bond like security.

     @return bond like terms
   */
  const shared_ptr<BondLikeTerms>& GetBondLikeTerms() const
  {
    return m_pBondLikeTerms;
  }

  /**
     Gets the accrued interest value at valuation date

     @return accrued interest value at valuation date
   */
  virtual double GetAccruedInterestValue() const;

  //@}

  // implements base class pure virtuals
  virtual Date GetMaturityDate() const;
  
  virtual void ValidateWith(const SessionData& sessionData) const;

protected:
  
  /**
     Creates a convertible-like instrument using BondLikeTerms.

     @param pBondLikeTerms The main characteristics of the bond like.
   */
  ConvertibleLike(const shared_ptr<BondLikeTerms>& pBondLikeTerms);

  /**
     Writes myself to tag parent which can be ConvertibleBond
     or Reset or some other derived class tag

     @param tagParent tag of the derived class
   */
  void DumpMe(XML::Tag& tagParent) const;

  /// bond like terms
  shared_ptr<BondLikeTerms> m_pBondLikeTerms;

  /// is exchangeable bond
  bool m_bIsExchangeable;

  /// conversion on default (only relevant for exchangeable bonds).
  bool m_bExchangeableUponDefault;

  /// Trigger currency type ( if the convertible is a cross currency)
  TriggerInCurrencyOf m_triggerInCurrencyOf;

  /// New share flag.
  bool m_bConvertIntoNewShare;

  /// Fixed quanto flag.
  bool m_bIsFixedQuanto;

  /// The fixed FX Rate (for cross currency)
  mutable double m_dFixedFXRate;

  /// The FX rate volatility
  double m_dFXRateVolatility;

  /** 
     The correlation between the underlying share and the FX rate 
     (for fixed quanto).
   */
  double m_dCorrelation;

  /**
     Of percentage of which value (fixed or accreted conversion price)
     is the trigger rate. Default to accreted conversion price.
   */
  TriggerAsPercentageOf m_triggerAsPercentageOf;

}; // class ConvertibleLike


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CONVERTIBLELIKE_H_
