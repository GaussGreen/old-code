/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/option.h
// Purpose:     financial option class
// Author:      zhang
// Created:     2003/11/05
// RCS-ID:      $Id: option.h,v 1.37 2006/08/19 09:37:17 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/option.h
    @brief declaration of the financial option class
 */

#ifndef _ITO33_FINANCE_OPTION_H_
#define _ITO33_FINANCE_OPTION_H_

#include "ito33/date.h"

#include "ito33/finance/optionlike.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    Option represents the financial aspects of an option.

    Classes corresponding to more complicated options (barrier and other
    exotics) derive from this one.
 */
class ITO33_DLLDECL Option : public OptionLike
{
public:
  /**
     Creates an option object.

     Option object is immutable once it was created, i.e. its description
     can't be changed.

     @param dStrike the strike value of the option which must be positive
     @param maturityDate the maturity date of the option
     @param optionType the type of the option
     @param exerciseType the exercise type of the option
   */
  Option(double dStrike,
         Date maturityDate,
         OptionType optionType,
         ExerciseType exerciseType);

  /**
     Empty virtual destructor.

     Option can be derived by barrier option
   */
  virtual ~Option() { }

  /**
      @name Methods for accessing the Option.
   */
  //@{
 
  /**
     Get the strike of the option.

     @return strike of the option
   */
  double GetStrike() const
  {
    return m_dStrike;
  }

  /**
     Indicates whether the implied volatility has been set.

     @return true if the implied volatility of the option is set, 
             false otherwise.
   */
  bool IsImpliedVolSet() const
  {
    return m_dImpliedVol >= 0.;
  }

  /**
     The implied BS volatility of the option.

     @return implied volatility of the option
   */
  double GetImpliedVol() const;

  /**
     Gets the implied volatility of the option for a given market price.

     @param dPrice market price

     @return implied vol that best fit the market price
   */
  double GetImpliedVolFrom(double dPrice) const;

  /**
     Gets the price of the option for a given implied Black-Scholes vol.

     @param dVol the implied Black-Scholes vol of this option

     @return the option price corresponds to the implied Black-Scholes vol
   */
  double GetPriceFrom(double dVol) const;

  /**
     Gets the BS vega for a given implied BS vol.

     @param dVol the implied Black-Scholes vol of this option
     
     @return the BS vega at the given implied vol
   */
  double GetVegaFrom(double dVol) const;

  //@}
 
  /**
     The implied BS volatility of the option.

     @param dImpliedVol implied volatility
   */
  void SetImpliedVol(double dImpliedVol);

  // implement base class pure virtuals
  virtual void Visit(DerivativeVisitor& visitor) const;
  void Visit(DerivativeModifyingVisitor& visitor);
  virtual XML::Tag Dump(XML::Tag& tagParent) const;

  // overload MarketPrice related functions as it can be defined 
  // by implied volatility
  virtual bool HasMarketPrice() const;
  virtual void SetMarketPrice(double dPrice);

protected:
  virtual double DoGetMarketPrice() const;

protected:

  /// strike of the option
  double m_dStrike;

  /// implied volatility
  double m_dImpliedVol;

}; // class Option


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_OPTION_H_
