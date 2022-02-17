/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/varianceswaption.h
// Purpose:     financial class for variance swaption
// Created:     2006/07/21
// RCS-ID:      $Id: varianceswaption.h,v 1.2 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/varianceswaption.h
    @brief Declaration of the financial class for variance swaption
 */

#ifndef _ITO33_FINANCE_VARIANCESWAPTION_H_
#define _ITO33_FINANCE_VARIANCESWAPTION_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/optiontype.h"
#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL VarianceSwapTerms;

/**
    A variance swaption gives its holder the right, but not the obligation,
    to write at a future date (the maturity date of the swaption) a variance
    swap whose volatility strike is just the strike of the variance swaption.
 */
class ITO33_DLLDECL VarianceSwaption : public Derivative
{
public:
  /**
      Creates a variance swaption.
      
      @param pTerms The variance swap terms of the swaption
      @param optionType The type of the variance swaption (call/put)
      @param dStrike The strike of the swaption
      @param maturityDate The maturity date of the variance swaption
   */
  VarianceSwaption(const shared_ptr<VarianceSwapTerms>& pTerms,
                   OptionType optionType, double dStrike, Date maturityDate);

  // Default dtor is ok

  /// @name Accessors for variance swaption.
  //@{

  /**
      Gets the variance swap terms of the swaption.

      @return The variance swap terms of the swaption
   */
  const shared_ptr<VarianceSwapTerms>& GetTerms() const { return m_pTerms; }

  /**
      Gets the option type of the variance swaption.

      @return The option type of the variance swaption.
   */
  OptionType GetOptionType() const { return m_optionType; }

  /**
      Gets the strike of the variance swaption.

      @return the volatility strike of the swap
   */
  double GetStrike() const { return m_dStrike; }

  /**
      Gets the maturity date of the variance swaption.

      @return maturity date of the variance swaption
   */
  Date GetMaturityDate() const { return m_maturityDate; }

  //@}

  // base class virtual functions
  virtual void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// The variance swap terms of the swaption
  shared_ptr<VarianceSwapTerms> m_pTerms;

  /// The type of the variance swaption(call/put)
  OptionType m_optionType;

  /// The strike
  double m_dStrike;

  /// The maturity date
  Date m_maturityDate;

}; // class VarianceSwaption


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_VARIANCESWAPTION_H_
