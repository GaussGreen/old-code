/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/options.h
// Purpose:     contracts class for a list of European options
// Author:      Wang
// Created:     2004/02/27
// RCS-ID:      $Id: options.h,v 1.14 2006/08/07 15:12:44 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/options.h
    @brief The declaration of the Options class.

    Options contain a list of options 
 */

#ifndef _ITO33_PRICING_OPTIONS_H_
#define _ITO33_PRICING_OPTIONS_H_

#include "ito33/vector.h"

#include "ito33/finance/optiontype.h"

#include "ito33/pricing/option.h"

namespace ito33
{

namespace finance
{
  class ForwardOption;
}

namespace pricing
{

/// The declaration of the Options class.
class Options : public Contracts 
{
public:

  /**
      ctor

      @param forwardOption a forward option
   */
  Options(const finance::ForwardOption& forwardOption);

  // default dtor is ok

  /// Get maturity
  double GetMaturityTime() const;

  /// Get special times
  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

  /// Get vector or maturity dates
  const std::vector<Date>& GetMaturityDates() { return m_pMaturityDates; }

  /// Get vector of strikes
  const std::vector<double>& GetStrikes() { return m_pdStrikes; }

  /// Get vector of market prices
  const std::vector<double>& GetMarketPrices() { return m_pdMarketPrices; }

  /// Get vector of weights
  const std::vector<double>& GetWeights() const 
  { 
    return m_pdWeights; 
  }

  /// Get vector of option types
  const std::vector<finance::OptionType>& GetOptionTypes() const
  { 
    return m_pOptionTypes; 
  }

protected:

  std::vector<Date> m_pMaturityDates;

  std::vector<double> m_pdStrikes;

  std::vector<double> m_pdMarketPrices;

  std::vector<double> m_pdWeights;

  std::vector<finance::OptionType> m_pOptionTypes;

}; // class Options;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_OPTIONS_H_

