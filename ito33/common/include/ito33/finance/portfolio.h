/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/portfolio.h
// Purpose:     Portfolio class declaration
// Author:      Vadim Zeitlin
// Created:     13.08.03
// RCS-ID:      $Id: portfolio.h,v 1.4 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/portfolio.h
    @brief Declaration of Portfolio class.

    Portfolio represents all the financial instruments the client is interested
    in. It can contain arbitrarily many Derivatives but has only a single
    Currency and all the Derivatives must have the same Underlying which is
    also stored in this class.
 */

#ifndef _FINANCE_PORTFOLIO_H_
#define _FINANCE_PORTFOLIO_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

namespace ito33
{

class Numeraire;
class Derivative;
class Underlying;

/**
    Portfolio: collection of Derivatives using the same Underlying.
 */
class Portfolio
{
public:
  /**
      Create a new portfolio.

      This constructor creates an empty portfolio, use AddInstrument() later to
      add stuff to it.

      @param pUnderlying the underlying common to all our derivatives
      @param pNumeraire the base currency for this portfolio
   */
  Portfolio(const shared_ptr<Underlying>& pUnderlying,
            const shared_ptr<Numeraire>& pNumeraire)
    : m_pUnderlying(pUnderlying),
      m_pNumeraire(pNumeraire)
  {
    // we might want to call m_derivatives.reserve() here
  }


  /**
      Adds a new financial instrument to the portfolio.

      @param pDerivative the instrument wrapper in a ref counting pointer
   */
  void AddInstrument(const shared_ptr<Derivative>& pDerivative)
  {
    m_derivatives.push_back(pDerivative);
  }

private:
  /// the underlying for all the derivatives of this portfolio
  shared_ptr<Underlying> m_pUnderlying;

  /// the currency of this portfolio
  shared_ptr<Numeraire> m_pNumeraire;

  /// all the derivatives
  std::vector< shared_ptr<Derivative> > m_derivatives;
};

} // namespace ito33

#endif // _FINANCE_PORTFOLIO_H_

