/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/payoffoption.h
// Purpose:     payoff implementations for option 
// Author:      WANG Xuewen
// Created:     2003/09/25
// RCS-ID:      $Id: payoffoption.h,v 1.10 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/payoffoption.h
    @brief The implementation of the option payoff classes.
 */

#ifndef _ITO33_FINANCE_PAYOFFOPTION_H_
#define _ITO33_FINANCE_PAYOFFOPTION_H_

#include "ito33/finance/payoff.h"

namespace ito33 
{

namespace finance 
{


/// The base class for option payoff
class PayoffOption : public Payoff
{
public:

  /**
     constructor, set the strike
     
     @param dStrike The strike of the option
   */
  PayoffOption(double dStrike) : m_dStrike(dStrike) { }
 
  /// dummy virtual destructor
  virtual ~PayoffOption() { }

  /**
     Get the option payoff values at an array of spots

     @param pdS the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  virtual void Get(const double *pdS, double *pdPrices, size_t nNbS) const = 0;

protected:

  /// the strike of the option 
  double
    m_dStrike;
};

/// Payoff class for a call option
class PayoffCall : public PayoffOption
{
public:
 
  /**
     constructor, set the strike
     @param dStrike The strike of the call option
   */
  PayoffCall(double dStrike) : PayoffOption(dStrike) { }

  /**
     Get the option payoff values at an array of spots

     @param pdS the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double *pdS, double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdxS = 0; nIdxS < nNbS; ++nIdxS)
      pdPrices[nIdxS] = (pdS[nIdxS] < m_dStrike) ? 0.0 : pdS[nIdxS] - m_dStrike;
  }
};

/// Payoff class for a put option
class PayoffPut : public PayoffOption
{
public:

  /**
     constructor, set the strike
   
     @param dStrike The strike of the put option
   */  
  PayoffPut(double dStrike) : PayoffOption(dStrike) { }

  /**
     Get the put option payoff values at an array of spots

     @param pdS the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double *pdS, double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdxS = 0; nIdxS < nNbS; ++nIdxS)
      pdPrices[nIdxS] = (pdS[nIdxS] > m_dStrike) ? 0.0 : m_dStrike - pdS[nIdxS];
  }
};

/// payoff class for a digital option
class PayoffDigital : public PayoffOption
{
public:

  /**
     constructor, set the strike
     
     @param dStrike The strike of the digital option
   */  
  PayoffDigital(double dStrike) : PayoffOption(dStrike) { }

  /**
     Get the digital option payoff values at an array of spots

     @param pdS the spots
     @param pdPrices the payoff values at the spots,
     @param nNbS the number of spots
   */
  void Get(const double *pdS, double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdxS = 0; nIdxS < nNbS; ++nIdxS)
    {
      pdPrices[nIdxS] = (pdS[nIdxS] > m_dStrike) ? 1.0 : 0.0;
     
      if (pdS[nIdxS] == m_dStrike)
        pdPrices[nIdxS] = 0.5;
    }
  }
};


/// Payoff class for a fix call payoff
class PayoffFixCall : public PayoffOption
{
public:
 
  /**
     constructor, set the strike max(A-K,0.0)

     @param dStrike The strike of the fix call payoff
     @param dAverage the average of the fix call payoff

   */
  PayoffFixCall(double dStrike, double dAverage) : PayoffOption(dStrike),
    m_dAverage(dAverage) 
  { 
  }

  /**
     Get the option payoff values at an array of spots

     @param  the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double * , double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdxS = 0; nIdxS < nNbS; ++nIdxS)
      pdPrices[nIdxS] = m_dAverage < m_dStrike ? 0.: m_dAverage - m_dStrike;
  }

private:
  ///Average of the fix call payoff
  double m_dAverage;

};


/// Payoff class for a fix put payoff
class PayoffFixPut : public PayoffOption
{
public:
 
  /**
     constructor, set the strike max(K-A,0.0)

     @param dStrike The strike of the fix put payoff
     @param dAverage the average of the fix put payoff

   */
  PayoffFixPut(double dStrike, double dAverage) : PayoffOption(dStrike),
  m_dAverage(dAverage)
  { 
  }

  /**
     Get the option payoff values at an array of spots

     @param  the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double * , double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdxS = 0; nIdxS < nNbS; ++nIdxS)
      pdPrices[nIdxS] =  m_dStrike < m_dAverage ? 0.: m_dStrike - m_dAverage;
  }

private:
  ///Average of the fix call payoff
  double m_dAverage;

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_PAYOFFOPTION_H_
