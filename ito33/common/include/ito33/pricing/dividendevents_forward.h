/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/dividendevents_forward.h
// Purpose:     the dividend events for forward pricer
// Author:      ICARE
// Created:     2003/11/17
// RCS-ID:      $Id: dividendevents_forward.h,v 1.15 2005/06/29 20:02:35 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/dividendevents_forward.h
    @brief The declaration of dividend event classes for forward pricer.

    Examples include fixed dividends and proportinal dividends.
*/

#ifndef _ITO33_PRICING_DIVIDENDEVENTS_FORWARD_H_
#define _ITO33_PRICING_DIVIDENDEVENTS_FORWARD_H_

#include "ito33/pricing/dividendevent.h"
 
namespace ito33
{

namespace pricing
{


/// cash dividend class
class CashDividendEvent_Forward : public DividendEvent
{

public:

  /** 
      Constructor

      @param dTime time of the dividend
      @param dCash the dividend amount
   */
  CashDividendEvent_Forward(double dTime, double dCash)
    : DividendEvent(dTime), m_dCash(dCash)
  {
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const;


private:

  /// The dividend amount
  double m_dCash;

};

/**
    Proportinal dividend class

    Proportional dividends pay a percentage of the current stock price.
 */
class YieldDividendEvent_Forward : public DividendEvent
{

public:

  /**
      Constructor

      @param dTime time of the dividend
      @param dRate the proportional dividend rate
   */
  YieldDividendEvent_Forward(double dTime, double dRate)
    : DividendEvent(dTime), m_dRate(dRate)
  {
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const;

  void 
  ApplyToPrice(const double *pdS, double* pdValues, size_t nNbS) const;

  virtual numeric::InterpolationMatrix* 
  GetInterpolationMatrix
  (const double* pdS, size_t nNbS, size_t nNbSubSystem) const;


private:

  /// The dividend rate (e.g. 0.1 indicates a 10% dividend)
  double m_dRate;

};


/** Fixed dividend class

    Fixed dividends pay a fixed dividend amount. If the amount gets close
    to causing a negative stock price, a proportional payment is made.
*/
class PseudoCashDividendEvent_Forward : public DividendEvent
{

public:

  /**
      Constructor

      @param dTime time of the dividend
      @param dAmount amount of the fixed dividend
      @param dPseudoYield Between 0 and 1. Pay a proportional dividend of rate
             dPseudoYield if stock price is less than dAmount/dPseudoYield
  */
  PseudoCashDividendEvent_Forward
     (double dTime, double dAmount, double dPseudoYield)
    : DividendEvent(dTime), m_dAmount(dAmount), m_dPseudoYield(dPseudoYield)
  {
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const;

  void 
  ApplyToPrice(const double *pdS, double* pdValues, size_t nNbS) const;


private:

  /// Amount of the fixed dividend
  double m_dAmount;

  /// Proportional dividend rate to be used for small stock values
  double m_dPseudoYield;
};


} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_DIVIDENDEVENTS_FORWARD_H_
