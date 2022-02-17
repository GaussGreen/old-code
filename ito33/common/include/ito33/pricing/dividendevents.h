/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/dividendevents.h
// Purpose:     the dividend events 
// Author:      David Pooley
// Created:     2003/08/14
// RCS-ID:      $Id: dividendevents.h,v 1.13 2005/06/02 11:29:33 wang Exp $
// Copyright:   (c) 2003-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/dividendevents.h
    @brief The declaration of dividend event classes.

    Examples include cash dividends, pseudo cash dividend and yield dividends.
*/

#ifndef _ITO33_PRICING_DIVIDENDEVENTS_H_
#define _ITO33_PRICING_DIVIDENDEVENTS_H_

#include "ito33/pricing/dividendevent.h"
 
namespace ito33
{

namespace pricing
{


/** 
    Yield dividend class

    Yield dividends pay a percentage of the current stock price.
 */
class YieldDividendEvent : public DividendEvent
{

public:

  /** 
      Constructor

      @param dTime time of the dividend
      @param dRate the yield dividend rate
   */
  YieldDividendEvent(double dTime, double dRate)
                   : DividendEvent(dTime), m_dRate(dRate)
  {
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const;


protected:

  /// The dividend rate (e.g. 0.1 indicates a 10% dividend)
  double m_dRate;

};


/** 
    Pseudo cash dividend event class

    Pseudo cash dividends pay an absolute cash amount. If the amount gets close
    to causing a negative stock price, a yield payment is made.
*/
class PseudoCashDividendEvent : public DividendEvent
{

public:

  /** 
      Constructor

      @param dTime time of the dividend
      @param dAmount the absolute cash amount
      @param dPseudoYield Between 0 and 1. Pay a yield dividend of rate
             dPseudoYield if stock price is less than dAmount / dPseudoYield
   */
  PseudoCashDividendEvent(double dTime, double dAmount, double dPseudoYield)
    : DividendEvent(dTime), m_dAmount(dAmount), m_dPseudoYield(dPseudoYield)
  {
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const;


private:

  /// The absolute cash amount of the dividend
  double m_dAmount;

  /// Yield dividend rate to be used for small stock values
  double m_dPseudoYield;
};


/** 
    Cash dividend event class

    Cash dividends pay an absolute cash amount. 
 */
class CashDividendEvent : public DividendEvent
{

public:

  /** 
      Constructor

      @param dTime time of the dividend
      @param dAmount the absolute cash amount
   */
  CashDividendEvent(double dTime, double dAmount)
                  : DividendEvent(dTime), m_dAmount(dAmount)
  {
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const;


private:

  /// The absolute cash amount of the dividend
  double m_dAmount;
};


} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_DIVIDENDEVENTS_H_
