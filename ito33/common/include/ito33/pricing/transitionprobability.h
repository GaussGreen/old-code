/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/transitionprobability.h
// Purpose:     Contracts class for transition probability
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobability.h,v 1.1 2006/03/31 17:43:48 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/transitionprobability.h
    @brief Contracts class for transition probability.
 */

#ifndef _ITO33_PRICING_TRANSITIONPROBABILITY_H_
#define _ITO33_PRICING_TRANSITIONPROBABILITY_H_

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace pricing
{


/// Contracts class for transition probability.
class TransitionProbability : public Contract 
{
public:
  
  TransitionProbability(double dStrike, double dMaturityTime)
    : m_dStrike(dStrike)
  {
    m_dMaturityTime = dMaturityTime;
  }
  
  double GetStrike() const { return m_dStrike; }

private:

  double m_dStrike;

  NO_COPY_CLASS(TransitionProbability);

}; // class TransitionProbability;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_TRANSITIONPROBABILITY_H_
