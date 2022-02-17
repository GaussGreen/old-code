/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/dividendevents.cpp
// Purpose:     implementation of the dividend events
// Author:      David Pooley
// Created:     2003/10/06
// RCS-ID:      $Id: dividendevents.cpp,v 1.19 2005/06/02 11:58:13 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @todo Should the interpolation scheme be the same, this code can be greatly
         reduced so that ApplyToPrice can be implememted in the base class
 */

#include "ito33/pricing/dividendevents.h"

namespace ito33
{

namespace pricing
{
  

void 
YieldDividendEvent::ApplyToSpots
(const double* pdS, double* pdNewS, size_t nNbS) const
{
  // Shift the S mesh by the amount of the dividend
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdNewS[nIdx] = pdS[nIdx] * (1. - m_dRate);
}
  
void 
PseudoCashDividendEvent::ApplyToSpots
(const double* pdS, double* pdNewS, size_t nNbS) const
{
  // Shift the S mesh by the amount of the dividend
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    if (pdS[nIdx] > m_dAmount / m_dPseudoYield)
      pdNewS[nIdx] = pdS[nIdx] - m_dAmount;
    else
      pdNewS[nIdx] = pdS[nIdx] * (1. - m_dPseudoYield);
  }
}

void 
CashDividendEvent::ApplyToSpots
(const double* pdS, double* pdNewS, size_t nNbS) const
{
  // Shift the S mesh by the amount of the dividend
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbS; nIdx++)
    pdNewS[nIdx] = (pdS[nIdx] > m_dAmount) ? pdS[nIdx] - m_dAmount : 0;
}


} // namespace pricing

} // namespace ito33
