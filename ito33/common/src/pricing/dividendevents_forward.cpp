/////////////////////////////////////////////////////////////////////////////
// Name:        common/src//pricing/dividendevents_forward.cpp
// Purpose:     implementation of the dividend events for forward equation
// Author:      ICARE
// Created:     2003/10/06
// RCS-ID:      $Id: dividendevents_forward.cpp,v 1.18 2005/06/29 20:05:53 dave Exp $
// Copyright:   (c) 2003 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/array.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"

#include "ito33/pricing/dividendevents_forward.h"

namespace ito33
{

namespace pricing
{
  
  using namespace numeric;

void 
CashDividendEvent_Forward::ApplyToSpots
(const double* pdS, double* pdNewS, size_t nNbS) const
{
  // Shift the S mesh by the amount of the dividend
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdNewS[nIdx] = pdS[nIdx] + m_dCash;
}

void 
YieldDividendEvent_Forward::ApplyToSpots
(const double* pdS, double* pdNewS, size_t nNbS) const
{
  double dInvRate = 1. / (1. - m_dRate);

  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdNewS[nIdx] = pdS[nIdx] * dInvRate;
}

/// Function to apply the event to pdValues on grid pdS
void 
YieldDividendEvent_Forward::ApplyToPrice
(const double* pdS, double* pdValues, size_t nNbS) const
{
  // Create temporary arrays for interpolation
  Array<double> pdSTmp(nNbS);
  Array<double> pdValuesTmp(nNbS);

  ApplyToSpots( pdS, pdSTmp.Get(), nNbS );

  Interpolate(pdS, pdValues, nNbS, 
              pdSTmp.Get(), pdValuesTmp.Get(), nNbS, 
              m_emLeft, m_emRight, m_interpolationMethod);

  // use the interpolated values to calculate the real values
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdValues[nIdx] = pdValuesTmp[nIdx] * (1. - m_dRate);
}

InterpolationMatrix* 
YieldDividendEvent_Forward::GetInterpolationMatrix
(const double* pdS, size_t nNbS, size_t nNbSubSystem) const
{
  InterpolationMatrix* 
    pMatrix = DividendEvent::GetInterpolationMatrix(pdS, nNbS, nNbSubSystem);

  pMatrix->MultiplyBy(1. - m_dRate);

  return pMatrix;
}

void 
PseudoCashDividendEvent_Forward::ApplyToSpots
(const double* pdS, double* pdNewS, size_t nNbS) const
{
  const double dSTmp = m_dAmount / m_dPseudoYield;
  const double dInvRate = 1. / (1. - m_dPseudoYield);

  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    if (pdS[nIdx] > dSTmp)
      pdNewS[nIdx] = pdS[nIdx] + m_dAmount;
    else
      pdNewS[nIdx] = pdS[nIdx] * dInvRate;
  }
}

/// Function to apply the event to pdValues on grid pdS
void 
PseudoCashDividendEvent_Forward::ApplyToPrice
(const double* pdS, double* pdValues, size_t nNbS) const
{
  // Get the value at the special point m_dAmount / m_dPseudoYield
  double dSTmp;
  double dValueTmp;
   
  dSTmp = m_dAmount / m_dPseudoYield;

  Interpolate(pdS, pdValues, nNbS,
              &dSTmp, &dValueTmp, 1,
              m_emLeft, m_emRight, m_interpolationMethod);
   
  dValueTmp *= m_dPseudoYield;

  dSTmp *= (1. - m_dPseudoYield);
    
  // Create temporary arrays for interpolation
  Array<double> pdSTmp(nNbS);
  Array<double> pdValuesTmp(nNbS);

  ApplyToSpots( pdS, pdSTmp.Get(), nNbS );

  // Now interpolate from (pdS,pdValues) to (pdSTmp, pdValuesTmp)  
  Interpolate(pdS, pdValues, nNbS,
              pdSTmp.Get(), pdValuesTmp.Get(), nNbS,
              m_emLeft, m_emRight, m_interpolationMethod);

  // use the interpolated values to calculate the real values
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    if (pdS[nIdx] > dSTmp)
      pdValues[nIdx] = pdValuesTmp[nIdx];
    else
      pdValues[nIdx] = (1. - m_dPseudoYield) * pdValuesTmp[nIdx] + dValueTmp;
}


} // namespace pricing

} // namespace ito33
