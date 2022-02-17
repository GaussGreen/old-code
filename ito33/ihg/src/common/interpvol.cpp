/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/interpvol.cpp
// Purpose:     Paramterized volatility surface class using interpolation
// Author:      David
// Created:     04.02.02
// RCS-ID:      $Id: interpvol.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/interpvol.cpp
*/

#include <math.h>

#include "ito33/numeric/interpolation.h"

#include "ihg/interpvol.h"

using namespace ito33::ihg;

double InterpVol::AverageValue(double /* dSpotInitial */)
{
  return 0.3;
}

void InterpVol::GetVols(double dTime, const double *Spots, double *Vols,
                             size_t nNbSpots) const 
{

  // Get the appropriate parameters from the surface
  Array<double> pdC(m_nNbPoints);
  GetParams(dTime, pdC.Get(), 0., 0.);

  // Interpolate the data
  ito33::numeric::QuadraticInterpolate(m_pdPoints.Get(), pdC.Get(), m_nNbPoints,
    Spots, Vols, nNbSpots, 
    ito33::numeric::ExtrapolationMode_Constant, ito33::numeric::ExtrapolationMode_Constant);
}

void InterpVol::GetParams(double dTime, double* pdParamValues, double, double) const
{
  // Get the appropriate interval
  FindInterval(dTime);

  // linearly interpolate the data.
  double dWeight = (dTime - *m_iterTimesLeft) / (*m_iterTimesRight - *m_iterTimesLeft);
  double dOneMinusWeight = 1.0 - dWeight;

  for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
    pdParamValues[nIdx] = (*m_iterParamsLeft)[nIdx] * dOneMinusWeight 
                        + (*m_iterParamsRight)[nIdx] * dWeight;
}
  

double InterpVol::GetSmoothnessMeasure(double dTime, double dS) 
{

  // Get the time derivatives at points (dS/2, dS, 2*dS).
  // Use the fact that linear interpolation in time is used.
  FindInterval(dTime);

  double dSpots[3];
  dSpots[0] = 0.5*dS;
  dSpots[1] = dS;
  dSpots[2] = 1.5*dS;

  double dTimeRight = *m_iterTimesRight;
  double dTimeLeft = *m_iterTimesLeft;

  double dVolRight[3], dVolLeft[3];
  GetVols(dTimeLeft, dSpots, dVolLeft, 3);
  GetVols(dTimeRight, dSpots, dVolRight, 3);

  // time deriv averaged at 3 points, with middle point having higher weight
  double dVolByTime = fabs(dVolRight[0] - dVolLeft[0]) 
                    + 2.0*fabs(dVolRight[1] - dVolLeft[1]) 
                    + fabs(dVolRight[2] - dVolLeft[2]);
  dVolByTime /= 4.0*(dTimeRight - dTimeLeft);

  // Scale to make a vol drop of 0.1 in 1 month comparable to 0.01% accuracy
  // in a single price (at least 3 digits correct). ie if the calibration is 
  // only to 1 price, then an accuracy of 0.1% in the price has about the same
  // weight to the above vol drop. This means that the user adjustable 
  // regulariation factor should have a magnitude around 1.
  dVolByTime /= 100.0;

  // An attempt at smoothness in the asset direction
  size_t nIdx;
  Array<double> pdVols(m_nNbPoints);
  GetParams(dTime, pdVols.Get(), 0., 0.);
  double dSum = 0.0;
  for (nIdx = 0; nIdx < m_nNbPoints-1; nIdx++)
  {
    dSum += fabs( (pdVols[nIdx + 1] - pdVols[nIdx]) / 
                  (m_pdPoints[nIdx+1] - m_pdPoints[nIdx]) );
  }

  return dSum/m_nNbPoints + dVolByTime;
}
