/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/polynomialhr.cpp
// Purpose:     Paramaterized hazard rate surface using polynomials
// Author:      David
// Created:     04/01/26
// RCS-ID:      $Id: polynomialhr.cpp,v 1.12 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "ihg/polynomialhr.h"

namespace ito33
{

namespace ihg
{

void PolynomialHR::GetHazardRates(double dTime,
                                  const double* pdSpots, 
                                  double *pdValues, 
                                  size_t nNumber) const
{

  FindInterval(dTime);

  double dWeight = (dTime - *m_iterTimesLeft) / (*m_iterTimesRight - *m_iterTimesLeft);
  double dOneMinusWeight = 1.0 - dWeight;

  for (size_t nIdx = 0; nIdx < nNumber; nIdx++)
  {
    double dTmp = pdSpots[nIdx] + 1.e-16;
    double dValLeft = (*m_iterParamsLeft)[0] + (*m_iterParamsLeft)[1] / dTmp;
    double dValRight = (*m_iterParamsRight)[0] + (*m_iterParamsRight)[1] / dTmp;

    pdValues[nIdx] = dOneMinusWeight * dValLeft + dWeight * dValRight;

    if (pdValues[nIdx] < 0.0) pdValues[nIdx] = 0.0;
    if (pdValues[nIdx] > 1.0) pdValues[nIdx] = 1.0;
  }

}

void PolynomialHR::GetParams(double dTime, double* pdParamValues, 
                             double /* dLeftBound */, double /* dRightBound */) const
{
  // Get the appropriate interval
  FindInterval(dTime);

  // SEE THE CODE FOR POLYNOMIALVOL IF THE NUMBER OF PARAMS INCREASES!!!!
  // linearly interpolate the data. Not the same as interplating the surface,
  // but is hopefully close enough. In particular, this will work for constant
  // surfaces.
  double dWeight = (dTime - *m_iterTimesLeft) / (*m_iterTimesRight - *m_iterTimesLeft);
  double dOneMinusWeight = 1.0 - dWeight;

  for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
    pdParamValues[nIdx] = (*m_iterParamsLeft)[nIdx] * dOneMinusWeight 
                        + (*m_iterParamsRight)[nIdx] * dWeight;

}




double PolynomialHR::GetSmoothnessMeasure(double dTime, double /* dS */)
{
  // Make sure we are in the proper interval
  FindInterval(dTime);

  // Since we are using linear interpolation, the time component drops out.
  // We can just use the endpoints of the current interval
  double dTimeInterval = (*m_iterTimesRight) - (*m_iterTimesLeft);
  double dTimeIntervalInv = 1.0/dTimeInterval;

  // Get the derivatives w.r.t. the paramaters
  double dC1bydt = ( (*m_iterParamsRight)[0] - (*m_iterParamsLeft)[0] ) * dTimeIntervalInv;

  // Now apply the formula...Actually, it is not clear what S value to use,
  // so just return the sum of the derivatives of the parameters

  double dTmp = fabs(dC1bydt);
  return dTmp;

}

} // namespace ihg

} // namespace ito33
