/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/polynomialvol.cpp
// Purpose:     Paramterized volatility surface class using polynomials
// Author:      David
// Created:     04.01.24
// RCS-ID:      $Id: polynomialvol.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/polynomialvol.cpp
*/

#include <math.h>
#include <nag.h>
#include <nage04.h>

#include "ito33/array.h"

#include "ihg/polynomialvol.h"

using namespace ito33::ihg;

void __stdcall objfunPolyVol(Integer iN, double pdX[], double* objf, double pdG[], Nag_Comm* comm);

ito33::ihg::PolynomialVol* globalPolyVolPointer;

double PolynomialVol::AverageValue(double /* dSpotInitial */)
{
  return 0.3;
}

void PolynomialVol::GetVols(double dTime, const double *Spots, double *Vols,
                            size_t nNbSpots)
{

  FindInterval(dTime);

  double dWeight = (dTime - *m_iterTimesLeft) / (*m_iterTimesRight - *m_iterTimesLeft);
  double dOneMinusWeight = 1.0 - dWeight;

  for (size_t nIdx = 0; nIdx < nNbSpots; nIdx++)
  {
    double dTmp1 = pow(Spots[nIdx], (*m_iterParamsLeft)[4]);
    double dTmp2 = pow(Spots[nIdx], (*m_iterParamsLeft)[3]);
    double dVolLeft = (*m_iterParamsLeft)[0] + (*m_iterParamsLeft)[1] / dTmp1 
                    + (*m_iterParamsLeft)[2] * dTmp2;

    dTmp1 = pow(Spots[nIdx], (*m_iterParamsRight)[4]);
    dTmp2 = pow(Spots[nIdx], (*m_iterParamsRight)[3]);
    double dVolRight = (*m_iterParamsRight)[0] + (*m_iterParamsRight)[1] / dTmp1 
                     + (*m_iterParamsRight)[2] * dTmp2;

    Vols[nIdx] = dOneMinusWeight * dVolLeft + dWeight * dVolRight;

    if (Vols[nIdx] < 0.0) Vols[nIdx] = 0.0;
    if (Vols[nIdx] > 2.0) Vols[nIdx] = 2.0;
  }
}


void PolynomialVol::GetParams(double dTime, double* pdParamValues,
                              double dLeftBound, double dRightBound)
{
  // Get the appropriate interval
  FindInterval(dTime);

  // Special case.  At interval boundary
  if ( fabs(*m_iterTimesLeft - dTime) < 1.e-14 )
  {
    pdParamValues[0] = (*m_iterParamsLeft)[0];
    pdParamValues[1] = (*m_iterParamsLeft)[1];
    pdParamValues[2] = (*m_iterParamsLeft)[2];
    pdParamValues[3] = (*m_iterParamsLeft)[3];
    pdParamValues[4] = (*m_iterParamsLeft)[4];
    return;
  }
  else if ( fabs(*m_iterTimesLeft - dTime) < 1.e-14 )
  {
    pdParamValues[0] = (*m_iterParamsRight)[0];
    pdParamValues[1] = (*m_iterParamsRight)[1];
    pdParamValues[2] = (*m_iterParamsRight)[2];
    pdParamValues[3] = (*m_iterParamsRight)[3];
    pdParamValues[4] = (*m_iterParamsRight)[4];
    return;
  }


  // The requested time is within the interval.  Use NAG to do a least
  // squares fit to the interpolated volatility (NOTE: interpolating
  // the parameters directly is NOT the same as fitting the interpolated
  // volatility values).

  // Setup the points and vols to fit
  m_nNbPointsToFit = 10;
  m_pdVolsToFit = Array<double>(m_nNbPointsToFit);
  m_pdPointsToFit = Array<double>(m_nNbPointsToFit);
  size_t nIdx;
  double dInterval = (dRightBound - dLeftBound)/(m_nNbPointsToFit - 1);
  for (nIdx = 0; nIdx < m_nNbPointsToFit; nIdx++)
    m_pdPointsToFit[nIdx] = dLeftBound + nIdx * dInterval;

  GetVols(dTime, m_pdPointsToFit.Get(), m_pdVolsToFit.Get(), m_nNbPointsToFit);

  // Setup NAG
  globalPolyVolPointer = this;

  double dFinalError;
  Array<double> pdG(5);
  Array<double> pdLowerBounds(5);
  Array<double> pdUpperBounds(5);

  Nag_BoundType nagBound;
  NagError nagError;
  Nag_E04_Opt nagOptions;

  e04xxc(&nagOptions);
  SET_FAIL(nagError);
  nagOptions.optim_tol = 1.e-4;
  nagOptions.f_est = 0.0;
  nagOptions.local_search = FALSE;
  nagOptions.max_iter = 50;

  #ifndef NDEBUG
    //nagOptions.print_level = Nag_Soln_Iter_Full;
    nagOptions.print_level = Nag_Soln_Iter;
    nagOptions.list = TRUE;
    nagError.print = FALSE;
  #else
    nagOptions.print_level = Nag_NoPrint;
    nagOptions.list = FALSE;
    nagError.print = FALSE;
  # endif

  // Initial guess, lower bounds, upper bounds
  nagBound = Nag_Bounds;
    
  GetParamBounds(pdLowerBounds.Get(), pdUpperBounds.Get());
  pdParamValues[0] = (*m_iterParamsRight)[0];
  pdParamValues[1] = (*m_iterParamsRight)[1];
  pdParamValues[2] = (*m_iterParamsRight)[2];
  pdParamValues[3] = (*m_iterParamsRight)[3];
  pdParamValues[4] = (*m_iterParamsRight)[4];

  // Do the optimization
  Integer n = 5;
  e04jbc(n, objfunPolyVol, nagBound, pdLowerBounds.Get(), pdUpperBounds.Get(), 
      pdParamValues, &dFinalError, pdG.Get(), &nagOptions, NAGCOMM_NULL, &nagError);

}

void __stdcall PolynomialVol::ParamFitObjFun(Integer /* iN */, double pdX[], double* objf, 
                                      double* /* pdG[] */, Nag_Comm* /* comm */)
{

  // Get the vol values using the current guess
  double dTotalError = 0.0;
  for (size_t nIdx = 0; nIdx < m_nNbPointsToFit; nIdx++)
  {
    double dTmp1 = pow(m_pdPointsToFit[nIdx], pdX[4]);
    double dTmp2 = pow(m_pdPointsToFit[nIdx], pdX[3]);
    double dVol = pdX[0] + pdX[1] / dTmp1 + pdX[2] * dTmp2;

    if (dVol < 0.0) dVol = 0.0;
    if (dVol > 2.0) dVol = 2.0;

    // Do a least squares fit
    double dErr = dVol - m_pdVolsToFit[nIdx];
    dTotalError += dErr * dErr;
  }

  // Return the total error
  *objf = dTotalError;

}

/*
void PolynomialVol::GetParams(double dTime, double* pdParamValues) const
{
  // Get the appropriate interval
  FindInterval(dTime);

  // If the polynomial powers are large, just return the closest values
  // If needed in the future, this should be changed to some type
  // of least squares fit to compute the parameters.
  // If the powers are small, the surface is flat, and it should be safe
  // to just average the values
  const double dLargePower = 0.01;

  if ( (*m_iterParamsLeft)[3] > dLargePower || (*m_iterParamsLeft)[4] > dLargePower 
    || (*m_iterParamsRight)[3] > dLargePower || (*m_iterParamsRight)[4] > dLargePower )
  {
    // Just return the closest values
    if (*m_iterTimesRight - dTime < dTime - *m_iterTimesLeft)
    {
      pdParamValues[0] = (*m_iterParamsRight)[0];
      pdParamValues[1] = (*m_iterParamsRight)[1];
      pdParamValues[2] = (*m_iterParamsRight)[2];
      pdParamValues[3] = (*m_iterParamsRight)[3];
      pdParamValues[4] = (*m_iterParamsRight)[4];
    }
    else
    {
      pdParamValues[0] = (*m_iterParamsLeft)[0];
      pdParamValues[1] = (*m_iterParamsLeft)[1];
      pdParamValues[2] = (*m_iterParamsLeft)[2];
      pdParamValues[3] = (*m_iterParamsLeft)[3];
      pdParamValues[4] = (*m_iterParamsLeft)[4];
    }
  }
  else
  {
    // linearly interpolate the data. Not the same as interplating the surface,
    // but is hopefully slose enough. In particular, this will work for constant
    // surfaces.
    double dWeight = (dTime - *m_iterTimesLeft) / (*m_iterTimesRight - *m_iterTimesLeft);
    double dOneMinusWeight = 1.0 - dWeight;

    for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
      pdParamValues[nIdx] = (*m_iterParamsLeft)[nIdx] * dOneMinusWeight 
                          + (*m_iterParamsRight)[nIdx] * dWeight;
  }
}
*/
  
double PolynomialVol::GetSmoothnessMeasure(double dTime, double dS)
{

  // Get the time derivatives at points (dS/2, dS, 2*dS).
  // Use the fact that linear interpolation in time is used.
  FindInterval(dTime);

  double dSpots[3];
  dSpots[0] = 0.5*dS;
  dSpots[1] = dS;
  dSpots[2] = 2.0*dS;

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

  // Scale to make a vol drop of 0.1 in 1 month comparable to 0.1% accuracy
  // in a single price (at least 3 digits correct). ie if the calibration is 
  // only to 1 price, then an accuracy of 0.1% in the price has about the same
  // weight to the above vol drop. This means that the user adjustable 
  // regulariation factor should have a magnitude around 1.
  return dVolByTime/10.0;

  /*
  // Get the current params. Also sets the correct interval
  double pdC[5];
  GetParams(dTime, pdC);

  // Since we are using linear interpolation, the time component drops out.
  // We can just use the current point and an endpoint of the current interval
  // to get derivatives of params w.r.t. time
  double dC1bydt, dC2bydt, dC3bydt, dC4bydt, dC5bydt;
  if ( fabs(dTime - *m_iterTimesLeft) < 1.e-14)
  {
    // Special case...on the left boundary. In normal usage, expect
    // to be on the right boundary, but we never know.
    double dTimeInterval = (*m_iterTimesRight) - dTime;
    double dTimeIntervalInv = 1.0/dTimeInterval;

    dC1bydt = ( (*m_iterParamsRight)[0] - pdC[0] ) * dTimeIntervalInv;
    dC2bydt = ( (*m_iterParamsRight)[1] - pdC[1] ) * dTimeIntervalInv;
    dC3bydt = ( (*m_iterParamsRight)[2] - pdC[2] ) * dTimeIntervalInv;
    dC4bydt = ( (*m_iterParamsRight)[3] - pdC[3] ) * dTimeIntervalInv;
    dC5bydt = ( (*m_iterParamsRight)[4] - pdC[4] ) * dTimeIntervalInv;
  }
  else
  {
    double dTimeInterval = dTime - (*m_iterTimesLeft);
    double dTimeIntervalInv = 1.0/dTimeInterval;

    dC1bydt = ( pdC[0] - (*m_iterParamsLeft)[0] ) * dTimeIntervalInv;
    dC2bydt = ( pdC[1] - (*m_iterParamsLeft)[1] ) * dTimeIntervalInv;
    dC3bydt = ( pdC[2] - (*m_iterParamsLeft)[2] ) * dTimeIntervalInv;
    dC4bydt = ( pdC[3] - (*m_iterParamsLeft)[3] ) * dTimeIntervalInv;
    dC5bydt = ( pdC[4] - (*m_iterParamsLeft)[4] ) * dTimeIntervalInv;
  }
  double dTmp = 0.0;

  // Want to get derivative of surface w.r.t. time, but it is not clear
  // which S value to use, or how to average 
  // Just return the sum of the derivatives of the parameters
  dTmp = fabs(dC1bydt) + fabs(dC2bydt) + fabs(dC3bydt) + fabs(dC4bydt)
              + fabs(dC5bydt);

  // Add in a measure of the curviness w.r.t S. This is ad hoc, but should 
  // help get flat surfaces. Terms 3 and 4 are the polynomial powers,
  // so we want them small. Terms 1 and 2 are the polynomial coefficients,
  // so we want them small, but don't care as much.  Term 0 is the
  // constant part, so we don't care about it at all
  dTmp += ( fabs(pdC[1]) + fabs(pdC[2]) ) / 100.0 + (pdC[3] + pdC[4]) / 10.0;

  return dTmp;
  */
}

void __stdcall objfunPolyVol(Integer iN, double pdX[], double* objf, double pdG[], Nag_Comm* comm)
{
  globalPolyVolPointer->ParamFitObjFun(iN, pdX, objf, pdG, comm);
}
