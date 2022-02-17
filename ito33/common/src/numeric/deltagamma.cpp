/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/deltagamma.cpp
// Purpose:     Functions for delta and gamma computation
// Author:      Wang
// Created:     2003/12/12
// RCS-ID:      $Id: deltagamma.cpp,v 1.14 2006/01/18 14:09:15 yann Exp $
// Copyright:   (c) 2003 - 2004  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/debug.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/predicatedouble.h"
#include "ito33/numeric/deltagamma.h"

namespace ito33
{

namespace numeric
{


void ComputeDelta(const double *pdS, const double *pdValues, size_t nNbS,
                  double *pdDeltas)
{

  // Handle trivial cases
  ASSERT_MSG(pdS != 0 && pdValues != 0 && pdDeltas != 0 && nNbS > 0,
             "Null pointer or size in ComputeDelta");

  if (nNbS == 1)
  {
    pdDeltas[0] = 0.0;
    return;
  }

  // if only two points linearly interpolate
  // and set the delta of both points to be the same
  if ( nNbS == 2)
  {
    double dTmp =  ( pdValues[1] - pdValues[0] ) / ( pdS[1] - pdS[0] );
    pdDeltas[0] = pdDeltas[1] = dTmp;
    return;
  }

  // The first two points may have been added due to a barrier
  // we force delta to be zero and use forward differencing
  // to compute the value at the dirichlet/barrier node
  size_t nIdxLeft = 0;
  if ( IsEqual( pdValues[0], pdValues[1] ) 
    && IsEqual( pdValues[0], pdValues[2] ) )
  {
    pdDeltas[0] = pdDeltas[1] = 0.;
    nIdxLeft = 2;
  }

  if ( nIdxLeft < nNbS - 1 )
    pdDeltas[nIdxLeft] = ( pdValues[nIdxLeft + 1] - pdValues[nIdxLeft] ) / 
                         ( pdS[nIdxLeft + 1] - pdS[nIdxLeft] );

  // The last two points may have been added due to a barrier
  // we force delta to be zero and use backward differencing
  // to compute the value at the dirichlet/barrier node
  size_t nIdxRight = nNbS - 1;
  if ( IsEqual( pdValues[nNbS - 1], pdValues[nNbS - 2] )
    && IsEqual( pdValues[nNbS - 1], pdValues[nNbS - 3] ) )
  {
    pdDeltas[nNbS - 1] = pdDeltas[nNbS - 2] = 0.;
    nIdxRight = nNbS - 3;
  }
 
  //backward differencing on the last point or dirichlet node
  if (nIdxRight > 0 )
    pdDeltas[nIdxRight] = (pdValues[nIdxRight] - pdValues[nIdxRight - 1]) / 
      (pdS[nIdxRight] - pdS[nIdxRight - 1]);

  // Use central formula at interior points
  for (size_t nIdxS = nIdxLeft + 1; nIdxS < nIdxRight; nIdxS++) 
  {
    double dSi = pdS[nIdxS];
    double dSim = pdS[nIdxS - 1];
    double dSip = pdS[nIdxS + 1];

    double dx_minus = dSi - dSim;
    double dx_plus = dSip - dSi;

    double dValueI = pdValues[nIdxS];
    double dValueIm = pdValues[nIdxS - 1];
    double dValueIp = pdValues[nIdxS + 1];

    pdDeltas[nIdxS] = (  (dValueIp - dValueI) / dx_plus * dx_minus
                       + (dValueI - dValueIm) / dx_minus * dx_plus )
                    / (dx_plus + dx_minus);

  }

}

void ComputeDelta(const double *pdS, const int *piFD, const double *pdValues, 
                  size_t nNbS, 
                  double *pdDeltas)
{
  ASSERT_MSG(pdS != 0 && pdValues != 0 && pdDeltas != 0 && piFD != 0 && nNbS > 0,
             "Null pointer or size in ComputeDelta");

  pdDeltas[0] = (pdValues[1] - pdValues[0]) / (pdS[1] - pdS[0]);

  pdDeltas[nNbS - 1] = (pdValues[nNbS - 1] - pdValues[nNbS - 2])
                     / (pdS[nNbS - 1] - pdS[nNbS - 2]);

  for (size_t nIdxS = 1; nIdxS < nNbS - 1; nIdxS++) 
  {
    double dSi = pdS[nIdxS];
    double dSim = pdS[nIdxS - 1];
    double dSip = pdS[nIdxS + 1];

    double dValueI = pdValues[nIdxS];
    double dValueIm = pdValues[nIdxS - 1];
    double dValueIp = pdValues[nIdxS + 1];
    
    if (piFD[nIdxS] == 1) // central weighting 
      pdDeltas[nIdxS] = (dValueIp - dValueIm) / ( dSip - dSim );
    else if (piFD[nIdxS] == 0) // backward differencing
      pdDeltas[nIdxS] = (dValueI - dValueIm) / (dSi - dSim);
    else // forward differencing
      pdDeltas[nIdxS] = (dValueIp - dValueI) / (dSip - dSi);
  }
}

void ComputeGamma(const double *pdS, const double *pdValues, size_t nNbS, 
                  double *pdGammas)
{
  ASSERT_MSG(pdS != 0 && pdValues != 0 && pdGammas != 0 && nNbS > 0,
             "Null pointer or size in ComputeGamma");

  double dYP0 = (pdValues[1] - pdValues[0]) / (pdS[1] - pdS[0]);

  double dYPNm1 = (pdValues[nNbS - 1] - pdValues[nNbS - 2])
                / (pdS[nNbS - 1] - pdS[nNbS - 2]);

  Spline(pdS, pdValues, nNbS, dYP0, dYPNm1, pdGammas);
}


void ComputeGammaFD(const double *pdS, const double *pdValues, size_t nNbS, 
                    double *pdGammas)
{
  // Check for trivial cases before computing
  ASSERT_MSG(pdS != 0 && pdValues != 0 && pdGammas != 0 && nNbS > 0,
             "Null pointer or size in ComputeGammaFD");

  if (nNbS == 1)
  {
    pdGammas[0] = 0.0;
    return;
  }

  if (nNbS == 2)
  {
    pdGammas[0] = 0.0;
    pdGammas[1] = 0.0;
    return;
  }

  // The first two points may have been added due to a barrier
  // we force gamma to be zero and use extrapolation
  // to compute the value at the dirichlet/barrier node
  size_t nIdxLeft = 0;
  //Check if the first three points have the same value
  if ( IsEqual( pdValues[0], pdValues[1] ) 
    && IsEqual( pdValues[0], pdValues[2] ) )
  {
    pdGammas[0] = pdGammas[1] = 0.;
    nIdxLeft = 2;
  }

  // The first two points may have been added due to a barrier
  // we force gamma to be zero and use extrapolation
  // to compute the value at the dirichlet/barrier node
  size_t nIdxRight = nNbS - 1;
  //check if the last two three points have the same value
  if ( IsEqual( pdValues[nNbS - 1], pdValues[nNbS - 2] )
    && IsEqual( pdValues[nNbS - 1], pdValues[nNbS - 3] ) )
  {
    pdGammas[nNbS - 1] = pdGammas[nNbS - 2] = 0.;
    nIdxRight = nNbS - 3;
  }

 
  // handle Heaveside-like data
  if (nIdxRight == nIdxLeft + 1)
  {
    pdGammas[nIdxLeft] = 0.0;
    pdGammas[nIdxRight] = 0.0;
  }
  
  // Use central formula for the interior points
  for (size_t nIdxS = nIdxLeft + 1; nIdxS < nIdxRight; nIdxS++) 
  {
    double dSi = pdS[nIdxS];
    double dSim = pdS[nIdxS - 1];
    double dSip = pdS[nIdxS + 1];

    double dx_minus = dSi - dSim;
    double dx_plus = dSip - dSi;

    double dValueI = pdValues[nIdxS];
    double dValueIm = pdValues[nIdxS - 1];
    double dValueIp = pdValues[nIdxS + 1];

    pdGammas[nIdxS] = (  (dValueIp - dValueI)/dx_plus 
                       - (dValueI - dValueIm)/dx_minus)
                    / ((dx_plus + dx_minus) * 0.5);
  }

  // Check if we set gamma to zero at boundary points. Also make
  // sure we have enough data points to extrapolate.
  if ( nIdxLeft > 0 && nIdxRight > nIdxLeft + 2)
  {
    double dSlope = (pdGammas[nIdxLeft + 2] - pdGammas[nIdxLeft + 1]) 
                  / (pdS[nIdxLeft + 2] - pdS[nIdxLeft + 1]);
    pdGammas[nIdxLeft] = pdGammas[nIdxLeft + 1] 
                       - dSlope * (pdS[nIdxLeft + 1] - pdS[nIdxLeft]);
  }
  else if ( nIdxLeft < nNbS - 1 )
    pdGammas[nIdxLeft] = pdGammas[nIdxLeft + 1];

  // Check if we set gamma to zero at boundary points. Also make
  // sure we have enough data points to extrapolate.
  if ( nIdxRight < nNbS - 1 && nIdxRight > nIdxLeft + 2)
  {
    double dSlope = ( pdGammas[nIdxRight - 1] - pdGammas[nIdxRight - 2] ) / 
                    ( pdS[nIdxRight - 1] - pdS[nIdxRight - 2] );
    pdGammas[nIdxRight] = pdGammas[nIdxRight - 1] 
                        + dSlope * (pdS[nIdxRight] - pdS[nIdxRight - 1]);
  }
  else if ( nIdxRight > 0 )
    pdGammas[nIdxRight] = pdGammas[nIdxRight - 1];

}

void Compute2nd(const double* pdX, const double* pdPrices, size_t nNbS, 
                double* pdValues)
{
  // we don't care for now the values at boundary
  pdValues[0] = pdValues[nNbS - 1] = 0.;

  for (size_t nIdxS = 1; nIdxS < nNbS - 1; nIdxS++)
  {
    pdValues[nIdxS] = (pdPrices[nIdxS] - pdPrices[nIdxS - 1])
                    / (pdX[nIdxS] - pdX[nIdxS - 1])
                    + (pdPrices[nIdxS] - pdPrices[nIdxS + 1])
                    / (pdX[nIdxS + 1] - pdX[nIdxS]);
  }

}

void Compute1st(const double* /* pdX */, const double* pdPrices, size_t nNbS, 
                double* pdValues)
{
  // we don't care for now the values at boundary
  pdValues[0] = pdValues[nNbS - 1] = 0.;

  for (size_t nIdxS = 1; nIdxS < nNbS - 1; nIdxS++)
  {
    pdValues[nIdxS] = 0.5 * (pdPrices[nIdxS + 1] - pdPrices[nIdxS - 1]);
  }
}


} // namespace numeric

} // namespace ito33

