///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/gammadistribution.cpp
// Purpose:          useful functions related to gamma distributions
// Author:           Ito33 Canada
// Created:          April 24, 2006
// RCS-ID:           $Id: gammadistribution.cpp,v 1.3 2006/05/16 21:05:24 wang Exp $
// Copyright:        (c) 2006 -  Trilemma LLP all rights reserved
///////////////////////////////////////////////////////////////////////////////


/**

  @file ito33/numeric/gammadistribution.cpp

  @brief Useful functions related to gamma distributions

 */

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"
#include "ito33/debug.h"

#include "ito33/numeric/gammadistribution.h"
#include "ito33/numeric/bisecnewton.h"
#include "ito33/numeric/predicatedouble.h"

namespace ito33
{

namespace numeric
{


double GammaLn(double dX)
{
  ASSERT_MSG(dX > 0, 
    "GammaLn: dX must be striclty positive.");
  
  double dY, dTmp, dSer;

  static double pdCoef[6]={76.18009172947146, -86.50532032941677,
                           24.01409824083091, -1.231739572450155,
                           0.1208650973866179e-2, -0.5395239384953e-5};  
  dY = dX;

  dTmp = dX + 5.5;

  dTmp -= (dX + 0.5) * log(dTmp);

  dSer = 1.000000000190015;

  for ( size_t nIdx = 0; nIdx <= 5; nIdx++) 
    dSer += pdCoef[nIdx]/++dY;

  return -dTmp + log(2.5066282746310005 * dSer/dX);
}



double  GammaSeries(double dA, double dX)
{
  if ( numeric::IsEqual( dX, 0.0 ) )
    return dX;

  ASSERT_MSG( dX > 0, 
    "GammaSeries: dX must be strictly positive.");

  double dGammaLn = GammaLn( dA );
  double dAp = dA;
  double dDel, dSum;
  dDel = dSum = 1./dA;

  size_t nIdx;

  for (nIdx = 1; nIdx <= 100; nIdx++) 
  {
    ++dAp;
    dDel *= dX / dAp;
    dSum += dDel;
  
    if ( fabs( dDel) < fabs( dSum ) * 1.e-6) 
      return dSum * exp( -dX + dA * log(dX) - dGammaLn);
  }

  ASSERT_MSG(nIdx <= 100, "GammaSeries: Too many iterations.");

  return 0.0;
}

double IncompleteGammaP(double dA, double dX)
{
  if ( dX < ( dA + 1.0)) 
    return GammaSeries(dA, dX);

  return 1. - GammaContinuousFraction(dA, dX);
}

double IncompleteGammaQ(double dA, double dX)
{
  if ( dX < ( dA + 1.0)) 
    return 1. - GammaSeries(dA, dX);

  return GammaContinuousFraction(dA, dX);
}

double Chi2Cdf(double dA, double dX)
{
  return IncompleteGammaP(.5 * dA , .5 * dX );
}

  
class InverseChiSquaredCDF
{
protected: 
  double m_dA;

  double m_dY;


public:
   
  InverseChiSquaredCDF(double dA, double dY): m_dA(dA), m_dY(dY)  {}
    
  void operator() (double dX, double &dFValue, double &dDerivValue)  
  {
    dFValue     = Chi2Cdf(m_dA, dX);
    dDerivValue = ( Chi2Cdf(m_dA, dX + 1.e-6) - dFValue )/ 1.e-6;
    dFValue    -= m_dY;
  }

};

double InvChi2Cdf(double dA, double dX)
{
  InverseChiSquaredCDF invcdf(dA, dX);

  BisecNewton solver;

  double dResult = solver(invcdf, 0.0, 400.);

  return dResult;
}

double GammaContinuousFraction(double dA, double dX)
{
  double dB = dX + 1.0 - dA;
  double dC = 1.0/1.e-30;
  double dD = 1.0/dB;
  double dH = dD;

  size_t nIdx;

  for (nIdx = 1 ; nIdx <= 100 ; nIdx++) 
  { 
    double dAn = - (nIdx - dA) * nIdx;
    dB += 2.0;
    dD = dAn * dD + dB;

    if ( fabs(dD) < 1.e-30) 
      dD = 1.e-30; 
    
    dC = dB + dAn / dC;

    if (fabs(dC) < 1.e-30) 
      dC = 1.e-30;

    dD = 1.0/dD;

    double dDel = dD * dC;
    
    dH *= dDel;
    
    if (fabs( dDel - 1.0) < 1.e-7) 
      break;
  }

  ASSERT_MSG(nIdx <= 100, "Gamma continuous fraction: too many iterations.");

  return exp( -dX + dA * log( dX )- GammaLn(dA) )* dH; 
}



} // namespace numeric

} // namespace ito33
