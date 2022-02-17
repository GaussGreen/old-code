///////////////////////////////////////////////////////////////////////////////
// File:             common/src/finance/blackscholesutils.cpp                     //
// Purpose:          Useful functions for the resolution of black-scholes PDE// 
// Author:           ZHANG Yunzhi                                            //
// Created:          june 2000                                               //
// RCS-ID:           $Id: blackscholesutils.cpp,v 1.15 2006/06/12 14:48:40 wang Exp $                                                    //
// Copyright         (C) 2000 - 2006 Trilemma LLP                                                //
///////////////////////////////////////////////////////////////////////////////

/**
    @todo Find a way to fix the bounds in functions using newton methods
    @todo Some tests may need to be relative
 */

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"
#include "ito33/debug.h"

#include "ito33/numeric/normaldistribution.h"
#include "ito33/numeric/bisecnewton.h"
#include "ito33/numeric/nonlinearsolver.h"

#include "ito33/finance/blackscholes.h"
#include "ito33/finance/blackscholesutils.h"

using namespace ito33::numeric;

namespace ito33
{

namespace finance
{

  
const double MINVAL = 1.e-8;

//-----------------------------------------------------------------------------

double CalculateCallPrice(double dSpot, double dStrike, double dVolatility,
                          double dMaturity, double dRate, double dForeignRate)

{

  double
    dTmp, 
    dD1, 
    dD2, 
    dResult, 
    dCND1, 
    dCND2;

  ASSERT(dSpot >= 0);

  ASSERT(dStrike > 0);

  if ( fabs(dSpot) < MINVAL )
    return 0.;

  if ( fabs(dStrike) < MINVAL )
    return dSpot;

  if ( dVolatility < MINVAL )
    return 0.;

  dTmp= dVolatility * sqrt(dMaturity);

  dD1 = (log(dSpot / dStrike) + (dRate - dForeignRate) * dMaturity) / dTmp 
        + dTmp * 0.5;

  dD2 = dD1 - dTmp;

  dCND1 = CumulatedNormal(dD1);

  dCND2 = CumulatedNormal(dD2);

  dResult = dSpot * exp(- dForeignRate * dMaturity) * dCND1  
          - dStrike * exp(- dRate * dMaturity) * dCND2;

  return dResult;

}

//-----------------------------------------------------------------------------

double CalculatePutPrice(double dSpot, double dStrike, double dVolatility, 
                         double dMaturity, double dRate, double dForeignRate)

{
  double
    dTmp,
    dMinusD2,
    dMinusD1,
    dCNMinusD1,
    dCNMinusD2,
    dResult;

  ASSERT(dSpot >= 0);

  ASSERT(dStrike > 0);

  if (fabs(dSpot) < MINVAL)
    return dStrike * exp(-dRate * dMaturity);

  if (fabs(dStrike) < MINVAL)
    return 0.;

  if (dVolatility < MINVAL)
    return 0.;

  dTmp = dVolatility * sqrt(dMaturity);

  dMinusD2 = (log(dStrike / dSpot) + (dForeignRate - dRate) * dMaturity) / dTmp 
           + dTmp * 0.5;
 
  dMinusD1 = dMinusD2 - dTmp;

  dCNMinusD1 = CumulatedNormal(dMinusD1);

  dCNMinusD2 = CumulatedNormal(dMinusD2);

  dResult = dStrike * exp(- dRate * dMaturity) * dCNMinusD2 
          - dSpot * exp(- dForeignRate * dMaturity) * dCNMinusD1;

  return dResult;
}

//-----------------------------------------------------------------------------


double CalculateVegaCall(double dSpot, double dStrike, double dMaturity, 
                         double dRate, double dForeignRate, double dVolatility)
{
  double   
    dTmp,
    dSqrtMaturity,
    dD1,
    dND1,
    dResult;

  if (fabs(dSpot) < MINVAL)
    return 0.;

  if (fabs(dStrike) < MINVAL)
    return 0.;

  if (dVolatility < MINVAL)
    return 0.;

  dSqrtMaturity = sqrt(dMaturity);

  dTmp = dVolatility * dSqrtMaturity;

  dD1 = (log(dSpot / dStrike) + (dRate - dForeignRate) * dMaturity) / dTmp
      + dTmp * 0.5;

  dND1 = Normal(dD1);

  dResult = dSpot * exp(-dForeignRate * dMaturity) * dSqrtMaturity * dND1;

  return dResult;
}

//-----------------------------------------------------------------------------
class BSVolFromPrice : public BlackScholes
{

protected :
  double  m_dPrice;

public :

  BSVolFromPrice(double dPrice): m_dPrice(dPrice) { }
    
    
  void operator () (double dVol, double &dPrice, double &dVega)
  {
    dPrice = CalculateCallPrice(m_dSpot, m_dStrike, dVol, m_dMaturity, 
                                m_dRate, m_dForeignRate) - m_dPrice ;

    dVega = CalculateVegaCall(m_dSpot, m_dStrike, m_dMaturity, 
                              m_dRate, m_dForeignRate, dVol);
  }

};


//-----------------------------------------------------------------------------

double CalculateImpliedVolatility(double dSpot, double dStrike, double dPrice, 
                                  double dMaturity, double dRate, 
                                  double dForeignRate)
{
  double
    dResult;

  BSVolFromPrice
    BSValues(dPrice);
    
  BSValues.SetSpotSharePrice(dSpot);
  BSValues.SetStrike(dStrike);
  BSValues.SetForeignRate(dForeignRate);
  BSValues.SetMaturity(dMaturity);
  BSValues.SetRate(dRate);

  ASSERT( dPrice >= 0.);

  ASSERT( dPrice <= dSpot );

  BisecNewton solver;

  //The release version don't stand the lower bound at 0.0
  dResult = solver(BSValues, 0.000000, 10.);

  return dResult; 
}

//-----------------------------------------------------------------------------

double CalculateDelta(int iCallPut, double dSpot, double dStrike, double dRate, 
                      double dForeignRate, double dVolatility, 
                      double dMaturity)
{
  double
    dTmp,
    dSqrtMaturity,
    dD1,
    dND1,
    dDiscount;
  
  if (dVolatility < MINVAL)
    return 0.;

  if (fabs(dSpot) < MINVAL)
    return 0.;

  dDiscount = exp(- dForeignRate * dMaturity);

  if (fabs(dStrike) < 1.e-10)
    return dDiscount;
 
  else
  {
    dSqrtMaturity = sqrt(dMaturity);

    dTmp = dVolatility * dSqrtMaturity;

    dD1 = (log(dSpot / dStrike) + (dRate - dForeignRate) * dMaturity) / dTmp
        + dTmp * 0.5;

    dND1 = CumulatedNormal(dD1);
    
    if (iCallPut == 0)
      dND1 -= 1.0;

    return dDiscount * dND1;
  }
}

//-----------------------------------------------------------------------------

class BSStrikeFromDelta : public BlackScholes
{
protected:
  double  m_dDelta;
  
public:

  BSStrikeFromDelta(double dDelta): m_dDelta(dDelta) {}

  double DeltaDeriv(double dStrike)
  {
    //Compute the derivative of Delta with respect to strike
    //The formula (found with the above black-scholes formula for Delta derived
    //according to Strike with fixed other parameters) is
    // \f$ \frac{e^{-D T} Normal(D1)}{\sigma \sqrt{T} E}
    //where D is the foreignrate, T the relative maturity, d1 the usual one used
    //by wilmott, sigma the volatility, E the strike
    //This function has been validated with the numerical derivative.
    double 
      dND1,
      dD1,
      dTmp,
      dResult;

    if (m_dVolatility < MINVAL)
      return 0.;

    if (m_dSpot < MINVAL)
      return 0.;

    if (dStrike < MINVAL)
      return 0.;

    dTmp = m_dVolatility * sqrt(m_dMaturity);

    dD1 = (log(m_dSpot / dStrike) + m_dMaturity * (m_dRate - m_dForeignRate 
          + 0.5  * m_dVolatility * m_dVolatility))/ dTmp;

    dND1 = Normal(dD1);

    dResult = - dND1 * exp(-m_dForeignRate * m_dMaturity) / (dStrike * dTmp);

    return dResult; 
  }

  void operator () (double dStrike, double &dDelta, double &dDerivDelta)
  {
    dDelta =  CalculateDelta(1, m_dSpot, dStrike, m_dRate, m_dForeignRate, 
                             m_dVolatility, m_dMaturity)-m_dDelta;

    dDerivDelta = DeltaDeriv(dStrike);
  }

};

//-----------------------------------------------------------------------------

double GetStrikeFromDelta(double dMaturity, double dSpot, double dDelta, 
                          double dVolatility, double dRate, 
                          double dForeignRate)
{
  double
    dStrikeMin = 0.0,
    dStrikeMax = dSpot * pow(2.,dMaturity),
    dTolerance = 1e-7,
    dResult;

  BSStrikeFromDelta 
    BSValues(dDelta);
    
  BSValues.SetSpotSharePrice(dSpot);
  BSValues.SetMaturity(dMaturity);
  BSValues.SetRate(dRate);
  BSValues.SetForeignRate(dForeignRate);
  BSValues.SetVolatility(dVolatility);
 
  BisecNewton solver(dTolerance);
  
  dResult = solver(BSValues,dStrikeMin, dStrikeMax);

  return dResult;
}


} // namespace finance

} // namespace ito33
