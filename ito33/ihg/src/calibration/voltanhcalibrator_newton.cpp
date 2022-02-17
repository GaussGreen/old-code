/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/voltanhcalibrator.cpp
// Purpose:     calibrate tanh volatility
// Created:     2005/02/04
// RCS-ID:      $Id: voltanhcalibrator_newton.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/newton2d.h"
#include "ito33/numeric/numericerror.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/voltanhcalibrator_newton.h"
#include "ihg/impliedvol.h"

extern const ito33::Error ITO33_MAX_ITER;
extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace ihg
{


using numeric::Exception;

void VolTanhCalibratorNewton::operator () 
     (double dLeft, double dRight, double& dPrice1, double& dPrice2,
      double& dF)
{
  m_theoreticalModel.SetHazardRate(m_pHazardRate);
  m_theoreticalModel.SetVolatility(shared_ptr<Volatility>(
      new VolatilityTanh(dLeft, dRight, m_dScale, m_dS0)) );

  dPrice1 = m_theoreticalModel.Compute(*m_pDeriv1)->GetPrice();
  dPrice2 = m_theoreticalModel.Compute(*m_pDeriv2)->GetPrice();

  dPrice1 -= m_dMarketPrice1;
  dPrice2 -= m_dMarketPrice2;

  dF = 0.5 * (  dPrice1 * dPrice1 * m_dScale1 * m_dScale1 
              + dPrice2 * dPrice2 * m_dScale2 * m_dScale2 );

  // Save best values in case calibration fails but the user still wants 
  // the best guess
  if (dF < m_dObjectif)
  {
    m_dObjectif = dF;
    m_dLeft = dLeft;
    m_dRight = dRight;
  }
} 

shared_ptr<VolatilityTanh> 
VolTanhCalibratorNewton::Calibrate
(
  const finance::Derivative& deriv1,
  const finance::Derivative& deriv2,
  const shared_ptr<HazardRate>& pHazardRate,
  size_t nIterMax,
  bool bUsePreviousSolution
)
{ 
  // Init parameters used when the calibration fails
  m_dObjectif = 1.e99;

  // Init everything else
  m_pHazardRate   = pHazardRate;
  m_pDeriv1       = &deriv1;
  m_pDeriv2       = &deriv2;
  m_dMarketPrice1 = m_pDeriv1->GetMarketPrice();;
  m_dMarketPrice2 = m_pDeriv2->GetMarketPrice();

  m_dScale1 = fabs(m_dMarketPrice1);
  if (m_dScale1 > 1.e-6)
    m_dScale1 = 1.0 / m_dScale1;
  else
    m_dScale1 = 1.0;

  m_dScale2 = fabs(m_dMarketPrice2);
  if (m_dScale2 > 1.e-6)
    m_dScale2 = 1.0 / m_dScale2;
  else
    m_dScale2 = 1.0;

  double dRightBound1 = 0.0;
  double dLeftBound1 = 5.0;
  double dRightBound2 = 0.0;
  double dLeftBound2 = 5.0;

  // Check if the caller wants to use the previous solution as initial
  // guess, and this is not the first time the code was called.  
  // Typically, this should be set to true.
  double dLeft;
  double dRight;

  if ( !bUsePreviousSolution || (m_dLeft == 0.0 && m_dRight == 0.0) )
  {
    // Just calibrate the flat vols corresponding to the two prices
    // The first be used to approximate the left limit
    // and the second be used to approximatethe right limit
    ImpliedVol impliedVol1(deriv1, pHazardRate);
    dLeft = impliedVol1.Compute();

    ImpliedVol impliedVol2(deriv2, pHazardRate);
    dRight = impliedVol2.Compute();

    // The implied vol guesses are working for tests in testsuite
    // if newton sets derivatives to zero for values less than 1.e-5

    // Use the best guess from a discrete set of guesses. Also include
    // the implied vol guesses from above
    //double dPrice1, dPrice2, dF;
    //(*this)(dLeft, dRight, dPrice1, dPrice2, dF);
    //FindInitialGuess();
    //dLeft = m_dLeft;
    //dRight = m_dRight;
  } 
  else
  {
    dLeft = m_dLeft;
    dRight = m_dRight;
  }

  // newton solver
  numeric::Newton2D solverN(dRightBound1, dLeftBound1,
                            dRightBound2, dLeftBound2,
                            1.e-8, nIterMax);

  numeric::NumericError err = solverN(*this, dLeft, dRight);

  CHECK_COND(err == numeric::ITO33_NO_ERROR, ITO33_CALIBRATION_FAIL);

  // save the converged result
  m_dLeft = dLeft;
  m_dRight = dRight;

  return shared_ptr<VolatilityTanh>(
            new VolatilityTanh(m_dLeft, m_dRight, m_dScale, m_dS0));  
}

shared_ptr<VolatilityTanh> 
VolTanhCalibratorNewton::GetVolatility()
{
  return shared_ptr<VolatilityTanh>(
            new VolatilityTanh(m_dLeft, m_dRight, m_dScale, m_dS0));  
}

/*
void VolTanhCalibratorNewton::FindInitialGuess()
{
  // Check discrete points in the range small to large 
  double dSmall = 0.01;
  double dLarge = 0.81;
  double dStep = 0.2;

  // Let the () operator track the best guess so far
  double dPrice1, dPrice2, dF;
  for ( double dLeft = dSmall; dLeft < dLarge; dLeft+=dStep)
  {
    for ( double dRight = dSmall; dRight < dLarge; dRight+=dStep)
    {
      (*this)(dLeft, dRight, dPrice1, dPrice2, dF);
    }
  }

}
*/

} // namespace ihg

} // namespace ito33
