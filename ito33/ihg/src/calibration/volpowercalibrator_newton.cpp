/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/volpowercalibrator.cpp
// Purpose:     calibrate power volatility
// Author:      Ito33
// Created:     2004/12/27
// RCS-ID:      $Id: volpowercalibrator_newton.cpp,v 1.6 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004-  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/volpowercalibrator.cpp
   @brief calibration of power volatility

 */
#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/useexception.h"
#include "ito33/error.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/newton2d.h"
#include "ito33/numeric/numericerror.h"


#include "ito33/finance/error.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/volpowercalibrator_newton.h"
#include "ihg/impliedvol.h"

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;
extern const ito33::Error ITO33_MAX_ITER;

namespace ito33
{

namespace ihg
{

using numeric::Exception;

void VolPowerCalibratorNewton::operator () 
   (double dAlpha, double dBeta, double &dPrice1, double &dPrice2,
    double &dF)
{

  shared_ptr<Volatility> pVolatility( new VolatilityPower(dAlpha,dBeta,m_dS0) );

  m_theoreticalModel.SetHazardRate(m_pHazardRate);
  m_theoreticalModel.SetVolatility(pVolatility);

  dPrice1 = m_theoreticalModel.Compute(*m_pDeriv1)->GetPrice();
  dPrice2 = m_theoreticalModel.Compute(*m_pDeriv2)->GetPrice();

  dPrice1 -= m_dMarketPrice1;
  dPrice2 -= m_dMarketPrice2;

  dF = .5*( dPrice1*dPrice1*m_dScale1*m_dScale1 
          + dPrice2*dPrice2*m_dScale2*m_dScale2 );

  // Save best values in case calibration fails but the user still wants 
  // the best guess
  if (dF < m_dObjectif)
  {
    m_dObjectif = dF;
    m_dAlpha = dAlpha;
    m_dBeta = dBeta;
  }
} 


shared_ptr<VolatilityPower> 
VolPowerCalibratorNewton::Calibrate
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
  m_dS0           = m_pDeriv1->GetSessionData()->GetSpotSharePrice();
  m_dMarketPrice1 = m_pDeriv1->GetMarketPrice();;
  m_dMarketPrice2 = m_pDeriv2->GetMarketPrice();

  m_dScale1 = fabs(m_dMarketPrice1);
  if (m_dScale1 > 1.e-6)
    m_dScale1 = 1.0/m_dScale1;
  else
    m_dScale1 = 1.0;

  m_dScale2 = fabs(m_dMarketPrice2);
  if (m_dScale2 > 1.e-6)
    m_dScale2 = 1.0/m_dScale2;
  else
    m_dScale2 = 1.0;


  double dAlphaLowerBound = 0.0;
  double dAlphaUpperBound = 5.0;
  double dBetaLowerBound  = -2.0;
  double dBetaUpperBound  = 2.0;

  // Check if the caller wants to use the previous solution as initial
  // guess, and this is not the first time the code was called.  
  // Typically, this should be set to true.  The implied vol
  // is computed and used as the initial guess for alpha.  Initial
  // testing showed that the derivative w.r.t. alpha changes substantially
  // as alpha changes.  This causes problems for Newton.
  double dAlpha;
  double dBeta;
  if ( !bUsePreviousSolution || (m_dAlpha == 0.0 && m_dBeta == 0.0) )
  {
    // Assume the first deriv is closer to being at the money
    // If this fails, then let it throw an exception since there is no way
    // a vol power can be calibrated
    ImpliedVol impliedVol1(deriv1, pHazardRate);
    double dAlpha1 = impliedVol1.Compute();

    dAlpha = dAlpha1;
    dBeta  = 0.0;
  } 
  else
  {
    dAlpha = m_dAlpha;
    dBeta = m_dBeta;
  }

  //newton solver
  numeric::Newton2D solverN(dAlphaLowerBound, dAlphaUpperBound,
                           dBetaLowerBound, dBetaUpperBound,
                           1.e-8, nIterMax);

  numeric::NumericError err = solverN(*this, dAlpha, dBeta);
/*  
  if ( err != numeric::ITO33_CONVERGED )
  {
    // Try again with different initial guess
    dAlpha = 0.0;
    dBeta = 0.0;
    err = solverN(*this, dAlpha, dBeta);
  }
*/
  if ( err == numeric::ITO33_NOT_CONVERGED )
      throw EXCEPTION_MSG
            (
              ITO33_CALIBRATION_FAIL,
              TRANS("Can not calibrate the power volatility.")
            );

  if ( err == numeric::ITO33_TOO_MANY_ITERATION )
      throw EXCEPTION_MSG
            (
              ITO33_MAX_ITER,
              TRANS("Too many iterations.")
            );

  // save the converged result
  m_dAlpha = dAlpha;
  m_dBeta = dBeta;

  shared_ptr<VolatilityPower> 
    pVol(new VolatilityPower(m_dAlpha, m_dBeta, m_dS0));

  return pVol;  
}


shared_ptr<VolatilityPower> 
VolPowerCalibratorNewton::GetVolatility()
{

  shared_ptr<VolatilityPower> 
    pVol(new VolatilityPower(m_dAlpha, m_dBeta, m_dS0));

  return pVol;  
}


} // namespace ihg

} // namespace ito33
