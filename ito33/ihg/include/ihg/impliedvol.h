/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/impliedvol.h
// Purpose:     Class for computing implied volatility
// Author:      David
// Created:     2004/06/01
// RCS-ID:      $Id: impliedvol.h,v 1.20 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/impliedvol.h
    @brief Class for computing implied volatility

    Implemented as a class so that it can interface with the nonlinear
    root finding algorithms
 */

#ifndef _IHG_IMPLIEDVOL_H_
#define _IHG_IMPLIEDVOL_H_

#include <cmath>

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/derivative.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"

namespace ito33 
{
  
namespace ihg 
{


/**
   Class for computing flat implied volatility with known hazard rate.
 */
class ImpliedVol
{

public:

  /**
     Ctor initializes the calibrator with a derivative and a known hazard rate

     @param derivative the derivative to be calibrated
     @param pHazardRate the known hazard rate of the ihg model
     @param pQualityControl the quality control to be used for pricing
   */
  ImpliedVol(const finance::Derivative& derivative,
             const shared_ptr<HazardRate>& pHazardRate,
             const shared_ptr<finance::QualityControl> pQualityControl =
                shared_ptr<finance::QualityControl>())
           : m_derivative(derivative)
  { 
    // Setup a computationalFlags with computeVega activated
    shared_ptr<finance::ComputationalFlags> 
      pComputationalFlags(new finance::ComputationalFlags);

    pComputationalFlags->SetComputeVega(true);

    // Set up the TM to be used
    m_theoreticalModel.SetExternalFlags(pComputationalFlags);
    m_theoreticalModel.SetHazardRate(pHazardRate);
    
    if (pQualityControl)
      m_theoreticalModel.SetQualityControl(pQualityControl);

    // Save market price in case it is costly to compute (eg. if derivative
    // is an option with market price defined by an implied vol)
    m_dMarketPrice = m_derivative.GetMarketPrice();

    // The scale value of the function
    m_dScale = fabs( m_dMarketPrice );

    if ( m_dScale < 1.e-6 )
      m_dScale = 1.0;

    m_dScale = 1. / m_dScale;
  }

  // Default dtor is ok

  /**
      Compute the implied volatility.

      Essentially just sets up and calls a root searching algorithm.

      @return The implied volatility
  */
  double Compute();
  

  /** 
      Used by the non-linear root finding algorithm

      @param dVol The current guess of the root (implied vol)
      @param dError (output) The value of the non-linear function at dVol
      @param dVega (output) The gradient of the non-linear function at dVol
   */
  void operator () (double dVol, double& dError, double& dVega)
  {
    // set the new volatility guess   
    m_theoreticalModel.SetVolatility(shared_ptr<Volatility>(
            new VolatilityFlat(dVol)) );
    
    shared_ptr<finance::ModelOutput> 
      pOutput = m_theoreticalModel.Compute(m_derivative);

    // Return the relative difference
    dError = (pOutput->GetPrice() - m_dMarketPrice) * m_dScale;

    // Set the derivative with respect to volatility to help the
    // root search algorithm
    dVega = pOutput->GetVega() * m_dScale;
  }


protected:

  /// Helper for the root finding
  double m_dScale;

  /// Cached value of the derivative market price
  double m_dMarketPrice;

  /// The theoretical model used for the implied vol
  ihg::TheoreticalModel m_theoreticalModel;

  /// The derivative contract for which the implied vol is sought
  const finance::Derivative& m_derivative;

  NO_COPY_CLASS(ImpliedVol);

}; // class ImpliedVol


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_IMPLIEDVOL_H_
