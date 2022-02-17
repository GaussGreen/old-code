/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/checkfit.h
// Purpose:     util function for tests
// Author:      Yann d'Halluin 
// Created:     2004/06/13
// RCS-ID:      $Id: checkfit.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/utiltest.h
    @brief set of functions to help when writting tests results
**/

#ifndef _IHG_INCLUDE_TESTS_CHECKFIT_H_
#define _IHG_INCLUDE_TESTS_CHECKFIT_H_

#include "ito33/cppunit.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/termstructure.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace finance
{
  class Derivative;
  class Derivatives;
  class ModelOutput;
  template<class T> class TermStructure;
}
  
namespace ihg
{

class Volatility;
class HazardRate;
class TheorecticalModel;

template<class T>
void CheckFit(
  const ito33::shared_ptr<HazardRate>& pHR,
  const ito33::shared_ptr<Volatility>& pVol,
  const ito33::finance::TermStructure<T>& tsDerivative)
{
  // Create and setup the pricing model
  TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  // compare the prices
  const TermStructure<T>::Elements& pCDSElements = tsDerivative.GetAll();
  TermStructure<T>::Elements::const_iterator ppCDS;
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const Derivative& deriv = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(deriv);
    
    double dScale = fabs( deriv.GetMarketPrice() );
    if (dScale < 1.e-6)
      dScale = 1.0;

    double dError = fabs( output->GetPrice() - deriv.GetMarketPrice() );
    CPPUNIT_ASSERT( dError/dScale < 1.e-3);
  }
}

inline void CheckFit(
  const ito33::shared_ptr<HazardRate>& pHR,
  const ito33::shared_ptr<Volatility>& pVol,
  const finance::Derivative& deriv)
{
  // Create and setup the pricing model
  TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  shared_ptr<finance::ModelOutput> output = model.Compute(deriv);
    
  double dScale = fabs( deriv.GetMarketPrice() );
  if (dScale < 1.e-6)
    dScale = 1.0;

  double dError = fabs( output->GetPrice() - deriv.GetMarketPrice() );
  CPPUNIT_ASSERT( dError/dScale < 1.e-3);
}

inline void CheckFit(
  const ito33::shared_ptr<HazardRate>& pHR,
  const ito33::shared_ptr<Volatility>& pVol,
  const finance::Derivatives& derivs)
{
  // Create and setup the pricing model
  TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  // loop over the derivatives
  finance::Derivatives::Elements::const_iterator iter;
  const finance::Derivatives::Elements& elements = derivs.GetAll();

  for (iter = elements.begin(); iter != elements.end(); ++iter)
  {
    shared_ptr<finance::ModelOutput> output = model.Compute( *iter->first );
    
    double dScale = fabs( iter->first->GetMarketPrice() );
    if (dScale < 1.e-6)
      dScale = 1.0;

    double dError = fabs( output->GetPrice() - iter->first->GetMarketPrice() );
    CPPUNIT_ASSERT( dError/dScale < 1.e-3);
  }
}


} // namespace ihg

} // namespace ito33


#endif
