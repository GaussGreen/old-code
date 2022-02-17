///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/normaldistribution.cpp
// Purpose:          useful functions related to normal distributions
// Author:           laurence
// Created:          18/09/2003
// RCS-ID:           $Id: normaldistribution.cpp,v 1.6 2004/10/05 09:13:46 pedro Exp $
// Copyright:        (c) 2003 Trilemma LLP all rights reserved
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Functions previously implemeted in the "common/src/utile.cpp" source file
///////////////////////////////////////////////////////////////////////////////

/**

  @file ito33/numeric/normaldistribution.cpp

  @brief Useful functions related to normal distributions

 */


// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/numeric/normaldistribution.h"

namespace ito33{

namespace numeric{

//-----------------------------------------------------------------------------
// Normal law density probability
//-----------------------------------------------------------------------------

/**

  @brief calculate \f$ e^{-x^2/2}/ \sqrt{2 \pi} \f$
  
 */
double Normal(double dX) 

{
  // @todo use a better formula when fabs(dX) is large
  return exp(- dX * dX * 0.5) * INVSQRT2PI ;

}

//-----------------------------------------------------------------------------
// Cumulated normal distribution
//-----------------------------------------------------------------------------


/**
                                                
  @brief Calculate the cumulated normal 
  \f$ \int_{-\infty}^x e^{-t^2/2}/\sqrt{2 \pi}\f$
  using the formulae found in "handbook of mathematical functions" by Milton
  Abramowitz and Irene A. Stegun, p931-932, according to expressions 26.2.2
  and 26.2.17

  @todo the precision of the CumulatedNormal is only of 1e-7. But using it 
  to compute the BS price of a call with a spot around 100 make the precision
  fall down to 1e-4 for the price. We shall try to use a more precise
  method (power series have been tried but are not very reliable)

 */
double CumulatedNormal(double dX)

{
  double
    dR,
    dT,
    dResult,
    dTmp;

  dT = 1.0 / (1.0 + 0.2316419 * (dX < 0.0 ? - dX : dX));

  dTmp = ((((1.330274429 * dT - 1.821255978) * dT + 1.781477937) * dT  
          - 0.356563782) * dT + 0.31938153) * dT;

  dR = 1.0 - exp(- dX * dX * 0.5) * INVSQRT2PI * dTmp;


  dResult = dX < 0.0 ? 1.0 - dR : dR;

  return dResult;

}

} // namespace numeric

} // namespace ito33
