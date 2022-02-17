/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatility.cpp
// Purpose:     implementations of volatility class
// Author:      z
// Created:     03.09.23
// RCS-ID:      $Id: volatility.cpp,v 1.20 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/common/volatility.cpp
 */

#include "ito33/sharedptr.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/ihg/volatility.h"
#include "ihg/volatilityperturb.h"

namespace ito33
{

namespace ihg
{


void Volatility::GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const
{
  specialTimes.clear();
}

void Volatility::GetVolsSquared(double dTime, 
                                const double *pdS, 
                                double *pdVolsSquared,
                                size_t nNbS) const
{
  GetVols(dTime, pdS, pdVolsSquared, nNbS);
  size_t n;

  for ( n = 0; n < nNbS; n++)
    pdVolsSquared[n] *= pdVolsSquared[n];
}

shared_ptr<Volatility> Volatility::Perturb(double dShift) const
{
  return shared_ptr<Volatility>( new VolatilityPerturb(*this, dShift) );
}

// empty implementation so that derived classes are not forced to implement it
void 
Volatility::GetModelParameters(finance::ModelParametersConsumer& /*visitor*/) const
{
}


} // namespace ihg

} // namespace ito33
