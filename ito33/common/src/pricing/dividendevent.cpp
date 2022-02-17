/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/dividendevent.cpp
// Purpose:     implementation of the base dividend event
// Created:     2005/06/02
// RCS-ID:      $Id: dividendevent.cpp,v 1.7 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
  @todo: In case of Dirichlet boundary condition, the extrapolation should 
         really use constant mode. This means that for now EDS with dividend
         is priced wrongly. 
*/

#include "ito33/array.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/pricing/dividendevent.h"

namespace ito33
{

namespace pricing
{
   
  using namespace numeric;


InterpolationMatrix* 
DividendEvent::GetInterpolationMatrix
(const double* pdS, size_t nNbS, size_t nNbSubSystem) const
{
  // Create temporary arrays for interpolation
  Array<double> pdSTmp(nNbS);

  ApplyToSpots(pdS, pdSTmp.Get(), nNbS);

  return new InterpolationMatrix(pdSTmp.Get(), nNbS, pdS, nNbS, nNbSubSystem, 
                            true, m_emLeft, m_emRight, m_interpolationMethod);
}

void 
DividendEvent::ApplyToPrice
(const double* pdS, double* pdValues, size_t nNbS) const
{
  // Create temporary arrays for interpolation
  Array<double> pdSTmp(nNbS);
  Array<double> pdValuesTmp(nNbS);

  ApplyToSpots(pdS, pdSTmp.Get(), nNbS);

  // Now interpolate from (pdS,pdValues) to (pdSTmp, pdValuesTmp)  
   
  Interpolate(pdS, pdValues, nNbS, pdSTmp.Get(), pdValuesTmp.Get(), nNbS, 
              m_emLeft, m_emRight, m_interpolationMethod);  

  // Copy the interpolated values back into the original value array
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdValues[nIdx] = pdValuesTmp[nIdx];
}


} // namespace pricing

} // namespace ito33
