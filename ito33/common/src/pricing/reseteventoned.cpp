/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/reseteventoned.cpp
// Purpose:     implementation of the reset events
// Author:      Yann d'Halluin David Pooley
// Created:     2004/10/25
// RCS-ID:      $Id: reseteventoned.cpp,v 1.2 2004/11/12 17:02:19 zhang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/array.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/pricing/reseteventoned.h"



using namespace ito33::numeric;
using ito33::pricing::ResetEvent1D;


/// Function to apply the event to prices pdValues on grid pdS
void ResetEvent1D::ApplyToPrice(
                              const double* pdS, 
                              double* pdValues, 
                              size_t nNbS) const
{
 
  // Create temporary arrays for interpolation
  Array<double> pdSTmp(nNbS);
  Array<double> pdValuesTmp(nNbS);

  // Shift the S mesh by the amount of the dividend
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbS; nIdx++)
  {
    pdSTmp[nIdx]
      = std::max
              ( 
                m_dConversionRatioFloorRate*pdS[nIdx],
                std::min(m_dConversionRatioCapRate*pdS[nIdx],
                         (m_dNominal / m_dMultiplier))
              );
  }
  //Now interpolate from (pdS,pdValues) to (pdSTmp, pdValuesTmp)  
  
  /*
  Interpolate(pdS, pdValues, nNbS, pdSTmp.Get(), pdValuesTmp.Get(), nNbS, 
              ExtrapolationMode_Constant, ExtrapolationMode_Linear);  
  */
  
  QuadraticInterpolate(pdS, pdValues, nNbS, 
                       pdSTmp.Get(), pdValuesTmp.Get(), nNbS, 
                       ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  // Copy the interpolated values back into the original value array
  for (nIdx = 0; nIdx < nNbS; nIdx++)
    pdValues[nIdx] = pdValuesTmp[nIdx];
    
}

