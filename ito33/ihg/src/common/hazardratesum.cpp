/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/hazardratesum.cpp
// Purpose:     implementation of HazardRateSum class
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: hazardratesum.cpp,v 1.11 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

#include "ito33/ihg/hazardratesum.h"

#include "ito33/xml/write.h"
#include "ihg/xml/common.h"
#include "ihg/xml/hazardrate.h"

namespace ito33
{

namespace ihg
{

void HazardRateSum::GetHazardRates(double dTime, 
                                   const double *pdSpots, 
                                   double *pdValues,
                                   size_t nNbS) const
{
  double dAlpha = (*m_pTimeComponent)(dTime);

  m_pSpotComponent->GetValues(pdSpots, pdValues, nNbS);
  
  for (size_t n = 0; n < nNbS; n++) 
  {
    pdValues[n] += dAlpha;
  }
}


} // namespace ihg

} // namespace ito33
