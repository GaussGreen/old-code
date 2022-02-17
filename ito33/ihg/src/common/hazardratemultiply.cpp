/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/hazardratemultiply.cpp
// Purpose:     implementation of HazardRateSum class
// Author:      Wang
// Created:     2004/06/04
// RCS-ID:      $Id: hazardratemultiply.cpp,v 1.4 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/arraycheckers.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

#include "ito33/ihg/hazardratemultiply.h"

#include "ito33/xml/write.h"
#include "ihg/xml/common.h"
#include "ihg/xml/hazardrate.h"

namespace ito33
{

namespace ihg
{


HazardRateMultiply::HazardRateMultiply
                    (
                      shared_ptr<SpotComponent> &pSpotcomponent,
                      const Date* pDates, 
                      const double *pdValues, 
                      size_t nNbTimes
                    )
  : HazardRateWithTimeComponent(pDates, pdValues, nNbTimes),
    m_pSpotComponent(pSpotcomponent)
{
  CheckNonNegativity(pdValues, nNbTimes);
}

void HazardRateMultiply::GetHazardRates(double dTime, 
                                   const double *pdSpots, 
                                   double *pdValues,
                                   size_t nNbS) const
{
  double dAlpha = (*m_pTimeComponent)(dTime);

  m_pSpotComponent->GetValues(pdSpots, pdValues, nNbS);
  
  for (size_t n = 0; n < nNbS; n++) 
  {
    pdValues[n] *= dAlpha;
  }
}


} // namespace ihg

} // namespace ito33
