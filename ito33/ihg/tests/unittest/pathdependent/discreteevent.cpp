///////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/common/pathdependent/discreteevent.cpp
// Purpose:     discrete  event
// Author:      ITO 33 
// Created:     04/02/2005
// RCS-ID:      $Id: discreteevent.cpp,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/binarysearch.h"

#include "discreteevent.h"

#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"

#include "ito33/numeric/interpolation.h"

namespace ito33
{

namespace pricing
{
  
void DiscreteEvent::ApplyAtStartTime(PathDepStructure& /*pathDepStruct*/) const
{

}

void DiscreteEvent::ApplyAtEndTime(PathDepStructure& /*pathDepStruct*/) const
{

}

void DiscreteEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{
   //copy pde1 into pde0 if S is above trigger
  //assume that after each event Vold=Voldold=V
  size_t nGridSize0     = pathDepStruct.m_path[0].meshes->GetNbS();
  const double* pdGrid0 = pathDepStruct.m_path[0].meshes->GetS();
  double* pdValues0     = pathDepStruct.m_path[0].instdata->m_pdPrices.Get();
  
  size_t nGridSize1     = pathDepStruct.m_path[1].meshes->GetNbS();
  const double* pdGrid1 = pathDepStruct.m_path[1].meshes->GetS();
  double* pdValues1     = pathDepStruct.m_path[1].instdata->m_pdPrices.Get();

  Array<double> pdInterpValues(nGridSize0);

  Interpolate(pdGrid1, pdValues1, nGridSize1,
      pdGrid0, pdInterpValues.Get(), nGridSize0,
      numeric::ExtrapolationMode_Linear,
      numeric::ExtrapolationMode_Linear);

  size_t nIdx = 0;

  for (nIdx = 0; nIdx < nGridSize0; nIdx++)
   {
     pdValues0[nIdx] = pdInterpValues[nIdx];
   }
}


} // namespace pricing

} // namespace ito33
