/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/conditionalswapevent.cpp
// Purpose:     conditional variance swap event
// Created:     2006/07/14
// RCS-ID:      $Id: conditionalswapevent.cpp,v 1.3 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/binarysearch.h"

#include "ito33/pricing/conditionalswapevent.h"
#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"

#include "ito33/numeric/predicatedouble.h"

namespace ito33
{

namespace pricing
{
  

void ConditionalSwapEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{
  size_t nPath = pathDepStruct.m_path.size();
  size_t nIdPath;

  // Updating in place is not possible (subsequent updates may need 
  // previously altered data). Create a new temporary solution matrix.
  std::vector< std::vector<double> > V(nPath);
  
  for ( nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    size_t nGridSize = pathDepStruct.m_path[nIdPath].meshes->GetNbS();

    V[nIdPath].resize(nGridSize);        
  } // end loop over paths

  // Get the corridor values. If not set, then make the corridor
  // wider than the grid so that the problem reduces to the usual
  // variance swap
  double dUpCorridorBarrier = m_varianceSwap.GetUpCorridorBarrier();
  double dDownCorridorBarrier = m_varianceSwap.GetDownCorridorBarrier();

  if (dDownCorridorBarrier == 0.0)
    dDownCorridorBarrier = 1.e99;

  // loop over each path. The top path never changes.
  for (nIdPath = 0; nIdPath < nPath-1; nIdPath++)  
  { 
    // Get the mesh of the current path
    size_t nGridSize  = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    const double* pdS = pathDepStruct.m_path[nIdPath].meshes->GetS();
    const double* pdValues = pathDepStruct.GetPriceData(nIdPath);

    for (size_t nIdS = 0; nIdS < nGridSize; nIdS++)
    { 
      double dS = pdS[nIdS];

      // If within the corridor, get the new value from the path with
      // one higher setting of Y.  Otherwise, no change in the value.
      if ( dS > dUpCorridorBarrier && dS < dDownCorridorBarrier )
        V[nIdPath][nIdS] = pathDepStruct.QuadraticInterpolate(nIdPath+1, dS);
      else
        V[nIdPath][nIdS] = pdValues[nIdS];

    } // end loop grid

  } // end loop path
  
  // Copy the solution back into the path structure
  for ( nIdPath = 0 ; nIdPath < nPath-1 ; nIdPath++)
  {
    size_t nGridSize= pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    double* pdValues= pathDepStruct.GetPriceData(nIdPath);

    for (size_t nIdS = 0 ; nIdS < nGridSize; nIdS++)
      pdValues[nIdS] = V[nIdPath][nIdS]; 
  
  } //end loop path
  
   // last observation turn off all paths
  if ( m_bRecurse == false && m_bIsLastEvent == true )
  {
 
    // To accurately compute theta, we need correct values at the previous
    // timestep.  Thus, also apply the event to the previous solution
    for (nIdPath = 0; nIdPath < nPath; nIdPath++)
    {
      swap(pathDepStruct.m_path[nIdPath].instdata->m_pdPrices,
           pathDepStruct.m_path[nIdPath].instdata->m_pdOldPrices);
    }

    // Set the recurse flag to avoid endless recursion
    m_bRecurse = true;
    this->ApplyAtTime(pathDepStruct);
    m_bRecurse = false;

    for (nIdPath = 0; nIdPath < nPath; nIdPath++)
    {
      swap(pathDepStruct.m_path[nIdPath].instdata->m_pdPrices,
           pathDepStruct.m_path[nIdPath].instdata->m_pdOldPrices);
    }
  }
 
} // ApplyAtTime(PathDepStructure& pathDepStruct) const

} // namespace pricing

} // namespace ito33
