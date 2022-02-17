/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/varianceswapevent_special.cpp
// Purpose:     In the case of a variance swap with no cap
//              where the volatility and hazard rate are not
//              spot dependent it is then possible to solve
//              a variance swap using only two paths.
// Created:     June 8, 2006
// RCS-ID:      $Id: varianceswapevent_special.cpp,v 1.4 2006/07/25 20:09:42 dave Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/binarysearch.h"

#include "ito33/pricing/varianceswapevent_special.h"
#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"


namespace ito33
{

namespace pricing
{
  

void VarianceSwapEventSpecial::DoEvent(PathDepStructure& pathDepStruct) const
{
  size_t nPath = pathDepStruct.m_path.size();
  size_t nIdPath;

  // Get references to the average squared return grid and previous spot grid.
  // They are needed in function calls below.
  const std::vector<double>& pdAvgSqrReturnGrid = 
    pathDepStruct.m_pppdGrids[0][0];
  const std::vector<double>& pdPreviousGrid = pathDepStruct.m_pppdGrids[0][1];
  size_t nNbPrevious = pdPreviousGrid.size();

  ASSERT_MSG( pdAvgSqrReturnGrid.size() == 2, 
    "VarianceSwapEventSpecial: The number of Z paths must be equal to 2.");

  // since updating in place is not possible (subsequent updates may need 
  // previously altered data), create a new solution matrix ordered 
  // as [return * previous][S] 
  std::vector< std::vector<double> > V(nPath);

  // In the case of variance swap, 
  // V(S,P,Z,t) = V(S/P Pstar, t; PStar, Z)
  // However the updating rules says S^+ = P^+
  // So we can save time by saving the interpolation value
  // at Pstar for each path. See documentation for more
  // details.
  m_pdSimilarityValues.resize(nPath, 0.0);

  for ( nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    size_t nGridSize = pathDepStruct.m_path[nIdPath].meshes->GetNbS();

    V[nIdPath].resize(nGridSize);        
    
    size_t nIdxP = nIdPath % nNbPrevious;
    m_pdSimilarityValues[nIdPath] 
       = pathDepStruct.QuadraticInterpolate(nIdPath, pdPreviousGrid[nIdxP] );

  } // end loop over paths

  // Check if corridors will limit the returns
  double dUpCorridorBarrier = m_varianceSwap.GetUpCorridorBarrier();
  double dDownCorridorBarrier = m_varianceSwap.GetDownCorridorBarrier();

  // loop over each path. Interpolate to get the new price at each mesh point
  for (nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    // Get the state variables of the current path. Could cast the params
    // and get from varianceswap params, but the construction of paths
    // is internal, so this method should be safe
    size_t nIdxReturn = nIdPath / nNbPrevious;
    size_t nIdxPrevious = nIdPath % nNbPrevious;

    double dAvgSqrReturn = pdAvgSqrReturnGrid[nIdxReturn];
    double dPrevious = pdPreviousGrid[nIdxPrevious];

    // Get the mesh of the current path
    size_t nGridSize  = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    const double* pdS = pathDepStruct.m_path[nIdPath].meshes->GetS();

    // Check if outside a corridor
    bool bComputeReturn = true;
    if ( dUpCorridorBarrier > 0.0 && dPrevious <= dUpCorridorBarrier )
      bComputeReturn = false;    
    if ( dDownCorridorBarrier > 0.0 && dPrevious >= dDownCorridorBarrier )
      bComputeReturn = false;

    for (size_t nIdS = 0; nIdS < nGridSize; nIdS++)
    { 
      double dS = pdS[nIdS];
      
      double dPreviousNew = dS;

      // Only compute returns if within a corridor
      double dR = 0.0;
      if ( bComputeReturn )
      {
        if (m_returnType == finance::Return_Actual)
          dR = (dS - dPrevious) / dPrevious;
        else
          dR = log(dS / dPrevious);
      }

      double dF = dR*dR;

      if ( m_varianceSwap.GetSwapPayoffType() == finance::SwapPayoff_Gamma )
        dF *= dS / m_varianceSwap.GetStartSharePrice(); 

      double dAvgSqrReturnNew = 
        dAvgSqrReturn + (dF - dAvgSqrReturn) / double(m_nObservation);

      V[nIdPath][nIdS] = Interpolate(dPreviousNew, dAvgSqrReturnNew,
        pdPreviousGrid, pdAvgSqrReturnGrid);

      
    } // end loop grid

  } // end loop path
  
  // Copy the solution back into the path structure
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {
    size_t nGridSize= pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    double* pdValues= pathDepStruct.GetPriceData(nIdPath);

    for (size_t nIdS = 0 ; nIdS < nGridSize; nIdS++)
      pdValues[nIdS] = V[nIdPath][nIdS];
  
  } //end loop path
 
} // ApplyAtTime(PathDepStructure& pathDepStruct) const


double VarianceSwapEventSpecial::Interpolate(
          double dP, double dZ, 
          const std::vector<double>& pdPreviousGrid,
          const std::vector<double> &pdAvgSqrReturnGrid) const
{
  size_t nNbPrevious = pdPreviousGrid.size();

  // Find closest paths in Z direction. 2nd upstream point found later.
  size_t nZUp  = 0;
  size_t nZLow = 1;
 
  double dValLow = 0.0;
  double dValUp  = 0.0;

  if ( m_bHasSimilarityReduction )
  {     
    dValLow = m_pdSimilarityValues[nZLow]; 
    dValUp = m_pdSimilarityValues[nZUp];
  }
  else
  {
    // Find closest paths in previous direction.  The diagonal interpolation
    // only uses paths above and below the target P value.
    size_t nIdxPrevious2 = BinSearch(&pdPreviousGrid[0], nNbPrevious, dP);
    size_t nIdxPrevious1 = nIdxPrevious2 - 1;

    size_t nIdx11 = nZLow * nNbPrevious + nIdxPrevious1;    
    size_t nIdx12 = nZLow * nNbPrevious + nIdxPrevious2;
  
    size_t nIdx21 = nZUp * nNbPrevious + nIdxPrevious1;
    size_t nIdx22 = nZUp * nNbPrevious + nIdxPrevious2;    

    // Linearly interpolate the similarity values
    double dVal11 = m_pdSimilarityValues[nIdx11];
    double dVal12 = m_pdSimilarityValues[nIdx12];
    double dVal21 = m_pdSimilarityValues[nIdx21];
    double dVal22 = m_pdSimilarityValues[nIdx22];
    
    double dDist = pdPreviousGrid[nIdxPrevious2] 
                    - pdPreviousGrid[nIdxPrevious1];
    double dWeight = (dP - pdPreviousGrid[nIdxPrevious1]) / dDist;

    dValLow = (1.0 - dWeight) * dVal11 + dWeight * dVal12;
    dValUp  = (1.0 - dWeight) * dVal21 + dWeight * dVal22;
  }

  // Now interpolate in Z direction
  // Only linear in Z direction
  double dZLow = pdAvgSqrReturnGrid[ nZLow ];
  double dZUp  = pdAvgSqrReturnGrid[ nZUp ];

  double dWeight = (dZ - dZLow) / (dZUp - dZLow);

  double dTmp = (1.0 - dWeight) * dValLow  + dWeight * dValUp;

  return dTmp;
}

} // namespace pricing

} // namespace ito33
