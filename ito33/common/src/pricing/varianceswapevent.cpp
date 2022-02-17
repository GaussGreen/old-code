/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/varianceswapevent.cpp
// Purpose:     variance swap event
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapevent.cpp,v 1.11 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/binarysearch.h"

#include "ito33/pricing/varianceswapevent.h"
#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"

#include "ito33/numeric/predicatedouble.h"
#include "ito33/numeric/interpolation.h"

namespace ito33
{

namespace pricing
{
  

void VarianceSwapEvent::DoEvent(PathDepStructure& pathDepStruct) const
{
  size_t nPath = pathDepStruct.m_path.size();
  size_t nIdPath;

  // Get references to the average squared return grid and previous spot grid.
  // They are needed in function calls below.
  const std::vector<double>& pdAvgSqrReturnGrid = 
    pathDepStruct.m_pppdGrids[0][0];
  const std::vector<double>& pdPreviousGrid = pathDepStruct.m_pppdGrids[0][1];
  size_t nNbPrevious = pdPreviousGrid.size();
  
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

      // Compute the new state variables
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

      //V[nIdPath][nIdS] = LinearInterpolate(dS, dPreviousNew, dAvgSqrReturnNew,
      //  pdPreviousGrid, pdAvgSqrReturnGrid, pathDepStruct);

      //V[nIdPath][nIdS] = QuadraticInterpolate(dS, dPreviousNew, dAvgSqrReturnNew,
      //  dR * dR, pdPreviousGrid, pdAvgSqrReturnGrid, pathDepStruct);

      V[nIdPath][nIdS] = QuadZDiagPInterpolate(dPreviousNew, dAvgSqrReturnNew,
        dF, pdPreviousGrid, pdAvgSqrReturnGrid);

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
}


void VarianceSwapEvent::CopyFromDiagonal(PathDepStructure& pathDepStruct) const
{
  // Get the S grid and values for the path to be saved
  size_t nPathToSave = pathDepStruct.GetPathToSave();
  size_t nGridSize   = pathDepStruct.m_path[nPathToSave].meshes->GetNbS();
  double* pdValues   = pathDepStruct.GetPriceData(nPathToSave);
  const double* pdS  = pathDepStruct.m_path[nPathToSave].meshes->GetS();

  // Get reference to the previous spot grid.
  const std::vector<double>& pdPreviousGrid = pathDepStruct.m_pppdGrids[0][1];
  size_t nNbPrevious = pdPreviousGrid.size();

  if ( nNbPrevious == 1 ) 
  {
    // With a similarity reduction in P, the price is constant
    double dSpotSharePrice 
      = pathDepStruct.m_path[nPathToSave].params->GetSpotSharePrice();

    double dVal = pathDepStruct.QuadraticInterpolate(nPathToSave, 
      dSpotSharePrice );

    for ( size_t nIdS = 0; nIdS < nGridSize ; nIdS++)  
      pdValues[nIdS] = dVal;      
  }
  else
  {
    // Interpolate along the diagonal S=P for the plane where Z
    // corresponds to the path to save.
    // The P grid may be different from the S grid, so store the data
    // in a temp array, then intepolate onto the S grid.
    size_t nIdxReturn = nPathToSave / nNbPrevious;
    size_t nOffSet = nIdxReturn * nNbPrevious;
    std::vector<double> pdTmpValues(nNbPrevious);

    for ( size_t nIdx = 0; nIdx < nNbPrevious; nIdx++ )  
    {    
      double dPrevious = pdPreviousGrid[nIdx];

      pdTmpValues[nIdx] = 
        pathDepStruct.QuadraticInterpolate(nOffSet + nIdx, dPrevious );
    } // end loop over paths

    // Interpoate onto the S grid
    numeric::QuadraticInterpolate(&pdPreviousGrid[0], &pdTmpValues[0], 
                                  nNbPrevious,
                                  pdS, pdValues, nGridSize);
  } //end if
  
}

void VarianceSwapEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{

  if ( m_nObservation > 0 )
    DoEvent(pathDepStruct);
  else
  {
    // At start of sampling date.  Copy from diagonal.
    CopyFromDiagonal(pathDepStruct);
  }

   // support for theta
  if ( m_bRecurse == false &&
      (m_nObservation == 0 || m_bIsLastEvent == true) )
  {
    size_t nIdPath;
    size_t nPath = pathDepStruct.m_path.size();

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


double VarianceSwapEvent::LinearInterpolate(double dS, double dP, 
    double dZ, const std::vector<double>& pdPreviousGrid,
    const std::vector<double> &pdAvgSqrReturnGrid,  
    PathDepStructure& pathDepStruct) const
{
  size_t nNbReturn = pdAvgSqrReturnGrid.size();
  size_t nNbPrevious = pdPreviousGrid.size();

  // Find closest paths in average sqr return (Z) direction
  size_t nIdxReturn2 = BinSearch(&pdAvgSqrReturnGrid[0], nNbReturn, dZ);
  size_t nIdxReturn1 = nIdxReturn2 - 1;

  if ( m_bHasSimilarityReduction )
  { 
    double dValLow = m_pdSimilarityValues[nIdxReturn1];

    double dValHigh = m_pdSimilarityValues[nIdxReturn2]; 

    double dWeight = (dZ - pdAvgSqrReturnGrid[nIdxReturn1])
      / (pdAvgSqrReturnGrid[nIdxReturn2] - pdAvgSqrReturnGrid[nIdxReturn1]);

    double dVal = (1.0 - dWeight) * dValLow + dWeight * dValHigh;

    return dVal;
  } 

  // Find closest paths in average sqr return (Z) direction
  size_t nIdxPrevious2 = BinSearch(&pdPreviousGrid[0], nNbPrevious, dP);
  size_t nIdxPrevious1 = nIdxPrevious2 - 1;

  // Get the indices of the paths, then get the values
  size_t nIdx11 = nIdxReturn1 * nNbPrevious + nIdxPrevious1;
  size_t nIdx12 = nIdxReturn1 * nNbPrevious + nIdxPrevious2;
  size_t nIdx21 = nIdxReturn2 * nNbPrevious + nIdxPrevious1;
  size_t nIdx22 = nIdxReturn2 * nNbPrevious + nIdxPrevious2;

  double dVal11 = pathDepStruct.LinearInterpolate(nIdx11, dS);
  double dVal12 = pathDepStruct.LinearInterpolate(nIdx12, dS);
  double dVal21 = pathDepStruct.LinearInterpolate(nIdx21, dS);
  double dVal22 = pathDepStruct.LinearInterpolate(nIdx22, dS);
  
  // linear iterpolation in Z and P.  P is done first.  
  double dVal1 = 0.0;
  double dVal2 = 0.0;

  // if off the grid, use a similarity reduction to get the value
  if ( dP > pdPreviousGrid[nNbPrevious - 1] )
  {
    size_t nIdxMax1 = nIdxReturn1 * nNbPrevious + nNbPrevious - 1; 
    dVal1 = m_pdSimilarityValues[nIdxMax1];
    
    size_t nIdxMax2 = nIdxReturn2 * nNbPrevious + nNbPrevious - 1; 
    dVal2 = m_pdSimilarityValues[nIdxMax2];
  }
  else if ( dP < pdPreviousGrid[0] )
  {
    size_t nIdxMin1 = nIdxReturn1 * nNbPrevious; 
    dVal1 = m_pdSimilarityValues[nIdxMin1];

    size_t nIdxMin2 = nIdxReturn2 * nNbPrevious; 
    dVal2 = m_pdSimilarityValues[nIdxMin2];
  }
  else
  {
    double dWeight = (dP - pdPreviousGrid[nIdxPrevious1])
      / (pdPreviousGrid[nIdxPrevious2] - pdPreviousGrid[nIdxPrevious1]);

    dVal1 = (1.0 - dWeight) * dVal11 + dWeight * dVal12;
    dVal2 = (1.0 - dWeight) * dVal21 + dWeight * dVal22;
  }

  double dWeight = (dZ - pdAvgSqrReturnGrid[nIdxReturn1])
    / (pdAvgSqrReturnGrid[nIdxReturn2] - pdAvgSqrReturnGrid[nIdxReturn1]);

  double dVal = (1.0 - dWeight) * dVal1 + dWeight * dVal2;

  return dVal;  

} // LinearInterpolate


double VarianceSwapEvent::QuadraticInterpolate(double dS, double dP, double dZ,
    double dRSquared, const std::vector<double>& pdPreviousGrid,
    const std::vector<double> &pdAvgSqrReturnGrid, 
    PathDepStructure& pathDepStruct) const
{
  int nPath = pathDepStruct.m_path.size();
  size_t nNbReturn = pdAvgSqrReturnGrid.size();
  size_t nNbPrevious = pdPreviousGrid.size();

  // Interpolate
  int nIdxPrevious2 = 0;
  int nIdxPrevious1 = 0;

  if ( !m_bHasSimilarityReduction )
  {
    nIdxPrevious2 = BinSearch(&pdPreviousGrid[0], nNbPrevious, dP);
    nIdxPrevious1 = nIdxPrevious2 - 1;
  }

  // do a binary search to find dZ
  int nZUp   = BinSearch(&pdAvgSqrReturnGrid[0], nNbReturn, dZ);
  int nZLow  = nZUp - 1;
  int nZ2Up  = -1;

  //Check if we are above or below the diagonal
  //if below the diag the flow of info is flowing
  //away from the diagonal to this point
  if ( dZ < dRSquared ) 
  {     
    nZ2Up = nZUp + 1;
         
    if( int(nZ2Up * nNbPrevious + nIdxPrevious1) > nPath - 1  
      || int(nZ2Up * nNbPrevious + nIdxPrevious2) > nPath - 1)    
      nZ2Up = -1;    
    else
    { 
      // sec upstream pt above diag
      if( pdAvgSqrReturnGrid[nZ2Up] > dRSquared)      
        nZ2Up = -1;      
    }// end inner else  

  }// end dZ < dRSquared    
  else // dZ >= dRSquared 
  { 
    size_t nTmp = nZUp;          
    nZUp = nZLow ;   // info flows from diag to this pt
                     // from dZ=dRSquared to here 
   
    nZ2Up = nZLow -1;
    nZLow = nTmp;
    
    if( nZ2Up < 0)              
      nZ2Up = -1;      
    else
    { 
      // sec upstream pt below diag
      if( pdAvgSqrReturnGrid[nZ2Up] < dRSquared) 
        nZ2Up = -1;                     
    }// end inner else    

  }// end dZ >= dRSquared
  
  if ( nZ2Up != -1 ) //quadratic interpolation
  {
      
    int nIdx11 = nZLow * nNbPrevious + nIdxPrevious1;    
    int nIdx12 = nZLow * nNbPrevious + nIdxPrevious2;
    
    int nIdx21 = nZUp * nNbPrevious + nIdxPrevious1;
    int nIdx22 = nZUp * nNbPrevious + nIdxPrevious2;

    int nIdx31 = nZ2Up * nNbPrevious + nIdxPrevious1;
    int nIdx32 = nZ2Up * nNbPrevious + nIdxPrevious2;
    
    //We can use quadratic approximation
    double dZ2Up = pdAvgSqrReturnGrid[ nZ2Up ];
    double dZLow = pdAvgSqrReturnGrid[ nZLow ];
    double dZUp  = pdAvgSqrReturnGrid[ nZUp ];

    double dValUp   = 0.0;
    double dValLow  = 0.0; 
    double dVal2Up  = 0.0;

    if ( m_bHasSimilarityReduction )
    {
      dValLow = m_pdSimilarityValues[nIdx11];
      dValUp  = m_pdSimilarityValues[nIdx21];
      dVal2Up = m_pdSimilarityValues[nIdx31]; 
    }
    else
    {
      //Still Linear in P direction
      double dVal11 = pathDepStruct.LinearInterpolate(nIdx11, dS);
      double dVal12 = pathDepStruct.LinearInterpolate(nIdx12, dS);
      double dVal21 = pathDepStruct.LinearInterpolate(nIdx21, dS);
      double dVal22 = pathDepStruct.LinearInterpolate(nIdx22, dS);
      double dVal31 = pathDepStruct.LinearInterpolate(nIdx31, dS);
      double dVal32 = pathDepStruct.LinearInterpolate(nIdx32, dS);

      double dWeight = (dP - pdPreviousGrid[nIdxPrevious1])
      / (pdPreviousGrid[nIdxPrevious2] - pdPreviousGrid[nIdxPrevious1]);

      dValLow = (1.0 - dWeight) * dVal11 + dWeight * dVal12;
      dValUp  = (1.0 - dWeight) * dVal21 + dWeight * dVal22;
      dVal2Up = (1.0 - dWeight) * dVal31 + dWeight * dVal32;

    }
    
    double dW2Up = (dZ - dZUp) /  ( dZ2Up - dZUp) * (dZ - dZLow)/(dZ2Up - dZLow);
    double dWUp  = (dZ - dZLow) / ( dZUp - dZLow)*(dZ - dZ2Up)/(dZUp - dZ2Up);
    double dWLow = (dZ - dZ2Up) / ( dZLow - dZ2Up)*(dZ - dZUp)/(dZLow - dZUp);

    double dTmp = dValUp * dWUp + dValLow * dWLow + dVal2Up * dW2Up;
   
    if ( dTmp > m_pdSimilarityValues[ m_pdSimilarityValues.size() - 1] )
      dTmp = m_pdSimilarityValues[ m_pdSimilarityValues.size() - 1];

    return dTmp;
    
  }
            
  return LinearInterpolate(dS, dP, dZ, pdPreviousGrid, pdAvgSqrReturnGrid,
                           pathDepStruct);

} // Quadratic Interpolate


double VarianceSwapEvent::QuadZDiagPInterpolate(
          double dP, double dZ, double dRSquared, 
          const std::vector<double>& pdPreviousGrid,
          const std::vector<double> &pdAvgSqrReturnGrid) const
{
  size_t nNbReturn = pdAvgSqrReturnGrid.size();
  size_t nNbPrevious = pdPreviousGrid.size();

  // Find closest paths in Z direction. 2nd upstream point found later.
  size_t nZUp  = BinSearch(&pdAvgSqrReturnGrid[0], nNbReturn, dZ);
  size_t nZLow = nZUp - 1;
  size_t nZ2Up = 0;
  bool bHas2Up = false;
  
  // Find the 2nd upstream point. Check if we are above or below R^2. 
  // The flow is away from this point.
  if ( dZ < dRSquared ) 
  {     
    nZ2Up = nZUp + 1;
    bHas2Up = true;
         
    if( nZ2Up >= nNbReturn )
      bHas2Up = false;
    else
    { 
      // sec upstream pt above diag
      if( pdAvgSqrReturnGrid[nZ2Up] > dRSquared )
        bHas2Up = false;   
    }// end inner else  

  }// end dZ < dRSquared    
  else // dZ >= dRSquared 
  { 
    // info flows from dZ=dRSquared to here 
    size_t nTmp = nZUp;          
    nZUp = nZLow ;
    nZLow = nTmp;
        
    if ( nZUp > 0 && pdAvgSqrReturnGrid[nZUp - 1] >= dRSquared )
    {
      nZ2Up = nZUp - 1;
      bHas2Up = true;
    }
    else
      bHas2Up = false;    

  }// end dZ >= dRSquared

  // Interpolate in P direction first    
  double dValUp   = 0.0;
  double dValLow  = 0.0; 
  double dVal2Up  = 0.0;
    
  if ( m_bHasSimilarityReduction )
  {     
    dValLow = m_pdSimilarityValues[nZLow]; 
    dValUp = m_pdSimilarityValues[nZUp];

    if (bHas2Up)
      dVal2Up = m_pdSimilarityValues[nZ2Up];
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

    if ( bHas2Up )
    {
      size_t nIdx31 = nZ2Up * nNbPrevious + nIdxPrevious1;
      size_t nIdx32 = nZ2Up * nNbPrevious + nIdxPrevious2;

      double dVal31 = m_pdSimilarityValues[nIdx31];
      double dVal32 = m_pdSimilarityValues[nIdx32];

      dVal2Up = (1.0 - dWeight) * dVal31 + dWeight * dVal32;      
    }
  
  }

  // Now interpolate in Z direction
  double dTmp = 0.0;
  if ( bHas2Up )
  {
    // Quadratic in Z direction
    double dZ2Up = pdAvgSqrReturnGrid[ nZ2Up ];
    double dZLow = pdAvgSqrReturnGrid[ nZLow ];
    double dZUp  = pdAvgSqrReturnGrid[ nZUp ];

    double dW2Up = (dZ - dZUp) / ( dZ2Up - dZUp)*(dZ - dZLow)/(dZ2Up - dZLow);
    double dWUp  = (dZ - dZLow) / ( dZUp - dZLow)*(dZ - dZ2Up)/(dZUp - dZ2Up);
    double dWLow = (dZ - dZ2Up) / ( dZLow - dZ2Up)*(dZ - dZUp)/(dZLow - dZUp);

    dTmp = dValUp*dWUp + dValLow*dWLow + dVal2Up*dW2Up;

    // Try to get a monotone interpolation.  Needed, for example, with call 
    // or put payoffs when the three values for interpolation are 0, 0, x.
    // The interpolated value can be negative in this case, which is wrong.
    double dMax = (dValLow < dValUp) ? 
                  ( ( dValUp < dVal2Up ) ? dVal2Up : dValUp ) : 
                  ( (dValLow < dVal2Up) ? dVal2Up : dValLow );
    double dMin = (dValLow < dValUp) ? 
                  ( ( dValLow < dVal2Up ) ? dValLow : dVal2Up ) : 
                  ( (dValUp < dVal2Up) ? dValUp : dVal2Up ); 

    if (dTmp < dMin)
      dTmp = dMin;

    if (dTmp > dMax)
      dTmp = dMax;

  }
  else
  {
    // Only linear in Z direction
    double dZLow = pdAvgSqrReturnGrid[ nZLow ];
    double dZUp  = pdAvgSqrReturnGrid[ nZUp ];

    double dWeight = (dZ - dZLow) / (dZUp - dZLow);

    dTmp = (1.0 - dWeight) * dValLow  + dWeight * dValUp;
  }

  if ( dTmp > m_pdSimilarityValues[ m_pdSimilarityValues.size() - 1] )
    dTmp = m_pdSimilarityValues[ m_pdSimilarityValues.size() - 1];

  return dTmp;

}

} // namespace pricing

} // namespace ito33
