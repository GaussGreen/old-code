/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/resetevent.cpp
// Purpose:     reset event
// Author:      Ito33 team
// Created:     05/10/2004
// RCS-ID:      $Id: resetevent.cpp,v 1.10 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/binarysearch.h"


#include "ito33/pricing/resetevent.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"
#include "ito33/pricing/conversionprovisions.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"
#include "ito33/numeric/predicatedouble.h"

namespace ito33
{

namespace pricing
{

    
void ResetEvent::ApplyAtStartTime(PathDepStructure& ) const
{
  // Nothing to do for discrete event
}

void ResetEvent::ApplyAtEndTime(PathDepStructure& ) const
{
  // Nothing to do for discrete event
 
}

void ResetEvent::TurnOffPaths(PathDepStructure& pathDepStruct) const
{
 
  // For non-zero cap and floor rates, the outer active paths need to 
  // obtain information from paths that do not exisit.  Hence, these
  // paths cannot be correctly updated, and can be turned off. At the
  // first reset date (going forward in time), only the current 
  // conversion price path should remain active
  size_t nPath = pathDepStruct.m_path.size();
  size_t nIdPath;

  // The last reset date (backward direction) is special
  if ( m_bIsLastReset)
  {
    for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++ )
    {
      double dConvPrice = pathDepStruct.m_pppdGrids[0][0][nIdPath];

      if ( numeric::IsEqual(dConvPrice, m_dCurrentConvPrice) )
      {
        pathDepStruct.m_pbIsActive[nIdPath] = true;
        pathDepStruct.SetPathToSave(nIdPath);
      }
      else
      {
        pathDepStruct.m_pbIsActive[nIdPath] = false;
      } //end if

    } // end loop over each path

    return;
  } // if last reset date

  // Find the outer active paths for this reset period
  size_t nFloorIndex = 0;
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++ )
  {
    if (pathDepStruct.m_pbIsActive[nIdPath] = true)
    {
      nFloorIndex = nIdPath;
      break;
    }
  }

  size_t nCapIndex = 0;
  for ( nIdPath = nPath-1 ; nIdPath > 0 ; nIdPath-- )
  {
    if (pathDepStruct.m_pbIsActive[nIdPath] = true)
    {
      nCapIndex = nIdPath;
      break;
    }
  }

  // Find the new cap and floor conversion prices 
  double dCap = pathDepStruct.m_pppdGrids[0][0][nCapIndex] / m_dCapRate;

  double dFloor = pathDepStruct.m_pppdGrids[0][0][nFloorIndex];
  if ( m_bIsResetFlooredBy_Prevailing )
    dFloor /= m_dFloorRate;

  // Turn off paths for which we cannot get valid data
  for ( nIdPath = nFloorIndex ; nIdPath <= nCapIndex ; nIdPath++ )
  {
    double dConvPrice = pathDepStruct.m_pppdGrids[0][0][nIdPath];

    if ( !numeric::IsEqualOrGreater(dConvPrice, dFloor) 
      || !numeric::IsEqualOrLess(dConvPrice, dCap) )
    {
      pathDepStruct.m_pbIsActive[nIdPath] = false;
    }

  } // end loop over each path

} // ResetEvent::TurnOffPaths(PathDepStructure& pathDepStruct)


void ResetEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{

  //ApplySecondMethod(pathDepStruct);
  //return;

  // turn off the paths which need not be updated
  TurnOffPaths(pathDepStruct); 

  size_t nPath   = pathDepStruct.m_path.size();
  size_t nIdPath;

  // Only one path. This can only happen when the cap rate and floor rate 
  // are the same and equal to 1.0
  // TODO: maybe we could have an exception
  // message for this. However, this is handy for debugging

  if ( nPath == 1 && m_bIsLastReset )
  {
    pathDepStruct.SetPathToSave(0);
    return;
  }

  if ( nPath == 1 )
    return;

  //create a local copy of the
  //conversion ratio grid
  //easier for future used
  //in particular for binary search
  std::vector<double> pGridY(nPath);

  //contains the solution matrix 
  std::vector< std::vector<double> > V(nPath);

  //initialization of solution matrix and y grid
  for ( nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    //local copy of y grid
     pGridY[nIdPath] = pathDepStruct.m_pppdGrids[0][0][nIdPath];

     //std::cout << nIdPath << " " << pGridY[nIdPath] << std::endl;

    size_t nGridSize  = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    size_t nIdS;
    
    V[nIdPath].resize(nGridSize);  
    
    for ( nIdS = 0 ; nIdS < nGridSize; nIdS++ )
    {
      V[nIdPath][nIdS] = 0.0;
    } //end loop over each grid
  
  } //end loop over paths


  //---------------------------------------------------------------------------
  //
  // Do the communication
  //
  //---------------------------------------------------------------------------
  //loop over each path
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {
   
    // No point in updating a path that is not active
    if (pathDepStruct.m_pbIsActive[nIdPath] == false)
      continue;

    size_t nGridSize     = pathDepStruct.m_path[nIdPath].meshes->GetNbS();    
    const double* pdGrid = pathDepStruct.m_path[nIdPath].meshes->GetS();

    // get the conversion ratio and price of this path
    double dConvPrice = pathDepStruct.m_pppdGrids[0][0][nIdPath];

    // Compute the max/min (cap/floor) conversion price values
    double dCapPrice   = m_dCapRate*dConvPrice;
    double dFloorPrice = 0.0;

    if ( m_bIsResetFlooredBy_Prevailing  )
      dFloorPrice = m_dFloorRate * dConvPrice;
    else
      dFloorPrice = m_dFloorRate * m_dInitialConvPrice;


    // loop over each grid
    size_t nIdS;
    for ( nIdS = 0 ; nIdS < nGridSize; nIdS++ )
    {
      // get the spot in the current grid
      double dS = pdGrid[nIdS];      

      // calculate a new (capped/floored) conversion price
      double dNewConvPrice = std::max( dFloorPrice, 
                             std::min( dCapPrice, dS * m_dMultiplier ) );

      //Step 1 find the path that corresponds to that
      //conversion price
      //Recall:
      //if x value lies between A[l] <= x <= A[u] returns index = u
      //if x < A[0], index = 1
      //if x > A[n-1], index = n-1
      size_t nIdConv = ito33::BinSearch(&pGridY[0], nPath, dNewConvPrice);
      
      //conversion rate k_i is \in [nIdConv-1,nIdConv]

      //Step 2 find the S value corresponding to pdGrid[nIdS]
      //as grids are going to be different

      //value of the path below
      const size_t nGridSizeBelow  
        = pathDepStruct.m_path[nIdConv-1].meshes->GetNbS();

      const double* pdGridBelow   
            = pathDepStruct.m_path[nIdConv-1].meshes->GetS();
      
      const double* pdValuesBelow     
            = pathDepStruct.m_path[nIdConv-1].instdata->m_pdPrices.Get();
      
      double dSBelow = dS * (pGridY[nIdConv-1] / dNewConvPrice);
      //double dSBelow = dS;
      size_t nIdSPathBelow = BinSearch(pdGridBelow, nGridSizeBelow, dSBelow);

      //value of the path above
      const size_t nGridSizeAbove  
        = pathDepStruct.m_path[nIdConv].meshes->GetNbS();

      const double* pdGridAbove
            = pathDepStruct.m_path[nIdConv].meshes->GetS();
      
      const double* pdValuesAbove     
            = pathDepStruct.m_path[nIdConv].instdata->m_pdPrices.Get();
      
      double dSAbove = dS * (pGridY[nIdConv] / dNewConvPrice);
      //double dSAbove = dS;
      size_t nIdSPathAbove = BinSearch(pdGridAbove, nGridSizeAbove, dSAbove);

     // 2D linear interpolation

     // horizontal interpolation between (S1, k1) and (S2, k1)
     double dS1 = pdGridAbove[ nIdSPathAbove - 1 ];
     double dS2 = pdGridAbove[ nIdSPathAbove ];
     double dV1 = pdValuesAbove[ nIdSPathAbove - 1 ];
     double dV2 = pdValuesAbove[ nIdSPathAbove ];

     double dVAbove = dV1 + (dV2 - dV1) / (dS2 - dS1) * (dSAbove - dS1);

     // horizontal interpolation between (S3, k2) and (S4, k2)
     double dS3 = pdGridBelow[ nIdSPathBelow - 1 ];
     double dS4 = pdGridBelow[ nIdSPathBelow] ;
     double dV3 = pdValuesBelow[ nIdSPathBelow - 1 ];
     double dV4 = pdValuesBelow[ nIdSPathBelow ];

     double dVBelow = dV3 + (dV4 - dV3) / (dS4 - dS3) * (dSBelow - dS3);

     V[nIdPath][nIdS]
        = dVBelow + (dVAbove - dVBelow) / (dSAbove - dSBelow) * (dS - dSBelow);

     // Old interpolation scheme. Need to set dSBelow and dSAbove to dS
     //V[nIdPath][nIdS] = dVBelow + 
     //  (dVAbove - dVBelow) / (pGridY[nIdConv] - pGridY[nIdConv-1]) * 
     //  (dNewConvPrice - pGridY[nIdConv-1]);

    } //end loop over nodes

  } //end loop over path


  //---------------------------------------------------------------------------
  //
  // Copy the solution appropriately
  //
  //---------------------------------------------------------------------------
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {

    size_t nGridSize     = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    size_t nIdS;
    double* pdValues     = pathDepStruct.m_path[nIdPath].instdata->m_pdPrices.Get();

    for ( nIdS = 0 ; nIdS < nGridSize; nIdS++)
    {
      // On the resets date values get reset no matter what
      pdValues[nIdS] = V[nIdPath][nIdS];
    } //end loop grid

  } //end loop path

}//end ResetEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const


void ResetEvent::ApplySecondMethod(PathDepStructure& pathDepStruct) const
{

  // Turn off paths which need not be updated
  TurnOffPaths(pathDepStruct);

  size_t nPath   = pathDepStruct.m_path.size();
  size_t nIdPath;

  // Only one path.  This can only happen when the cap rate and floor rate 
  // are the same and equal to 1.0. Path to save is updated in TurnOffPaths.
  // TODO: maybe we could have an exception message for this. However, this 
  // is handy for debugging
  if ( nPath == 1 )
    return;

  // Create a local copy of the conversion ratio grid for future use (in
  // particular, binary search)
  std::vector<double> pGridY(nPath);

  // Contains the solution matrix 
  std::vector< std::vector<double> > V(nPath);

  // Initialization of solution matrix and y grid
  size_t nGridSizeMax = 0;
  for ( nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    //local copy of y grid
    pGridY[nIdPath] = pathDepStruct.m_pppdGrids[0][0][nIdPath];

    size_t nGridSize  = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    if (nGridSize > nGridSizeMax)
      nGridSizeMax = nGridSize;
    
    V[nIdPath].resize(nGridSize);  

  } //end loop over paths

  // Will need some work arrays below
  std::vector<double> pdWorkArray1(nGridSizeMax);
  std::vector<double> pdWorkArray2(nGridSizeMax);

  // Ignoring caps and floors, the new conversion price is S * multiplier.
  // Create a grid with the appropriate values from the paths.  Each path
  // then needs to interpolate from this "master" grid, and apply the
  // appropriate caps and floors.  The grid size is simply the number
  // of paths. If S * multiplier is out of range, it must be in the
  // cap or floor region.
  std::vector<double> pdMasterS(nPath);
  std::vector<double> pdMasterValues(nPath);

  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {
    // The conversion price of this path
    double dConvPrice = pGridY[nIdPath];

    // The S value that needs to get data from this path
    double dS = dConvPrice / m_dMultiplier;

    // Interpolate within this path for the value at S
    const double* pdGrid = pathDepStruct.m_path[nIdPath].meshes->GetS();
    const double* pdValues = pathDepStruct.m_path[nIdPath].instdata->m_pdPrices.Get();
    const size_t nGridSize = pathDepStruct.m_path[nIdPath].meshes->GetNbS();

    size_t nIdS = BinSearch(pdGrid, nGridSize, dS);

    double dSlope = (pdValues[nIdS] - pdValues[nIdS-1]) / 
                    (pdGrid[nIdS] - pdGrid[nIdS-1]);
    double dValue = pdValues[nIdS-1] + dSlope * (dS - pdGrid[nIdS-1]);

    // Store the values
    pdMasterS[nIdPath] = dS;
    pdMasterValues[nIdPath] = dValue;
    
  } // loop over paths to construct master grid


  // loop over each path. For values of S that lead to conversion prices
  // over the cap, the new values come from the same path.  Ditto for
  // values below the floor.  Intermediate prices come from the master
  // grid constructed above.
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {

    // No point in updating a path that is not active
    if (pathDepStruct.m_pbIsActive[nIdPath] == false)
      continue;

    // Get the conversion ratio and price of this path
    double dConvPrice = pGridY[nIdPath];

    // Compute the max/min (cap/floor) conversion price values
    double dCapPrice   = m_dCapRate * dConvPrice;
    double dFloorPrice = 0.0;

    if ( m_bIsResetFlooredBy_Prevailing  )
      dFloorPrice = m_dFloorRate * dConvPrice;
    else
      dFloorPrice = m_dFloorRate * m_dInitialConvPrice;

    // Get the grid of the current path
    size_t nGridSize     = pathDepStruct.m_path[nIdPath].meshes->GetNbS();    
    const double* pdGrid = pathDepStruct.m_path[nIdPath].meshes->GetS();

    // Interpolate for the floor values
    double dSGridFloor = dFloorPrice / m_dMultiplier;
    size_t nIdGridFloor = BinSearch(pdGrid, nGridSize, dSGridFloor);

    size_t nIdFloorPath = BinSearch(&pGridY[0], nPath, dFloorPrice);
    const size_t nGridSizeFloorBelow 
      = pathDepStruct.m_path[nIdFloorPath-1].meshes->GetNbS();
    const double* pdGridFloorBelow 
      = pathDepStruct.m_path[nIdFloorPath-1].meshes->GetS();
    const double* pdValuesFloorBelow 
      = pathDepStruct.m_path[nIdFloorPath-1].instdata->m_pdPrices.Get();

    const size_t nGridSizeFloorAbove 
      = pathDepStruct.m_path[nIdFloorPath].meshes->GetNbS();
    const double* pdGridFloorAbove 
      = pathDepStruct.m_path[nIdFloorPath].meshes->GetS();
    const double* pdValuesFloorAbove 
      = pathDepStruct.m_path[nIdFloorPath].instdata->m_pdPrices.Get();

    // Do the horizontal interpolation in work arrays
    numeric::Interpolate(pdGridFloorBelow, pdValuesFloorBelow, 
      nGridSizeFloorBelow, pdGrid, &pdWorkArray1[0], nIdGridFloor, 
      numeric::ExtrapolationMode_Linear, numeric::ExtrapolationMode_Linear);

    numeric::Interpolate(pdGridFloorAbove, pdValuesFloorAbove, 
      nGridSizeFloorAbove, pdGrid, &pdWorkArray2[0], nIdGridFloor, 
      numeric::ExtrapolationMode_Linear, numeric::ExtrapolationMode_Linear);

    // Now do vertical interpolation between the paths
    size_t nIdS;
    for ( nIdS = 0; nIdS < nIdGridFloor; nIdS++)
    {
      double dVBelow = pdWorkArray1[nIdS];
      double dVAbove = pdWorkArray2[nIdS];

      V[nIdPath][nIdS] = dVBelow + 
        (dVAbove - dVBelow) / (pGridY[nIdFloorPath] - pGridY[nIdFloorPath-1]) 
        * (dFloorPrice - pGridY[nIdFloorPath-1]);

    } // floor prices


    // Now do the cap prices
    double dSCap = dCapPrice / m_dMultiplier;
    size_t nIdGridCap = BinSearch(pdGrid, nGridSize, dSCap);

    size_t nIdCapPath = BinSearch(&pGridY[0], nPath, dCapPrice);
    const size_t nGridSizeCapBelow 
      = pathDepStruct.m_path[nIdCapPath-1].meshes->GetNbS();
    const double* pdGridCapBelow 
      = pathDepStruct.m_path[nIdCapPath-1].meshes->GetS();
    const double* pdValuesCapBelow 
      = pathDepStruct.m_path[nIdCapPath-1].instdata->m_pdPrices.Get();

    const size_t nGridSizeCapAbove 
      = pathDepStruct.m_path[nIdCapPath].meshes->GetNbS();
    const double* pdGridCapAbove 
      = pathDepStruct.m_path[nIdCapPath].meshes->GetS();
    const double* pdValuesCapAbove 
      = pathDepStruct.m_path[nIdCapPath].instdata->m_pdPrices.Get();

    // Do the horizontal interpolation in work arrays
    numeric::Interpolate(pdGridCapBelow, pdValuesCapBelow, nGridSizeCapBelow,
      &pdGrid[nIdGridCap], &pdWorkArray1[nIdGridCap], nGridSize - nIdGridCap, 
      numeric::ExtrapolationMode_Linear, numeric::ExtrapolationMode_Linear);

    numeric::Interpolate(pdGridCapAbove, pdValuesCapAbove, nGridSizeCapAbove,
      &pdGrid[nIdGridCap], &pdWorkArray2[nIdGridCap], nGridSize - nIdGridCap, 
      numeric::ExtrapolationMode_Linear, numeric::ExtrapolationMode_Linear);

    for ( nIdS = nIdGridCap; nIdS < nGridSize; nIdS++)
    {

      double dVBelow = pdWorkArray1[nIdS];
      double dVAbove = pdWorkArray2[nIdS];

      V[nIdPath][nIdS] = dVBelow + 
        (dVAbove - dVBelow) / (pGridY[nIdCapPath] - pGridY[nIdCapPath-1]) * 
        (dCapPrice - pGridY[nIdCapPath-1]);

    } // cap prices


    // now do the intermediate prices
    numeric::Interpolate(&pdMasterS[0], &pdMasterValues[0], nPath,
      &pdGrid[nIdGridFloor], &(V[nIdPath][nIdGridFloor]),
      nIdGridCap-nIdGridFloor,
      numeric::ExtrapolationMode_Linear, numeric::ExtrapolationMode_Linear);

  } //end loop over path


  //---------------------------------------------------------------------------
  //
  // Copy the solution appropriately
  //
  //---------------------------------------------------------------------------
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {

    if (pathDepStruct.m_pbIsActive[nIdPath] == false)
      continue;

    size_t nGridSize = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    double* pdValues = pathDepStruct.m_path[nIdPath].instdata->m_pdPrices.Get();

    size_t nIdS;
    for ( nIdS = 0 ; nIdS < nGridSize; nIdS++)
    {
      // On the resets date values get reset no matter what
      pdValues[nIdS] = V[nIdPath][nIdS];
    } //end loop grid

  } //end loop path


} //end ResetEvent::ApplySecondMethod


} // namespace pricing

} // namespace ito33
