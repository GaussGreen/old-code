/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/ndayscall.cpp
// Purpose:     instrument can be called at time 
//              t only when the stock has been traded 
//              above trigger during N days prior to t.
//              Copy the solution from above if S is above trigger
//              Copy the solution from below is S is below the trigger.
//              Solution is copied below because we are considering
//              consecutive days and not cumulative.
// Author:      Yann and David
// Created:     08/09/2004
// RCS-ID:      $Id: ndayscallevent.cpp,v 1.15 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/binarysearch.h"
#include "ito33/constants.h"

#include "ito33/pricing/ndayscallevent.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"
#include "ito33/pricing/conversionprovisions.h"
#include "ito33/pricing/cblikeparams.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"

namespace ito33
{

namespace pricing
{
    
void NDaysCallEvent::ApplyAtStartTime(PathDepStructure& pathDepStruct) const
{
  //Start time event, so solving backward in time
  //this event is applied last
  //basically we deactivating all uncessary block.
  //Block corresponds to a type of contingent conversion
  //          0..n-1   pde  corresponds to cannot convert at all
  //          n..2n-1  pdes corresponds to conversion next quarter
  //         2n...3n-1 pdes corresponds to conversion anytime forever
  // So a block corresponds to the number of n contingent call days
  // deactivate all pdes except the base path in each block 0,n,2n.

  size_t nIdPath  = 0;
  size_t nPath = pathDepStruct.m_path.size();

  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {
    pathDepStruct.m_pbIsActive[nIdPath] = false; 
  } //end for loop

  //loop over the different block and activate
  //the appropriate pde
  size_t nNbBlocks = nPath / m_nCallPaths;
  for (size_t nIdxBlock = 0; nIdxBlock < nNbBlocks; nIdxBlock++)
    pathDepStruct.m_pbIsActive[nIdxBlock*m_nCallPaths] = true; 

} // end ApplyAtStartTime


void NDaysCallEvent::ApplyAtEndTime(PathDepStructure& pathDepStruct) const
{
  //this events is called first when going backward in time
  // We need to ensure that at the base pde is active.

  // Check that the base pde in each conversion block is active
  ASSERT_MSG( pathDepStruct.m_pbIsActive[0] == true,
   "Error the base path is not activated");

  // Activate all intermediate pdes
  size_t nIdPath = 0;
  size_t nPath   = pathDepStruct.m_path.size();

  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {
    pathDepStruct.m_pbIsActive[nIdPath] = true;
  }

  // copy the lowest pde in each block to all other pdes in the block
  // that is for
  // pdes 0..n-1    copy pde[0] into 1...n-1
  // pdes n..2n-1   copy pde[n] into n+1..2n-1
  // pdes 2n...3n-1 copy pde[2n] into 2n+1..3n-1
  size_t nNbBlocks = nPath / m_nCallPaths;
  for (size_t nIdxBlock = 0; nIdxBlock < nNbBlocks; nIdxBlock++)
  {

    // Get the base offset into the current block
    size_t nOffset = nIdxBlock * m_nCallPaths;

    // Get the base path data that is pde[0], pde[n], pde[2n]
    size_t nGridSize0     = pathDepStruct.m_path[nOffset].meshes->GetNbS();
    const double* pdGrid0 = pathDepStruct.m_path[nOffset].meshes->GetS();
    double* pdValues0    
      = pathDepStruct.m_path[nOffset].instdata->m_pdPrices.Get();

    //loop over each path inside the block
    for ( nIdPath = 1 ; nIdPath < m_nCallPaths ; nIdPath++)
    {
      size_t nGridSizenIdx 
        = pathDepStruct.m_path[ nOffset+nIdPath ].meshes->GetNbS();
      
      const double* pdGridnIdx  
        = pathDepStruct.m_path[ nOffset+nIdPath ].meshes->GetS();
      
      double* pdValuesnIdx = 
        pathDepStruct.m_path[ nOffset+nIdPath ].instdata->m_pdPrices.Get();

      Interpolate(pdGrid0, pdValues0, nGridSize0, 
        pdGridnIdx, pdValuesnIdx, nGridSizenIdx,
        numeric::ExtrapolationMode_Linear,
        numeric::ExtrapolationMode_Linear);

    }//end loop over different paths

  } // end loop over blocks
  
  // If this call event spans the valuation time, then set the path to save
  // to be the trigger history
  double dValuationTime = pathDepStruct.m_path[0].params->GetValuationTime();
  if ( numeric::IsEqualOrBefore(m_dStartTime, dValuationTime) )
    pathDepStruct.SetPathToSave(m_nPathToSave);
}


void NDaysCallEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{

  ASSERT_MSG(pathDepStruct.m_path.size() % m_nCallPaths == 0,
    "The total number paths is not a multiple of the call paths.");

  // Compute the call trigger  
  // Same for each path
  double dTime    = pathDepStruct.m_path[0].meshes->GetTime();
    
  finance::TriggerAsPercentageOf 
    trigger = m_pConversionProvisions->GetTriggerAsPercentageOf();

  double dTrigger = m_pConversionProvisions->GetConversionPrice(dTime, trigger);
  dTrigger *= m_dTriggerRate;

  // In each block, copy from path i+1 to path i for values above the trigger.
  // If below the trigger, copy from the base base of the block (to
  // effectively reset the counter)
  size_t nNbBlocks = pathDepStruct.m_path.size() / m_nCallPaths;
  for (size_t nIdxBlock = 0; nIdxBlock < nNbBlocks; nIdxBlock++)
  {

    // Get the base offset into the current block
    size_t nOffset = nIdxBlock * m_nCallPaths;

    size_t nGridSize0     = pathDepStruct.m_path[nOffset].meshes->GetNbS();
    const double* pdGrid0 = pathDepStruct.m_path[nOffset].meshes->GetS();
    double* pdValues0     
                = pathDepStruct.m_path[nOffset].instdata->m_pdPrices.Get();

    size_t nIdPath  = 0;
    //need to go up to m_nCallPaths+1 since we
    //also want the last pde of each block to be resetted
    for ( nIdPath = 1; nIdPath < m_nCallPaths+1; nIdPath++)
    {   
      // copy pdeIplus into pdeI
      size_t nGridSizeI 
        = pathDepStruct.m_path[ nOffset+nIdPath-1 ].meshes->GetNbS();
      
      const double* pdGridI
        = pathDepStruct.m_path[ nOffset+nIdPath-1 ].meshes->GetS();
      
      double* pdValuesI = 
        pathDepStruct.m_path[ nOffset+nIdPath-1].instdata->m_pdPrices.Get();
     

      size_t nPathAboveIndex = nOffset+nIdPath ;
     
      //check to avoid going above the last pde.
      if ( nIdPath >= m_nCallPaths )
      {
         nPathAboveIndex = nOffset + m_nCallPaths-1;
      }
     
      size_t nGridSizeIPlus 
        = pathDepStruct.m_path[ nPathAboveIndex ].meshes->GetNbS();
      
      const double* pdGridIPlus 
        = pathDepStruct.m_path[ nPathAboveIndex].meshes->GetS();
      
      double* pdValuesIPlus = 
          pathDepStruct.m_path[ nPathAboveIndex ].instdata->m_pdPrices.Get();

      

      Array<double> pdInterpValuesAbove(nGridSizeI);

      Interpolate(pdGridIPlus, pdValuesIPlus, nGridSizeIPlus,
                  pdGridI, pdInterpValuesAbove.Get(), nGridSizeI,
                  numeric::ExtrapolationMode_Linear,
                  numeric::ExtrapolationMode_Linear);
 
      Array<double> pdInterpValuesBelow(nGridSizeI);

      Interpolate(pdGrid0, pdValues0, nGridSize0,
                  pdGridI, pdInterpValuesBelow.Get(), nGridSizeI,
                  numeric::ExtrapolationMode_Linear,
                  numeric::ExtrapolationMode_Linear);
      
      size_t nIdxTrigger = BinSearch(pdGridI,nGridSizeI,dTrigger);
      size_t nIdx = 0;
 
      if (    nIdxTrigger > 0 
           && !numeric::LessThanTrigger(pdGridI[nIdxTrigger - 1], dTrigger) ) 
        nIdxTrigger--;

      for ( nIdx = 0 ; nIdx < nGridSizeI; nIdx++)
      {

        if ( nIdx >= nIdxTrigger ) 
        {
          // Copy above the trigger
          pdValuesI[nIdx] = pdInterpValuesAbove[nIdx];
        }
        else
        {
          //copy below trigger
          //effectively resetting the counter
          pdValuesI[nIdx] = pdInterpValuesBelow[nIdx];   
        }

      } //end loop over grid

    } //end loop over different path


  } // loop over blocks

} //end ApplyAtTime


} // namespace pricing

} // namespace ito33
