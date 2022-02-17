/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/sharedependentconversionevent.cpp
// Purpose:     Share dependent conversion event
// Author:      Ito33 team
// Created:     2005/01/14
// RCS-ID:      $Id: sharedependentconversionevent.cpp,v 1.5 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/binarysearch.h"

#include "ito33/pricing/sharedependentconversionevent.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/predicatedouble.h"

namespace ito33
{

namespace pricing
{


void ShareDependentConversionEvent::TurnOffPaths(PathDepStructure& pathDepStruct) const
{
 
  size_t nPath = pathDepStruct.m_path.size();
  size_t nIdPath;

  ASSERT_MSG( numeric::IsEqual(m_dBaseRatio, pathDepStruct.m_pppdGrids[0][0][0]),
              "Error when turning off path, could not find the base ratio path.");

  //deactivate all path except the path 0
   //who has by default the right conversion stuff
  for ( nIdPath = 1 ; nIdPath < nPath ; nIdPath++ )
  {
    pathDepStruct.m_pbIsActive[nIdPath] = false;
  } // end loop over each path

  //ensure that path 0 is activated
  pathDepStruct.m_pbIsActive[0] = true;
  pathDepStruct.SetPathToSave(0);
 
} // ShareDependendConversionEvent::TurnOffPaths(PathDepStructure& pathDepStruct)


void ShareDependentConversionEvent::ApplyAtTime(
     PathDepStructure& pathDepStruct) const
{

  ASSERT_MSG( numeric::IsEqual(m_dBaseRatio, pathDepStruct.m_pppdGrids[0][0][0]),
              "Error could not find the base ratio path.");

  size_t nPath   = pathDepStruct.m_path.size();
  size_t nIdPath;

  //create a local copy of the conversion ratio grid
  //easier for future used in particular for binary search
  std::vector<double> pGridY(nPath);

  //initialization of solution matrix and y grid
  for ( nIdPath = 0; nIdPath < pGridY.size(); nIdPath++)  
  { 
    //local copy of y grid
     pGridY[nIdPath] = pathDepStruct.m_pppdGrids[0][0][nIdPath];
  } //end loop over paths

  //find the strike value in the first grid which corresponds
  //to the base conversion ratio
  const double* pdGridBase = pathDepStruct.m_path[0].meshes->GetS();
  double* pdValuesBase     = pathDepStruct.m_path[0].instdata->m_pdPrices.Get();
  size_t nGridSize         = pathDepStruct.m_path[0].meshes->GetNbS(); 

  std::vector<double> pTmpVal(nGridSize);

  //helper 
  //loop over the path 0 : 
  //     if ( dS < dStrike ) get info from 0
  //     if ( c_r < capratio ) get info from the different paths
  //     if ( c_r >= capratio) get info from the upmost paths

  //since at some point the cap ratio is going to be reached
  //which corresponds to the value of the up most path
  //we precompute the interpolation of the up most path
  //onto path0 
  const double *pdGridCapRatio = pathDepStruct.m_path[nPath-1].meshes->GetS();
  
  const double *pdValuesCapRatio = 
    pathDepStruct.m_path[nPath-1].instdata->m_pdPrices.Get();
  
  size_t nGridSizeCapRatio  = pathDepStruct.m_path[nPath-1].meshes->GetNbS();

  std::vector<double> pdValCapRatioInterp(nGridSize);

  numeric::Interpolate(pdGridCapRatio, pdValuesCapRatio, nGridSizeCapRatio,
    pdGridBase, &pdValCapRatioInterp[0], nGridSize,
    numeric::ExtrapolationMode_Linear, numeric::ExtrapolationMode_Linear);

  //loop over the base grid
  size_t nIdS;

  double dStrike = m_pShareDependentConversion->GetStrike();

  for (  nIdS = 0 ; nIdS < nGridSize; nIdS++)
  {

    double dS = pdGridBase[nIdS];

    double dCurrentRatio = m_dBaseRatio; 
    
    if ( dS > dStrike )
      dCurrentRatio += ( dS - dStrike)/dS 
       * m_dIncrementalShareFactor;

    dCurrentRatio = std::min(dCurrentRatio, m_dCapRatio);

    if ( numeric::IsEqual(dCurrentRatio, m_dBaseRatio) )
     pTmpVal[nIdS]  = pdValuesBase[nIdS];
    else if ( numeric::IsEqual(dCurrentRatio, m_dCapRatio) ) // c_r  = c_max   
     pTmpVal[nIdS] = pdValCapRatioInterp[nIdS];  
    else   
    {
      //Step1 find the conversion ratio path
      //that corresponds to the current conversion ratio
      size_t nIdConv = ito33::BinSearch(&pGridY[0], nPath, dCurrentRatio);

      //value of the path below    
      const size_t nGridSizeBelow  
          = pathDepStruct.m_path[nIdConv-1].meshes->GetNbS();

      const double* pdGridBelow   
          = pathDepStruct.m_path[nIdConv-1].meshes->GetS();
      
      const double* pdValuesBelow     
          = pathDepStruct.m_path[nIdConv-1].instdata->m_pdPrices.Get();
   
      size_t nIdSPathBelow = ito33::BinSearch(pdGridBelow, nGridSizeBelow, dS);

      //value of the path above
      const size_t nGridSizeAbove  
          = pathDepStruct.m_path[nIdConv].meshes->GetNbS();

      const double* pdGridAbove
          = pathDepStruct.m_path[nIdConv].meshes->GetS();
      
      const double* pdValuesAbove     
          = pathDepStruct.m_path[nIdConv].instdata->m_pdPrices.Get();
     
      size_t nIdSPathAbove = ito33::BinSearch(pdGridAbove, nGridSizeAbove, dS);
    
      //2D linear interpolation
      double dS1 = pdGridAbove[ nIdSPathAbove - 1 ];
      double dS2 = pdGridAbove[ nIdSPathAbove ];
      double dV1 = pdValuesAbove[ nIdSPathAbove - 1 ];
      double dV2 = pdValuesAbove[ nIdSPathAbove ];
  
      // horizontal interpolation between (S3, k2) and (S4, k2)
      double dS3 = pdGridBelow[ nIdSPathBelow - 1 ];
      double dS4 = pdGridBelow[ nIdSPathBelow] ;
      double dV3 = pdValuesBelow[ nIdSPathBelow - 1 ];
      double dV4 = pdValuesBelow[ nIdSPathBelow ];
              
      pTmpVal[nIdS] = numeric::LinearInterpolate2D(dS, dCurrentRatio, dS1, dS2,
                      dS3, dS4, pGridY[nIdConv], pGridY[nIdConv-1], dV1, dV2, 
                      dV3, dV4);
    
    }//end if
    
  } //end loop over grid

  //copy solution
  for ( nIdS = 0 ; nIdS < nGridSize; nIdS++)
    pdValuesBase[nIdS] = pTmpVal[nIdS];
  

  //---------------------------------------------------------------------------
  //
  // turn off all the paths except the one that
  // corresponds to the base ratio.
  // The conversion type is changed as the path are turned off.
  //
  //---------------------------------------------------------------------------
  TurnOffPaths(pathDepStruct); 


}//end ShareDependendConversionEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const


} // namespace pricing

} // namespace ito33
