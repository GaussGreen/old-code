/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/averagingevent.cpp
// Purpose:     Averaging Event
// Author:      Ito33 team Canada
// Created:     April 13, 2005
// RCS-ID:      $Id: averagingevent.cpp,v 1.9 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/binarysearch.h"

#include "ito33/pricing/averagingevent.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"

#include "ito33/numeric/predicatedouble.h"

namespace ito33
{

namespace pricing
{


void AveragingEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const
{
 
  size_t nPath = pathDepStruct.m_path.size();
  size_t nIdPath;

  //create a local copy of the averaging grid
  std::vector<double> pdGridY(nPath);

  //contains the solution matrix 
  std::vector< std::vector<double> > V(nPath);

  //initialization of solution matrix and y grid
  for ( nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    //local copy of averaging grid
    pdGridY[nIdPath] = pathDepStruct.m_pppdGrids[0][0][nIdPath];

    size_t nGridSize = pathDepStruct.m_path[nIdPath].meshes->GetNbS();

    V[nIdPath].resize(nGridSize);  
      
  } //end loop over paths


  //loop over each path
  for ( nIdPath = 0; nIdPath < nPath; nIdPath++)  
  { 
    double dA = pdGridY[nIdPath];
    
    size_t nGridSize  = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    size_t nIdS;
    const double* pdS = pathDepStruct.m_path[nIdPath].meshes->GetS();

    for ( nIdS = 0; nIdS < nGridSize; nIdS++)
    {   
   
      double dS = pdS[nIdS];

      //Compute the new average
      double dANew = -1.0;
      
      dANew = GetNewAverage(dA, dS);

      ASSERT_MSG( dANew >= 0.0, "Error computing the average.");

       if ( m_bHasSimilarityReduction )
         V[nIdPath][nIdS] =  pathDepStruct.SimilarityInterpolation(0, dANew,
                              dS, pdGridY[0], 1, false );
       else
          V[nIdPath][nIdS] = Interpolate(dANew, dS, pathDepStruct, pdGridY);

    } // end loop grid
  } // end loop path

  //---------------------------------------------------------------------------
  //
  // Copy the solution appropriately
  //
  //---------------------------------------------------------------------------
  for ( nIdPath = 0 ; nIdPath < nPath ; nIdPath++)
  {

    size_t nGridSize = pathDepStruct.m_path[nIdPath].meshes->GetNbS();
    size_t nIdS;
    double* pdValues = pathDepStruct.GetPriceData(nIdPath);

    for ( nIdS = 0 ; nIdS < nGridSize; nIdS++)
    {
      pdValues[nIdS] = V[nIdPath][nIdS];
    } //end loop grid

  } //end loop path

   //last observation turn off all paths
  if ( m_bRecurse == false &&
      (m_nObservation == 1 || m_bIsLastEvent == true) )
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
 
}//end AveragingEvent::ApplyAtTime(PathDepStructure& pathDepStruct) const

double AveragingEvent::GetNewAverage(double dA, double dS) const
{ 
  double dANew = dA + ( dS - dA )/ double( m_nObservation ); 
       
  return dANew;
} // AveragingEvent::GetNewAverage(double dA, double dS) const

double AveragingEvent::Interpolate(double dY, double dS, 
    PathDepStructure& pathDepStruct, std::vector<double> &pdGridY) const
{
  //do a binary search to find dY
  int nPath = pathDepStruct.m_path.size();
  int nUp   = BinSearch(&pdGridY[0], pdGridY.size(), dY);
  int nLow  = nUp - 1;
  int n2Up     = -1;

  //Check if we are above or below the diagonal
  //if below the diag the flow of info is flowing
  //away from the diagonal to this point
  if( dY < dS ) 
  { //  dY < dS
    // info from diag to this pt (below y=x)
    // from y=x to here
                
    n2Up = nUp + 1;
         
    if( n2Up > nPath -1 )    
      n2Up = -1;    
    else
    { // inner else       
      if( pdGridY[n2Up] > dS) // sec upstream pt above diag     
        n2Up = -1;      
    }// end inner else     
  }// end dY < dS    
  else // dY >= dS 
  { 
    size_t nTmp = nUp;          
    nUp = nLow ;   // info flows from diag to this pt
                   // from y=x to here 
   
    n2Up = nLow -1;
    nLow = nTmp;
    
    if( n2Up < 0)              
      n2Up = -1;      
    else
    { 
      if( pdGridY[n2Up] < dS) // sec upstream pt below diag
        n2Up = -1;                     
    }// end inner else    

  }// end y >= x
  
  if( n2Up != -1 ) //quadratic interpolation
  {
    //We can use quadratic approximation
    double dY2Up = pdGridY[ n2Up ];
    double dYLow = pdGridY[ nLow ];
    double dYUp  = pdGridY[ nUp ];

    double dValUp   = pathDepStruct.LinearInterpolate(nUp, dS);
    double dValLow  = pathDepStruct.LinearInterpolate(nLow,dS); 
    double dVal2Up  = pathDepStruct.LinearInterpolate(n2Up,dS);

    double dW2Up = (dY - dYUp)/(dY2Up - dYUp)*(dY - dYLow)/(dY2Up - dYLow);
    double dWUp  = (dY - dYLow)/( dYUp - dYLow)*(dY - dY2Up)/(dYUp - dY2Up);
    double dWLow = (dY - dY2Up)/( dYLow - dY2Up)*(dY - dYUp)/(dYLow - dYUp);

    double dTmp = dValUp*dWUp + dValLow*dWLow + dVal2Up*dW2Up;
   
    if (dTmp > 0.0)
      return dTmp;
  }
            
  // Use linear interpolation
  double dWHigh = ( dY - pdGridY[nLow])/(pdGridY[nUp] - pdGridY[nLow]);
  double dWLow  = 1. - dWHigh;

  //linear interpolation in x
  double dValHigh = pathDepStruct.LinearInterpolate(nUp,dS); 
  double dValLow  = pathDepStruct.LinearInterpolate(nLow,dS);
      
  return( dValLow * dWLow + dValHigh * dWHigh);


} // AveragingEvent::Interpolate

} // namespace pricing

} // namespace ito33
