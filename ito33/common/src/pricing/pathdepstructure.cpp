/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/pathdepstructure.cpp
// Purpose:     path dependent strcuture class implementation
// Author:      Yann and David
// Created:     2004/07/07
// RCS-ID:      $Id: pathdepstructure.cpp,v 1.9 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/binarysearch.h"

#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/backwardmeshmanager_withspacemesh.h"
#include "ito33/pricing/instdata.h"

#include "ito33/numeric/interpolation.h"

namespace ito33
{

namespace pricing
{

void PathDepStructure::SetPathToSave(size_t nIdx)
{
  // Make sure the path is within range
  ASSERT_MSG(nIdx < m_path.size(), "Invalid path to save");

  m_nPathToSave = nIdx;
}

double PathDepStructure::LinearInterpolate(size_t nPath, double dS)
{

  size_t nGridSize = m_path[nPath].meshes->GetNbS();
  const double* pdGrid = m_path[nPath].meshes->GetS();
  const double* pdVal = GetPriceData(nPath);

  double dInterpolatedValue = 0.0;

  numeric::LinearInterpolate(pdGrid, pdVal, nGridSize, &dS,
    &dInterpolatedValue, 1, numeric::ExtrapolationMode_Linear, 
    numeric::ExtrapolationMode_Linear); 

  return dInterpolatedValue;
}

double PathDepStructure::QuadraticInterpolate(size_t nPath, double dS)
{
  size_t nGridSize = m_path[nPath].meshes->GetNbS();
  const double* pdGrid = m_path[nPath].meshes->GetS();
  const double* pdVal = GetPriceData(nPath);

  double dInterpolatedValue = 0.0;
  
  numeric::QuadraticInterpolate(pdGrid, pdVal, nGridSize, &dS, 
    &dInterpolatedValue, 1, numeric::ExtrapolationMode_Linear, 
    numeric::ExtrapolationMode_Linear); 
                          
  return dInterpolatedValue;

}

double PathDepStructure::SimilarityInterpolation(size_t nIdPath, 
      double dY, double dS, double dYStar, size_t nHomogeneousDegree,
      bool bIsLinearInterpolation)
{
  double dSInterpolated = dYStar * ( dS / (dY + 1.e-100) );

  double dVal = 0.0;
  
  if ( bIsLinearInterpolation )
    dVal = LinearInterpolate(nIdPath, dSInterpolated); 
  else
    dVal = QuadraticInterpolate(nIdPath, dSInterpolated); 
  
  return pow(dY / dYStar, int(nHomogeneousDegree) ) * dVal;
}

double PathDepStructure::SimilarityInterpolation(size_t nIdPathLow, 
       size_t nIdPathHigh, double dY, double dS, 
       const std::vector<double>& pdGridY, size_t nHomogeneousDegree,
       bool bIsLinearInterpolation)
{

  ASSERT_MSG( pdGridY.size() == 2,  "Grid error for similarity reduction.");
  
  double dYStar = pdGridY[1];

  size_t nGridSize  = m_path[nIdPathHigh].meshes->GetNbS();
  const double* pdS = m_path[nIdPathHigh].meshes->GetS();

  double dSMax = pdS[ nGridSize -1 ];

  double dSInterpolated = dYStar * ( dS / (dY + 1.e-100) );

  if( dSInterpolated > dSMax)
  { //boundary condition interpolation 
    //by using the other path
    //linear interpolation in x
    double dWHigh = ( dY - pdGridY[0]) / (pdGridY[1] - pdGridY[0]);
    double dWLow  = 1. - dWHigh;

    double dValHigh = LinearInterpolate(nIdPathHigh, dS); 
    double dValLow  = LinearInterpolate(nIdPathLow, dS);
      
    return( dValLow * dWLow + dValHigh * dWHigh);
  }
 
  double dVal = 0.0;
  
  if ( bIsLinearInterpolation )
    dVal = LinearInterpolate(nIdPathHigh, dSInterpolated); 
  else
    dVal = QuadraticInterpolate(nIdPathHigh, dSInterpolated); 
  
  return pow(dY / dYStar, int(nHomogeneousDegree) ) * dVal;
}

} // namespace pricing

} // namespace ito33
