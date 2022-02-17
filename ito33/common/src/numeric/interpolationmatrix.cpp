/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/interpolationmatrix.cpp
// Purpose:     Implementation of InterpolationMatrix
// Created:     2005/06/04
// RCS-ID:      $Id: interpolationmatrix.cpp,v 1.7 2005/09/15 12:38:57 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @todo only linear interpolation is implemented
 */

#include "ito33/autoptr.h"
#include "ito33/array.h"
#include "ito33/debug.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::InterpolationMatrix);
}

namespace ito33
{

namespace numeric
{


double 
WeightedPoints::InterpolateWith
(const double* pdValues) const
{
  double dTmp = 0.;
  for (const_iterator i = begin(); i != end(); ++i)
    dTmp += pdValues[i->m_nIdx] * i->m_dWeight;

  return dTmp;
}

void 
WeightedPoints::TransposeInterpolate
(double dValue, double* pdValues) const
{
  for (const_iterator i = begin(); i != end(); ++i)
    pdValues[i->m_nIdx] += dValue * i->m_dWeight;
}

InterpolationMatrix::InterpolationMatrix
(const double* pdX, size_t nNbPoints, const double* pdS, size_t nNbS,
 size_t nNbSubSystem, bool bIsBackward,
 ExtrapolationMode emLeft, ExtrapolationMode emRight,
 InterpolationMethod interpMethod)
{
  m_nNbPoints = nNbPoints;
  m_nNbS = nNbS,
  m_nNbSubSystem = nNbSubSystem;
  m_bIsBackward = bIsBackward;

  Array<int> piIntervals(nNbPoints);

  SearchIntervals(pdS, nNbS, pdX, piIntervals.Get(), nNbPoints);

  const int iNbIntervals = (int)nNbS - 1;

  m_points.resize(nNbPoints);

  for (size_t nIdxP = 0; nIdxP < nNbPoints; nIdxP++)
  {
    int iInterval = piIntervals[nIdxP];

    int iI = iInterval;
    
    if ( iInterval == -1)
    {  
      if ( emLeft == ExtrapolationMode_Linear )
      {
        iI = 0;
        double dWeight = (pdX[nIdxP] - pdS[iI]) / (pdS[iI + 1] - pdS[iI]);

        m_points[nIdxP].push_back( WeightedPoint(iI, 1 - dWeight) );
        m_points[nIdxP].push_back( WeightedPoint(iI + 1, dWeight) );
      }
      else if ( emLeft == ExtrapolationMode_Constant )
        m_points[nIdxP].push_back(WeightedPoint(0, 1));
    }
    else if ( iInterval == iNbIntervals )
    {
      if ( emRight == ExtrapolationMode_Linear )
      {
        iI = iNbIntervals - 1;
      
        double dWeight = (pdX[nIdxP] - pdS[iI]) / (pdS[iI + 1] - pdS[iI]);

        m_points[nIdxP].push_back( WeightedPoint(iI, 1 - dWeight) );
        m_points[nIdxP].push_back( WeightedPoint(iI + 1, dWeight) );
      }
      else if ( emRight == ExtrapolationMode_Constant )
        m_points[nIdxP].push_back(WeightedPoint(nNbS - 1, 1));
    }
    else
    {
      if ( interpMethod == InterpolationMethod_Linear ||  iNbIntervals == 1)
      {
        double dWeight = (pdX[nIdxP] - pdS[iI]) / (pdS[iI + 1] - pdS[iI]);

        m_points[nIdxP].push_back( WeightedPoint(iI, 1 - dWeight) );
        m_points[nIdxP].push_back( WeightedPoint(iI + 1, dWeight) );
      }
      else if ( interpMethod == InterpolationMethod_Quadratic )
      {
        if ( iInterval == iNbIntervals - 1 )
          iI--; 

        double dX1 = pdS[iI];
        double dX2 = pdS[iI + 1];
        double dX3 = pdS[iI + 2];

        double dX = pdX[nIdxP];
        double dW1 = (dX - dX2) * (dX - dX3) / ( (dX1 - dX2) * (dX1 - dX3) );
        double dW2 = (dX - dX1) * (dX - dX3) / ( (dX2 - dX1) * (dX2 - dX3) );
        double dW3 = (dX - dX1) * (dX - dX2) / ( (dX3 - dX1) * (dX3 - dX2) );
      
        m_points[nIdxP].push_back( WeightedPoint(iI, dW1) );
        m_points[nIdxP].push_back( WeightedPoint(iI + 1, dW2) );
        m_points[nIdxP].push_back( WeightedPoint(iI + 2, dW3) );
      }
    }
  }
}

void InterpolationMatrix::MultiplyBy(double dScalar)
{
  for (size_t nIdxP = 0; nIdxP < m_nNbPoints; nIdxP++)
  {
    WeightedPoints& elements = m_points[nIdxP];
    WeightedPoints::iterator i;
    for ( i = elements.begin(); i != elements.end(); ++i )
      i->m_dWeight *= dScalar;
  }
}

void InterpolationMatrix::Add(size_t nIdxP1, size_t nIdxP2, double dWeight)
{
  ASSERT_MSG( nIdxP1 < m_points.size(), "Element requested does not exist." );
  m_points[nIdxP1].push_back( WeightedPoint(nIdxP2, dWeight) );
}

void 
InterpolationMatrix::ProductMatrixVector
(const double* pdX, double* pdF, bool bAddOn) const
{
  if ( !bAddOn )
  {
    size_t nNbValues;

    if ( m_bIsBackward )
      nNbValues = m_nNbPoints * m_nNbSubSystem;
    else
      nNbValues = m_nNbPoints;

    for (size_t nIdx = 0; nIdx < nNbValues; nIdx++)
      pdF[nIdx] = 0.;
  }

  if ( m_bIsBackward )
  {
    for (size_t nIdxP = 0; nIdxP < m_nNbPoints; nIdxP++)
    {
      const WeightedPoints& elements = m_points[nIdxP];
      for (size_t nIdxSS = 0; nIdxSS < m_nNbSubSystem; nIdxSS++)
        pdF[nIdxP + nIdxSS * m_nNbPoints]
            += elements.InterpolateWith(pdX + nIdxSS * m_nNbS);
    }
  }
  else
  {
    for (size_t nIdxP = 0; nIdxP < m_nNbPoints; nIdxP++)
    {
      const WeightedPoints& elements = m_points[nIdxP];
      for (size_t nIdxSS = 0; nIdxSS < m_nNbSubSystem; nIdxSS++)
        pdF[nIdxP] += elements.InterpolateWith(pdX + nIdxSS * m_nNbS);
    }
  }
}

void 
InterpolationMatrix::ProductTransposeMatrixVector
(const double* pdX, double* pdF, bool bAddOn) const
{
  if ( !bAddOn )
    for (size_t nIdx = 0; nIdx < m_nNbS * m_nNbSubSystem; nIdx++)
      pdF[nIdx] = 0.;

  if ( m_bIsBackward )
  {
    for (size_t nIdxP = 0; nIdxP < m_nNbPoints; nIdxP++)
    {
      const WeightedPoints& elements = m_points[nIdxP];
      for (size_t nIdxSS = 0; nIdxSS < m_nNbSubSystem; nIdxSS++)        
        elements.TransposeInterpolate
                 (pdX[nIdxP + nIdxSS * m_nNbPoints], pdF + nIdxSS * m_nNbS);
    }
  }
  else
  {
    for (size_t nIdxP = 0; nIdxP < m_nNbPoints; nIdxP++)
    {
      const WeightedPoints& elements = m_points[nIdxP];
      for (size_t nIdxSS = 0; nIdxSS < m_nNbSubSystem; nIdxSS++)
        elements.TransposeInterpolate(pdX[nIdxP], pdF + nIdxSS * m_nNbS);
    }
  }
}


} // namespace numeric

} // namespace ito33
