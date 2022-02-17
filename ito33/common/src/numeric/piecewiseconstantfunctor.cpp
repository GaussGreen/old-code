/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/piecewiseconstantfunctor.cpp
// Purpose:     Class for piecewise constant function
// Created:     2004/06/04
// RCS-ID:      $Id: piecewiseconstantfunctor.cpp,v 1.5 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/arraycheckers.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

using ito33::numeric::PiecewiseConstantFunctor;

PiecewiseConstantFunctor::PiecewiseConstantFunctor
                          (const double* pdX, const double* pdY, size_t nNbX)
{
  // For now, allow zero sized vectors, and assume the caller has a reason
  // for constructing and empty object

  m_pdX.resize(nNbX);
  m_pdY.resize(nNbX);
  
  if (nNbX > 0)
  {

    ito33::CheckIncreasingOrder(pdX, nNbX);
 
    for (size_t nIdx = 0; nIdx < nNbX; nIdx++)
    {
      m_pdX[nIdx] = pdX[nIdx];
      m_pdY[nIdx] = pdY[nIdx];
    }
  } // if nNbX is > 0

}

double PiecewiseConstantFunctor::operator()(double dX) const
{
  // If no data is available, then we could throw an excpetion.
  // However, for now, just return zero, and assume that the caller 
  // has a reason for constructing this object without any data
  if ( m_pdX.size() == 0)
    return 0.0;

  if (dX < m_pdX[0] || m_pdX.size() == 1)
    return m_pdY[0];

  size_t nIdx = 1;

  while ( nIdx < m_pdX.size() - 1 && m_pdX[nIdx] <= dX )
    nIdx++;

  return m_pdY[nIdx];
}
