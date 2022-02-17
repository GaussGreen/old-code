/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/mesh/gridwithcv.cpp
// Purpose:     Grid with change of variable class
// Author:      Nabil
// Created:     2003/11/19
// RCS-ID:      $Id: gridwithcv.cpp,v 1.4 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file numeric/mesh/gridwithcv.cpp
    @brief Grid with change of variable class

    Implementation of some functions of the GridWithCV class.
    
 */

#include "ito33/numeric/mesh/gridwithcv.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{

size_t GridWithCV::GetIndex(double dS)
{
  size_t
    nMin = 0,
    nMax = m_nNbS - 1,
    nEca = nMax,
    nMid;
  
  int
    iRes = -1;
  
  double
    dTolerance = dS * 1.e-10;
  
  if(m_nNbS <= 0)
    return iRes;

  if(dS < m_pdS[0] - dTolerance || dS > m_pdS[nMax] + dTolerance)
    return iRes;

  if(dS <= m_pdS[0] + dTolerance)
    return 0;

  if(dS >= m_pdS[nMax] - dTolerance)
    return nMax;

  // now, dSpot is always between p[nMin] et p[nMax]
  while(nEca > 1)
    {
    nEca >>= 1;
    nMid = nMax - nEca;
    if(dS < m_pdS[nMid] - dTolerance)
      {
      nMax = nMid;
      nEca = nMax - nMin;
      }
    else if(dS > m_pdS[nMid] + dTolerance)
      nMin = nMid;
    else
      {
      iRes = nMid;
      break;
      }
    }

  return iRes;
}


size_t GridWithCV::GetIndexLog(double dLogS)
{
  size_t
    nMin = 0,
    nMax = m_nNbS - 1,
    nEca = nMax,
    nMid;

  int
    iRes = -1;

  const double DTOLERANCE = 1.e-10;
  
  if(m_nNbS <= 0)
    return iRes;

  if(dLogS < m_pdLogS[0] - DTOLERANCE || dLogS > m_pdLogS[nMax] + DTOLERANCE)
    return iRes;

  if(dLogS <= m_pdLogS[0] + DTOLERANCE)
    return 0;

  if(dLogS >= m_pdLogS[nMax] - DTOLERANCE)
    return nMax;

  // now, dLogSpot is always between p[nMin] et p[nMax]
  while(nEca > 1)
    {
    nEca >>= 1;
    nMid = nMax - nEca;
    if(dLogS < m_pdLogS[nMid] - DTOLERANCE)
      {
      nMax = nMid;
      nEca = nMax - nMin;
      }
    else if(dLogS > m_pdLogS[nMid] + DTOLERANCE)
      nMin = nMid;
    else
      {
      iRes = nMid;
      break;
      }
    }

  return iRes;
}

} // namespace mesh

} // namespace numeric

} // namespace ito33
