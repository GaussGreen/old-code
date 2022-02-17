/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/minconstconstraint.cpp
// Purpose:     a constant minimum constraint class
// Author:      (z)
// Created:     2003/10/16
// RCS-ID:      $Id: minconstconstraint.cpp,v 1.8 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/minconstconstraint.cpp
    @brief The implementation of the minimum constant constraint class.
*/

#include <cmath>
#include "ito33/pricing/minconstconstraint.h"

using ito33::pricing::MinConstConstraint;

void MinConstConstraint::Apply(double* pdPrices, int* piFlags,
                               size_t nNbValues) const
{
  size_t 
    nIdx;

  // Initialize flags to 0. Avoid to initialize the flags 
  // in if(){}else(){} check.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
    piFlags[nIdx] = 0;

  if (!m_bOn)
    return;

  double dTol;
  
  dTol = m_dTol * fabs(m_dValue);
  if (dTol < m_dTol)
    dTol = m_dTol;

  // If the constraint applies at a given node, update both piFlags and
  // pdPrices.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  {
    if (pdPrices[nIdx] < m_dValue - dTol)
      piFlags[nIdx] = 1;
    
    if (pdPrices[nIdx] < m_dValue + dTol)
      pdPrices[nIdx] = m_dValue;
  }
}

void MinConstConstraint::ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
   size_t nNbValues) const
{ 
  size_t nIdx;

  // Initialize flags to 0. Avoid to initialize the flags in if(){}else(){} check.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
    piFlags[nIdx] = 0;

  if (!m_bOn)
    return;

  // If the constraint applies at a given node, update both piFlags and
  // pdPrices.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  { 
    if (pdPrices[nIdx] <= m_dValue)
    {
      piFlags[nIdx] = 1;
      pdPrices[nIdx] = m_dValue;
    }
  }
  
}


void MinConstConstraint::ApplyPenalty(const double* pdPrices, 
                                      int* piFlags, 
                                      double* pdConstraints, 
                                      size_t nNbValues) const
{ 
  size_t nIdx;

  // Initialize flags to 0. Avoid to initialize the flags in if(){}else(){}
  // check.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
    piFlags[nIdx] = 0;

  if (!m_bOn)
    return;

  // If the constraint applies at a given node, update both piFlags and
  // pdConstraints.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  { 
    if (pdPrices[nIdx] < m_dValue)
    {
      piFlags[nIdx] = 1;
      pdConstraints[nIdx] = m_dValue;
    }
  }
  
}
