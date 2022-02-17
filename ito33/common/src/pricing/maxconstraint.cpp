/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/maxconstraint.cpp
// Purpose:     implementation of the MaxConstraint class
// Author:      Nabil
// Created:     2003/12/12
// RCS-ID:      $Id: maxconstraint.cpp,v 1.11 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/maxconstraint.cpp
    @brief The implementation of the maximum constraint class.

    Classes to define maximum constraints (such as a call for a cb).
*/

#include "ito33/debug.h"

#include <cmath>
#include "ito33/pricing/maxconstraint.h"

using ito33::pricing::MaxConstraint;


void MaxConstraint::Apply
     (double* pdPrices, int* piFlags, size_t nNbValues) const
{
  ASSERT_MSG(!m_bOn || m_nNbValues == nNbValues, 
             "sizes not equal when applying constraint");

  size_t nIdx;

  // Initialize flags to 0 (default to no constraint) even if the 
  // constraint is not on, unless specifically set otherwise
  if ( m_bClearBeforeApply )
    for (nIdx = 0; nIdx < nNbValues; nIdx++)
      piFlags[nIdx] = 0;

  if (!m_bOn)
    return;

  // If the constraint applies at a given node, update both piFlags and
  // pdPrices.
  for (nIdx = 0; nIdx < m_nNbValues; nIdx++)
  {
    if (pdPrices[nIdx] > m_pdValues[nIdx] + m_pdTols[nIdx])
      piFlags[nIdx] = m_iFlagValue;
    
    if (pdPrices[nIdx] >= m_pdValues[nIdx] - m_pdTols[nIdx])
      pdPrices[nIdx] = m_pdValues[nIdx];
  }
}


void MaxConstraint::ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
    size_t nNbValues) const
{  
  ASSERT_MSG(!m_bOn || m_nNbValues == nNbValues, 
             "sizes not equal when applying constraint");

  size_t nIdx;

  // Initialize flags to 0 (default to no constraint) even if the 
  // constraint is not on, unless specifically set otherwise
  if ( m_bClearBeforeApply )
    for (nIdx = 0; nIdx < nNbValues; nIdx++)
      piFlags[nIdx] = 0;

  if (!m_bOn)
    return;

  // If the constraint applies at a given node, update both piFlags and
  // pdPrices.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  {
    if (pdPrices[nIdx] >= m_pdValues[nIdx])
    {
      piFlags[nIdx] = m_iFlagValue;
      pdPrices[nIdx] = m_pdValues[nIdx];
    }
  }

}


void MaxConstraint::ApplyPenalty(const double* pdPrices, int* piFlags, 
                                 double* pdConstraints, size_t nNbValues) const
{  
  ASSERT_MSG(!m_bOn || m_nNbValues == nNbValues, 
             "sizes not equal when applying constraint");

  size_t nIdx;

  // Initialize flags to 0 (default to no constraint) even if the 
  // constraint is not on, unless specifically set otherwise
  if ( m_bClearBeforeApply )
    for (nIdx = 0; nIdx < nNbValues; nIdx++)
      piFlags[nIdx] = 0;

  if (!m_bOn)
    return;

  // If the constraint applies at a given node, update both piFlags and
  // pdConstraints.
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  {
    if (pdPrices[nIdx] > m_pdValues[nIdx])
    {
      piFlags[nIdx] = m_iFlagValue;
      pdConstraints[nIdx] = m_pdValues[nIdx];
    }
  }

}
