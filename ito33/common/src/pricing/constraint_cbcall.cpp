/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/constraint_cbcall.cpp
// Purpose:     implementation of the ConstraintCBCall class
// Author:      Nabil
// Created:     2003/12/12
// RCS-ID:      $Id: constraint_cbcall.cpp,v 1.7 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/constraint_cbcall.cpp
    @brief The implementation of the maximum constraint class.

    Classes to define maximum constraints (such as a call for a cb).
*/

#include "ito33/debug.h"

#include <cmath>
#include "ito33/pricing/constraint_cbcall.h"

using ito33::pricing::ConstraintCBCall;


void ConstraintCBCall::Apply
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
  
  // for a normal max-constraint, when price is closed to constraint value
  // for sake of security, we consider that constraint is not applied.
  for (nIdx = 0; nIdx < m_nIdxStartConversion; nIdx++)
  {
    if (pdPrices[nIdx] > m_pdValues[nIdx] + m_pdTols[nIdx])
      piFlags[nIdx] = m_iFlagValue;
    
    if (pdPrices[nIdx] >= m_pdValues[nIdx] - m_pdTols[nIdx])
      pdPrices[nIdx] = m_pdValues[nIdx];
  }

  // This is different from normal max-constraint after m_nIdxStartConversion
  // when the price is a bit smaller than conversion value, it is still
  // greater than call strike, forced conversion is still available.  
  for (; nIdx < m_nNbValues; nIdx++)
  {
    if (pdPrices[nIdx] > m_pdValues[nIdx] - m_pdTols[nIdx])
    {
      piFlags[nIdx] = m_iFlagValue;
      pdPrices[nIdx] = m_pdValues[nIdx];
    }
  }

}


void ConstraintCBCall::Update(double *pdValues,
                              size_t nNbValues,
                              size_t nIdxStartConversion)
{
  Constraint::Update(pdValues, nNbValues);

  m_nIdxStartConversion = nIdxStartConversion;
  
  if (nIdxStartConversion > nNbValues)
    m_nIdxStartConversion = nNbValues;
}


void ConstraintCBCall::ApplyWithoutTolerance
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
  // pdPrices.  Treat as a normal max constraint.
  for (nIdx = 0; nIdx < m_nNbValues; nIdx++)
  {
    if (pdPrices[nIdx] >= m_pdValues[nIdx])
    {
      piFlags[nIdx] = m_iFlagValue;
      pdPrices[nIdx] = m_pdValues[nIdx];
    }
  }

}


void ConstraintCBCall::ApplyPenalty(const double* pdPrices, 
                                    int* piFlags, 
                                    double* pdConstraints, 
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
  // pdPrices.  Treat as a normal max constraint.
  for (nIdx = 0; nIdx < m_nNbValues; nIdx++)
  {
    if ( pdPrices[nIdx] > m_pdValues[nIdx] )
    {
      piFlags[nIdx] = m_iFlagValue;
      pdConstraints[nIdx] = m_pdValues[nIdx];
    }
  }

}
