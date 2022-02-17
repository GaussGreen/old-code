/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/constraint.cpp
// Purpose:     a base class for single constraint 
// Author:      Wang
// Created:     2003/10/14
// RCS-ID:      $Id: constraint.cpp,v 1.15 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/constraint.cpp
    @brief The implementation of the single constraint class.
*/

#include <cmath>

#include "ito33/debug.h"
#include "ito33/array.h"

#include "ito33/pricing/constraint.h"

using ito33::pricing::Constraint;

Constraint::Constraint(const double *pdValues, size_t nNbValues, 
                       double dTolerance, size_t nMaxSize,
                       int iFlagValue, bool bClearBeforeApply)
  : m_iFlagValue(iFlagValue), 
    m_bClearBeforeApply(bClearBeforeApply)    
{
  ASSERT_MSG(nNbValues, 
             "The size of the constraints array should be non-null!");
  
  ASSERT_MSG(dTolerance >= 0,
             "The tolerance should be positive or zero!");

  m_dTol = dTolerance;

  m_nMaxSize = (nMaxSize < nNbValues) ? nNbValues : nMaxSize;
  m_nNbValues = nNbValues;

  m_pdValues = Array<double>(m_nMaxSize);
  m_pdTols = Array<double>(m_nMaxSize);
  
  // Store the constraint values
  for (size_t nIdx = 0; nIdx < m_nNbValues; nIdx++)
  {
    m_pdValues[nIdx] = pdValues[nIdx];
    
    m_pdTols[nIdx] = fabs(m_pdValues[nIdx]) * m_dTol;    
    if (m_pdTols[nIdx] < m_dTol)
      m_pdTols[nIdx] = m_dTol;
  }

  m_bOn = true;
}


void Constraint::Update(double *pdValues, size_t nNbValues)
{
  ASSERT_MSG(nNbValues, 
             "The size of the constraints array must be non-null");

  if (m_nMaxSize < nNbValues)
  {
    m_nMaxSize = nNbValues;

    m_pdValues = Array<double>(m_nMaxSize);
    m_pdTols   = Array<double>(m_nMaxSize);
  }
  
  // Store the constraint values
  m_nNbValues = nNbValues;
  for (size_t nIdx = 0; nIdx < m_nNbValues; nIdx++)
  {
    m_pdValues[nIdx] = pdValues[nIdx];

    m_pdTols[nIdx] = fabs(m_pdValues[nIdx]) * m_dTol;
  }

  m_bOn = true;
}


void Constraint::TurnOn()
{
  ASSERT_MSG(m_nNbValues > 0,
             "Constraint:: TurnOn() called for a void constraint object.");
  
  m_bOn = true; 
}


void Constraint::SetTolerance(double dTolerance)
{
  ASSERT_MSG(dTolerance >= 0,
             "The tolerance must be positive or zero!");

  m_dTol = dTolerance;

  // have to do this even when m_bOn is off, because of the TurnOn function
  for (size_t nIdx = 0; nIdx < m_nNbValues; nIdx++)
  {
    m_pdTols[nIdx] = fabs(m_pdValues[nIdx]) * m_dTol; 
    if (m_pdTols[nIdx] < m_dTol)
      m_pdTols[nIdx] = m_dTol;
  }
}


void Constraint::Get(double *pdPrices, int *piFlags, size_t nNbValues) const
{
  ASSERT_MSG(!m_bOn || m_nNbValues == nNbValues, 
             "sizes not equal when applying constraint");

  size_t nIdx;

  if (m_bOn)
  {
    // Set the prices to the constraint values for the active constraints
    for (nIdx = 0; nIdx < m_nNbValues; nIdx++)
      if (piFlags[nIdx] == m_iFlagValue)
        pdPrices[nIdx] = m_pdValues[nIdx];
  }
  else
  {
    // Since the constraint is off, turn off the flags that are
    // set to this constraint
    for (nIdx = 0; nIdx < nNbValues; nIdx++)
      if (piFlags[nIdx] == m_iFlagValue)
        piFlags[nIdx] = 0;
  }
}
