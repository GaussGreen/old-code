/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/minconstraint.cpp
// Purpose:     implementation of the MinConstraint class
// Author:      Wang
// Created:     2003/10/15
// RCS-ID:      $Id: minconstraint.cpp,v 1.14 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/minconstraint.cpp
    @brief The implementation of the minimum constraint class.

    Classes to define minimum constraints (such as an American option).
*/

#include <cmath>

#include "ito33/debug.h"
#include "ito33/autoptr.h"

#include "ito33/pricing/minconstraint.h"

using ito33::pricing::MinConstraint;


namespace ito33
{

  // implement the autoptrdeleter of the MinConstraint class
  ITO33_IMPLEMENT_AUTOPTR(pricing::MinConstraint);

} // namespace ito33


void MinConstraint::Apply
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
    if (pdPrices[nIdx] < m_pdValues[nIdx] - m_pdTols[nIdx])
      piFlags[nIdx] = m_iFlagValue;
    
    if (pdPrices[nIdx] <= m_pdValues[nIdx] + m_pdTols[nIdx])
      pdPrices[nIdx] = m_pdValues[nIdx];
  }
}

 
void MinConstraint::ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
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
    if ( pdPrices[nIdx] <= m_pdValues[nIdx] )
    {
      piFlags[nIdx] = m_iFlagValue;
      pdPrices[nIdx] = m_pdValues[nIdx];
    }  
  }

}


void MinConstraint::ApplyPenalty(const double* pdPrices, int* piFlags, 
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
    if ( pdPrices[nIdx] < m_pdValues[nIdx] )
    {
      piFlags[nIdx] = m_iFlagValue;
      pdConstraints[nIdx] = m_pdValues[nIdx];
    }  
  }

}
