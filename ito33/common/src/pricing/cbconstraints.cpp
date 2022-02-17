/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/cbconstraints.cpp
// Purpose:     a cb constraints class
// Author:      Nabil
// Created:     2003/12/12
// RCS-ID:      $Id: cbconstraints.cpp,v 1.15 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/vector.h"

#include "ito33/pricing/cbconstraints.h"
#include "ito33/pricing/minconstraint.h"
#include "ito33/pricing/constraint_cbcall.h"

using namespace ito33::pricing;

void 
CBConstraints::Apply(double *pdPrices, int *piFlags, size_t nNbSpots) const
{
  size_t nIdx;

  // Initialize the flags here since it will not be done by the individual 
  // constraints.
  for (nIdx = 0; nIdx < nNbSpots; nIdx++)
    piFlags[nIdx] = State_none;

  // Call
  m_call.Apply(pdPrices, piFlags, nNbSpots);
  
  // Conversion
  m_conv.Apply(pdPrices, piFlags, nNbSpots);
  
  // Put
  m_put.Apply(pdPrices, piFlags, nNbSpots);

}


void CBConstraints::ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
                     size_t nNbSpots) const
{
  size_t nIdx;

  // Initialize the flags here since it will not be done by the individual 
  // constraints.
  for(nIdx = 0; nIdx < nNbSpots; nIdx++)
    piFlags[nIdx] = State_none;

  // Call
  m_call.ApplyWithoutTolerance(pdPrices, piFlags, nNbSpots);
  
  // Conversion
  m_conv.ApplyWithoutTolerance(pdPrices, piFlags, nNbSpots);

  // Put
  m_put.ApplyWithoutTolerance(pdPrices, piFlags, nNbSpots);

}


void CBConstraints::ApplyPenalty(const double* pdPrices, int* piFlags, 
                                 double* pdConstraints, size_t nNbSpots) const
{
  size_t nIdx;

  // Initialize the flags here since it will not be done by the individual 
  // constraints.
  for(nIdx = 0; nIdx < nNbSpots; nIdx++)
    piFlags[nIdx] = State_none;

  // Call
  m_call.ApplyPenalty(pdPrices, piFlags, pdConstraints, nNbSpots);
  
  // Conversion
  m_conv.ApplyPenalty(pdPrices, piFlags, pdConstraints, nNbSpots);
  
  // Put
  m_put.ApplyPenalty(pdPrices, piFlags, pdConstraints, nNbSpots);

}


void CBConstraints::Get(double *pdPrices, int *piFlags, size_t nNbSpots) const
{
  // Each individual constraint will only update pdPrices and/or piFlags
  // where piFlags is set to the individual constraint flag value. For
  // example, call only updates where piFlags[nIdx] == State_call.

  // Call
  m_call.Get(pdPrices, piFlags, nNbSpots);
  
  // Conversion
  m_conv.Get(pdPrices, piFlags, nNbSpots);

  // Put
  m_put.Get(pdPrices, piFlags, nNbSpots);
}


void CBConstraints::UpdatePut(double *pdValues, size_t nNbValues)
{
  m_put.Update(pdValues, nNbValues);
}


void CBConstraints::UpdateConv(double *pdValues, size_t nNbValues)
{
  m_conv.Update(pdValues, nNbValues);
}


void CBConstraints::UpdateCall(double *pdValues,
                               size_t nNbValues,
                               size_t nIdxStartConversion)
{
  m_call.Update(pdValues, nNbValues, nIdxStartConversion);
}
