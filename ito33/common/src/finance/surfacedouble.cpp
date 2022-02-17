/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/surfacedouble.cpp
// Purpose:     implementation of the interface surface class for double 
// Author:      Wang
// Created:     2004/05/04
// RCS-ID:      $Id: surfacedouble.cpp,v 1.9 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 -   Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/finance/surfacedouble.cpp
   @brief implementation of the interface surface class for double. Ex:prices
 */

#include "ito33/sharedptr.h"

#include "ito33/numeric/domain.h"
#include "ito33/numeric/surfacedouble.h"

#include "ito33/finance/domain.h"
#include "ito33/finance/surfacedouble.h"

namespace ito33
{
  
namespace finance
{


void SurfaceDouble::GetValuesAt(size_t nDate, Doubles &values) const
{
  m_pSurface->GetValuesAt
              (
                m_pSurface->GetDomain()->GetTimeIndexFromDateIndex(nDate), 
                GetDomain()->GetUnderlyingSharePrices(), 
                values
              );
}

double SurfaceDouble::GetValueAt(size_t nDate, double dSpot) const
{
  Domain::Spots pdS(1);

  pdS[0] = dSpot;

  Doubles valuesTmp;

  m_pSurface->GetValuesAt
              (
                m_pSurface->GetDomain()->GetTimeIndexFromDateIndex(nDate),
                pdS, 
                valuesTmp
              );

  return valuesTmp[0];
}

shared_ptr<Domain> SurfaceDouble::GetDomain() const
{
  return m_pSurface->GetDomain();
}
  
shared_ptr<numeric::SurfaceDouble> SurfaceDouble::GetImpl() const
{
  return m_pSurface;
}


} // namespace finance

} // namespace ito33
