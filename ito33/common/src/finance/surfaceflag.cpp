/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/surfaceflag.cpp
// Purpose:     implementation of the interface surface class for flag 
// Author:      Wang
// Created:     2004/05/04
// RCS-ID:      $Id: surfaceflag.cpp,v 1.7 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 -   Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/finance/surfaceflag.cpp
   @brief implementation of the interface surface class for flags.
 */

#include "ito33/sharedptr.h"

#include "ito33/numeric/domain.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/finance/domain.h"
#include "ito33/finance/surfaceflag.h"

namespace ito33
{
  
namespace finance
{


void SurfaceFlag::GetValuesAt(size_t nDate, Flags &values) const
{
  m_pSurface->GetValuesAt
              (
                m_pSurface->GetDomain()->GetTimeIndexFromDateIndex(nDate), 
                GetDomain()->GetUnderlyingSharePrices(), 
                values
              );
}

SurfaceFlag::Flag SurfaceFlag::GetValueAt(size_t nDate, double dSpot) const
{
  Domain::Spots pdS(1);

  pdS[0] = dSpot;

  Flags valuesTmp;

  m_pSurface->GetValuesAt
              (
                m_pSurface->GetDomain()->GetTimeIndexFromDateIndex(nDate),
                pdS, 
                valuesTmp
              );

  return valuesTmp[0];
}

shared_ptr<Domain> SurfaceFlag::GetDomain() const
{
  return m_pSurface->GetDomain();
}


} // namespace finance

} // namespace ito33
