/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/domain.cpp
// Purpose:     implementation of the interface domain class
// Created:     2004/05/04
// RCS-ID:      $Id: domain.cpp,v 1.9 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 - 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/finance/domain.cpp
   @brief implementation of the interface domain class
 */

#include "ito33/arraycheckers.h"

#include "ito33/finance/error.h"
#include "ito33/finance/domain.h"

extern const ito33::finance::Error ITO33_DOMAIN_SPOTS_NOT_SET;

using ito33::finance::Domain;

void Domain::SetUnderlyingSharePrices(const Spots& spots)
{
  CheckIncreasingOrder(spots);

  m_spots = spots;
}

const Domain::Spots& Domain::GetUnderlyingSharePrices() const
{
  CHECK_COND(m_spots.size() > 0, ITO33_DOMAIN_SPOTS_NOT_SET);

  return m_spots;
}
