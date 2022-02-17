/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/surfacezero.cpp
// Purpose:     implementation of the double value surface class
// Created:     2004/05/04
// RCS-ID:      $Id: surfacezero.cpp,v 1.9 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 - 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/numeric/surfacezero.cpp
   @brief implementation of the zero value surface class
 */

#include "ito33/sharedptr.h"

#include "ito33/numeric/surfacezero.h"

namespace ito33
{

namespace numeric
{

typedef shared_ptr<SurfaceDouble> SurfacePtr;

void SurfaceZero::GetValuesAt
      (size_t /* nIdxT */, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize(pdS.size(), 0);
}

void SurfaceZero::GetFirstValuesAt
      (size_t /* nIdxT */, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize(pdS.size(), 0);
}

void SurfaceZero::GetLastValuesAt
      (size_t /* nIdxT */, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize(pdS.size(), 0);
}

void SurfaceZero::GetFirstValuesAt
      (double /* dTime*/, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize(pdS.size(), 0);
}

void SurfaceZero::GetLastValuesAt
      (double /* dTime */, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize(pdS.size(), 0);
}

void SurfaceZero::GetDeltaAndGamma
                  (SurfacePtr& pDeltaSurface, SurfacePtr& pGammaSurface) const
{
  pDeltaSurface = SurfacePtr( new SurfaceZero(m_pDomain) );
  pGammaSurface = SurfacePtr( new SurfaceZero(m_pDomain) );
}

SurfacePtr 
SurfaceZero::ComputeFiniteDifference(const SurfaceDouble& shiftedSurface,
                                     double /*dInverseShift*/) const
{
  ASSERT_MSG( dynamic_cast<const SurfaceZero*>(&shiftedSurface), 
              "Surfaces do not have the same type");   

  // do trivial evalutation to remove compile warning C4100: non-reference
  shiftedSurface;

  return SurfacePtr( new SurfaceZero( GetDomain() ) );
}

} // namespace numeric

} // namespace ito33
