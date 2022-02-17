/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/surfacetimeonly.cpp
// Purpose:     Surface class with value depending only on time
// Author:      Wang
// Created:     2004/05/07
// RCS-ID:      $Id: surfacetimeonly.cpp,v 1.15 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 -   Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/numeric/surfacetimeonly.cpp
   @brief implementation of Surface class with value depending only on time
 */

#include "ito33/sharedptr.h"

#include "ito33/numeric/domain.h"
#include "ito33/numeric/surfacezero.h"
#include "ito33/numeric/surfacetimeonly.h"

namespace ito33
{

namespace numeric
{

typedef shared_ptr<SurfaceDouble> SurfacePtr;

void SurfaceTimeOnly::GetValuesAt
     (size_t nIdxT, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize( pdS.size(), GetValueAt(nIdxT));
}

void SurfaceTimeOnly::GetFirstValuesAt
     (size_t nIdxT, const Spots& pdS, Doubles& values) const
{
  values.clear(); 

  values.resize( pdS.size(), GetFirstValueAt(nIdxT) );
}

void SurfaceTimeOnly::GetLastValuesAt
     (size_t nIdxT, const Spots& pdS, Doubles& values) const
{
  values.clear(); 

  values.resize( pdS.size(), GetLastValueAt(nIdxT) );
}

void SurfaceTimeOnly::GetFirstValuesAt
     (double dTime, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize( pdS.size(), GetFirstValueAt(dTime) );
}

void SurfaceTimeOnly::GetLastValuesAt
      (double dTime, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize( pdS.size(), GetLastValueAt(dTime) );
}


// Not yet implemented
double SurfaceTimeOnly::GetFirstValueAt(double /* dTime */) const
{
  ASSERT_MSG(false, "Method not implemented in SurfaceTimeOnly"); 
  return 0.;
}


double SurfaceTimeOnly::GetLastValueAt(double /* dTime */) const
{
  ASSERT_MSG(false, "Method not implemented in SurfaceTimeOnly"); 
  return 0.;
}

void SurfaceTimeOnly::GetDeltaAndGamma
         ( SurfacePtr& pDeltaSurface, SurfacePtr& pGammaSurface ) const
{
  pDeltaSurface = SurfacePtr( new SurfaceZero(m_pDomain) );
  pGammaSurface = SurfacePtr( new SurfaceZero(m_pDomain) );
}

// Not implemented
void SurfaceTimeOnly::GetThetaBackwardOnly
                      ( SurfacePtr& pThetaSurface ) const
{
  double dOldTime = m_pDomain->GetTimeAt(m_ppdValues.size() - 2);
  double dTime;

  shared_ptr<SurfaceTimeOnly> 
    pSurfaceTmp( new SurfaceTimeOnly(m_pDomain) );

  pSurfaceTmp->Add(0.);

  for (size_t nIdxT = m_ppdValues.size() - 2; 
       nIdxT < m_ppdValues.size() - 1; 
       nIdxT--)
  {  
    dTime = m_pDomain->GetTimeAt(nIdxT + 1);

    size_t nSize = m_ppdValues[nIdxT].size();

    double dOldPrice = GetFirstValueAt(nIdxT);
    double dPrice = GetLastValueAt(nIdxT + 1);

    double dTheta = - (dPrice - dOldPrice) / (dTime - dOldTime);
      
    // there is an event
    if (nSize != 1)
      pSurfaceTmp->Add(0., 0.);
    else
      pSurfaceTmp->Add(dTheta);

    dOldTime = dTime;
  }

  pThetaSurface = pSurfaceTmp;
}

SurfacePtr
SurfaceTimeOnly::ComputeFiniteDifference(const SurfaceDouble& shiftedSurface,
                                         double dInverseShift) const
{
  ASSERT_MSG( dynamic_cast<const SurfaceTimeOnly*>(&shiftedSurface), 
              "Surfaces do not have the same type.");

  const SurfaceTimeOnly&
    shiftedTimeOnlySurface = 
      static_cast<const SurfaceTimeOnly&>(shiftedSurface);

  ASSERT_MSG
  ( 
    m_ppdValues.size() == shiftedTimeOnlySurface.m_ppdValues.size() &&
    m_ppdValues[0].size() == shiftedTimeOnlySurface.m_ppdValues[0].size(),
    "Surfaces do not have the same structure." 
  );

  // Create the return surface
  shared_ptr<SurfaceTimeOnly> pNewSurface( new SurfaceTimeOnly( GetDomain() ) );

  size_t
    nIdxTime,
    nNbValues,
    nIdxValue,
    nNbTimes = m_ppdValues.size();

  for( nIdxTime = 0; nIdxTime < nNbTimes; ++nIdxTime )
  {
    nNbValues = m_ppdValues[nIdxTime].size();
    std::vector<double> pdValues(nNbValues);

    for( nIdxValue = 0; nIdxValue < nNbValues; ++nIdxValue )
    {
      pdValues[nIdxValue] = 
        (shiftedTimeOnlySurface.m_ppdValues[nIdxTime][nIdxValue] -
         m_ppdValues[nIdxTime][nIdxValue]) * dInverseShift;
    }
    
    pNewSurface->m_ppdValues.push_back(pdValues);
  }
  
  return pNewSurface;

}

} // namespace numeric

} // namespace ito33

