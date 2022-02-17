/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/meshmanager.cpp
// Purpose:     base mesh manager for PDE problems
// Author:      David Pooley
// Created:     2004/01/05
// RCS-ID:      $Id: meshmanager.cpp,v 1.30 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <list>
#include <vector>
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/debug.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/mesh/genraf.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/model.h"

using namespace ito33;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

void MeshManager::SetInitialState()
{
  m_params.SetInitialState( GetTime() );

  m_model.SetInitialState( GetTime() );
}

void MeshManager::GetInstValues(double& dTimeStep,
                                double& dRate, double& dDerivativeRate,
                                double& dForeignRate,
                                SchemeType& schemeType) const
{
  dTimeStep = GetTimeStep();

  dRate = m_pdRates[m_nIdx];

  if( m_pdDerivativeRates.Get() )
    dDerivativeRate = m_pdDerivativeRates[m_nIdx];
  else
    dDerivativeRate = dRate;

  dForeignRate = m_pdForeignRates[m_nIdx];

  schemeType = m_pSchemeTypes[m_nIdx];
}

/// Try if another timestep can be made
bool MeshManager::TryGoAhead()
{
  if ( CanGoAhead() )
  {
    GoAhead();

    m_params.Update( GetTime() );

    m_model.Update( GetTime() );

    return true;
  }

  return false;
}

void MeshManager::SetupMe()
{
  // @todo, check that the special times are in the time mesh
  if (m_bUseOutsideTimeMesh)
    return;

  // Get the special times of contract for the construction of the time mesh
  SpecialTimes specialTimesAll;

  m_params.GetSpecialTimes(specialTimesAll);
  
  double dValuationTime = m_params.GetValuationTime();
  double dStoppingTime = m_params.GetStoppingTime();

  // Get the special times of model for the construction of the time mesh
  SpecialTimes specialTimesTmp;

  m_model.GetSpecialTimes(specialTimesTmp);

  // Add the useful special times in the model to the list
  SpecialTimes::const_iterator i;

  for (i = specialTimesTmp.begin(); i != specialTimesTmp.end(); ++i)
    if (i->GetTime() >= dValuationTime && i->GetTime() < dStoppingTime)
      specialTimesAll.push_back(*i);

  // Sort the list
  specialTimesAll.sort();

  // Remove duplicated times, manually written, can probably be improved 
  
  SpecialTimes::iterator iter1 = specialTimesAll.begin();
  
  // the final, well sorted, duplicated points removed list
  SpecialTimes specialTimes;
  
  specialTimes.push_back(*iter1++);

  for ( ; iter1 != specialTimesAll.end(); ++iter1)
    if ( *iter1 > specialTimes.back() ) // new entry
      specialTimes.push_back(*iter1);
    else // use highest RefineLevel among SpecialTimes at given time 
    { 
      if (iter1->GetRefineLevel() > specialTimes.back().GetRefineLevel())
      {
        specialTimes.pop_back();
        specialTimes.push_back(*iter1);
      }
    }

  // Construct the time mesh
  ConstructTimeMesh(specialTimes);

  ASSERT_MSG(m_nNbTimes > 2, "More than 2 time points are required");
  
  /* 
    Set up the scheme type for each time step. 

    Whenever a event is encountered, we need to use an implicite scheme
  */
  SetupSchemeTypes(specialTimes);

  /*
    Set up the rates so that constant will be an exact solution for the 
    discretized PDE. 
    
    It needs that scheme types are already determined for the time steps.
  */
  SetupRates();
}


void MeshManager::SetupTimeMesh(const double* pdTimeMesh, 
                                const SchemeType* pSchemeTypes, 
                                size_t nNbTimesteps)
{

  // Construct the time mesh by copying.  Also copy the scheme types
  m_nNbTimes = nNbTimesteps;
  m_pdTimes = Array<double>(m_nNbTimes);
  m_pSchemeTypes = Array<SchemeType>(m_nNbTimes);
  for (size_t nIdx = 0; nIdx < m_nNbTimes; nIdx++)
  {
    m_pdTimes[nIdx] = pdTimeMesh[nIdx];
    m_pSchemeTypes[nIdx] = pSchemeTypes[nIdx];
  }

  ASSERT_MSG(m_nNbTimes >= 2, "At least 2 time points are required");
  
  /*
    Set up the rates so that constant will be an exact solution for the
    discretized PDE. 
    
    It needs that scheme types are already determined for the time steps.
  */
  SetupRates();

  m_bUseOutsideTimeMesh = true;
}

void MeshManager::SetupTimeMesh(SpecialTimes& specialTimes)
{
  // Only construct the time mesh and scheme types. Assume that
  // this mesh is being used to construct a master list for path
  // dependent solving
  ConstructTimeMesh(specialTimes);

  SetupSchemeTypes(specialTimes);

  SetupRates();

  m_bUseOutsideTimeMesh = true;
}

void MeshManager::ConstructUniformTimeMesh(SpecialTimes& specialTimes)
{
  // Try to construct a uniform time mesh. It may be impossible if special
  // times are required. Also, since we don't know the final mesh size
  // (especially if two special times are closer than the average stepsize), 
  // construct a temp mesh and then copy at the end.
  size_t nNbRequestedTimesteps = m_params.GetNumParams()->GetNbTimeSteps();
  if ( nNbRequestedTimesteps < specialTimes.size() )
    nNbRequestedTimesteps = specialTimes.size();
  Array<double> pdTempMesh(nNbRequestedTimesteps * 10);

  // Start a counter for the actual number of timesteps 
  size_t nNbTimesteps = 0;

  // Iterate over the special points, constructing the mesh in intervals
  SpecialTimes::iterator iterTimes = specialTimes.begin();
  double dLeftTime = iterTimes->GetTime();
  double dTotalTime = specialTimes.back().GetTime()
                    - specialTimes.front().GetTime();

  // Add the first point to start the process
  pdTempMesh[nNbTimesteps] = dLeftTime;
  nNbTimesteps++;

  for (++iterTimes; iterTimes != specialTimes.end(); ++iterTimes)
  {
    // Determine how many points should be used in the current interval
    // If in the last interval, force the number of points to add to the 
    // requested number of points
    double dRightTime = iterTimes->GetTime();
    size_t nNbTimestepsInterval;
    if ( dRightTime == specialTimes.back().GetTime() )
      nNbTimestepsInterval = nNbRequestedTimesteps - nNbTimesteps + 1;
    else
    {
      double dRatio = (dRightTime - dLeftTime) / dTotalTime;
      nNbTimestepsInterval = (size_t) ceil(dRatio * nNbRequestedTimesteps);
    }

    if (nNbTimestepsInterval < 2)
      nNbTimestepsInterval = 2;

    // Get the timestep size, and construct this part of the mesh
    // Do not include the first (left) point which has already been added as the
    // right point of the previous interval. Recall that with n points, there
    // are n-1 intervals
    double dDeltaT = (dRightTime - dLeftTime) / (nNbTimestepsInterval - 1);
    for (size_t nIdx = 1; nIdx < nNbTimestepsInterval; nIdx++)
    {
      pdTempMesh[nNbTimesteps] = pdTempMesh[nNbTimesteps - 1] + dDeltaT;
      nNbTimesteps++;
    }
    //pdTempMesh[nNbTimesteps-1] = dRightTime;

    // Get ready for next interval
    dLeftTime = dRightTime;
  }
  
  // Construct the real mesh by copying
  m_nNbTimes = nNbTimesteps;
  m_pdTimes = Array<double>(m_nNbTimes);
  for (size_t nIdx = 0; nIdx < m_nNbTimes; nIdx++)
  {
    m_pdTimes[nIdx] = pdTempMesh[nIdx];
  }
}

void GeneralTimeMesh(SpecialTimes& specialTimes,
                     int iDirection,
                     size_t nMinNbTimeSteps,
                     std::vector<double>& pdTimes)
{
  size_t nNbEvents = specialTimes.size();
  
  ASSERT_MSG(nNbEvents >= 2,
             "pdSpecialTimes's dimension should be greater than 1.");

  std::vector<double> pdEvents(nNbEvents);
  std::vector<RefineLevel> pRefineLevels(nNbEvents);

  size_t nIdxEvent;

  // Time points are reversed for backward problem to use the same algorithm
  // as for forward problem. 
  if (iDirection)
  {
    SpecialTimes::iterator iterTimes;
    for (iterTimes = specialTimes.begin(), nIdxEvent = 0;
         iterTimes != specialTimes.end();
         ++iterTimes, nIdxEvent++)
    {
      pdEvents[nIdxEvent] = iterTimes->GetTime();
      pRefineLevels[nIdxEvent] = iterTimes->GetRefineLevel();
    }
  }
  else
  {
    SpecialTimes::reverse_iterator iterTimes;
    for (iterTimes = specialTimes.rbegin(), nIdxEvent = 0;
         iterTimes != specialTimes.rend();
         ++iterTimes, nIdxEvent++)
    {
      pdEvents[nIdxEvent] = - iterTimes->GetTime();
      pRefineLevels[nIdxEvent] = iterTimes->GetRefineLevel();
    }
  }
 
  // primary datas of the time mesh
  double dTotalTime = pdEvents[nNbEvents - 1] - pdEvents[0];
  double dRequestedMaxDeltaT = dTotalTime / nMinNbTimeSteps;

  // Genraf numbers for VeryHigh refine lvel
  const double dZoomStart = 20;
  const double dRhoStart = 1.5;
  
  // GenRaf numbers for High refine level
  const double dZoomMedium = 10; 
  const double dRhoMedium = 1.7;

  // GenRaf numbers for standard refine level
  const double dZoomNormal = 2;
  const double dRhoNormal = 2.;

  // Estimate the size of the to be generated time mesh 
  const size_t nNbJoinWorst = int(log(dZoomStart) / log(dRhoStart)) + 1;
  size_t nNbTimes = int(dTotalTime / dRequestedMaxDeltaT) + 2 * nNbEvents
                  + nNbEvents * nNbJoinWorst;
  
  // allocate the momory for time array according to the estimated size 
  pdTimes.resize(nNbTimes); 
  
  // Generate the time mesh for the first period. Can be merged?
  double dDeltaStart = dRequestedMaxDeltaT / dZoomStart;
  int iNbTmp;

  GenRaf(pdEvents[0], pdEvents[1], 
         dDeltaStart, dRequestedMaxDeltaT, dRhoStart,
         1, 
         &(pdTimes[0]), iNbTmp);

  nNbTimes = size_t(iNbTmp);

  // Generate the time mesh for the other periods.
  for (nIdxEvent = 1; nIdxEvent < nNbEvents - 1; nIdxEvent++)
  {
    // Translate RefineLevel into numbers required by GenRaf
    double dZoom, dRho;
    switch (pRefineLevels[nIdxEvent])
    {
    case RefineLevel_High:
      dZoom = dZoomMedium;
      dRho = dRhoMedium;
      break;
    case RefineLevel_VeryHigh:
      dZoom = dZoomStart;
      dRho = dRhoStart;
      break;  
    default:
      dZoom = dZoomNormal;
      dRho = dRhoNormal;
      break;
    }

    double dDelta = dRequestedMaxDeltaT / dZoom;

    double dRealNormalDelta = pdTimes[nNbTimes - 1] - pdTimes[nNbTimes - 2];
    if (dDelta >= dRealNormalDelta)
      dDelta = dRealNormalDelta;
    
    // Generate the time mesh between two events
    GenRaf(pdEvents[nIdxEvent], 
           pdEvents[nIdxEvent + 1],
           dDelta,
           dRequestedMaxDeltaT,
           dRho,
           1, 
           &(pdTimes[nNbTimes - 1]),
           iNbTmp);
    
    nNbTimes += iNbTmp - 1;
  }
  
  // change the vector to the real size
  pdTimes.resize(nNbTimes);

  if (!iDirection) // change the sign back for backward problem
  { 
    size_t nIdxTmp = nNbTimes - 1;
    double dTmp;
    for (nIdxEvent = 0; nIdxEvent < nIdxTmp; nIdxEvent++, nIdxTmp--)
    {
      dTmp = pdTimes[nIdxEvent];
      pdTimes[nIdxEvent] = -pdTimes[nIdxTmp];
      pdTimes[nIdxTmp] = -dTmp;
    }
    if (nIdxEvent == nIdxTmp)
      pdTimes[nIdxEvent] = - pdTimes[nIdxEvent];
  }

  // the first and last points can get slightly modified, change them back.
  // should do the same thing for other special times? or it's not necessary 
  // at all even for these two points since we use always a tolerance to 
  // compare two time points? 
  pdTimes[0] = specialTimes.front().GetTime();
  pdTimes[nNbTimes - 1] = specialTimes.back().GetTime();
}

void MeshManager::ConstructNonuniformTimeMesh
                  (SpecialTimes& specialTimes, int iDirection)
{
  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes,
                  iDirection,
                  GetNbRequestedTimes(),
                  pdTimes);

  // Construct the real mesh by copying
  m_nNbTimes = pdTimes.size();

  m_pdTimes = Array<double>(m_nNbTimes);
  for (size_t nIdx = 0; nIdx < m_nNbTimes; nIdx++)
    m_pdTimes[nIdx] = pdTimes[nIdx];
}

size_t MeshManager::GetNbRequestedTimes() const
{
  return m_params.GetNumParams()->GetNbTimeSteps();
}

