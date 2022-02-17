/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/onetouchmeshmanager.cpp
// Purpose:     OneTouch mesh manager for backward PDE problems
// Created:     2005/07/04
// RCS-ID:      $Id: onetouchmeshmanager.cpp,v 1.5 2006/02/23 13:12:18 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/genraf.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/onetouchparams.h"
#include "ito33/pricing/onetouchmeshmanager.h"

namespace ito33
{

namespace pricing
{
  using namespace numeric;
  using namespace numeric::mesh;


OneTouchMeshManager::OneTouchMeshManager
(OneTouchParams& params, Model& model)
  : BackwardMeshManagerFix(params, model), m_oneTouchParams(params)    
{ 
}

// todo: all these magic numbers should be documented and shared with others
void OneTouchMeshManager::ConstructSpaceMesh()
{  
  const double dSpot = m_params.GetSpotSharePrice();
  const MeshParams* pMeshParams = m_params.GetMeshParams();

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
                            ( m_params.GetStoppingTime(), dSpot );

  double dDiffusion = pMeshParams->GetGridSpan()
                    * m_params.GetDiffusionSize(dSquaredTotalVol);

  if ( dDiffusion < 0.6 )
    dDiffusion = 0.6;

  double dConvectionSize = m_params.GetBackwardConvectionSize(dSquaredTotalVol);
  if ( dConvectionSize < 0.)
    dConvectionSize = 0.;

  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_params.GetNumParams()->GetNbSpaceSteps();
    
  double dA = log( m_oneTouchParams.GetOneTouch().GetBarrier() );

  double dB = dA + dDiffusion + dConvectionSize;
  double dTmp = log(dSpot * 2.0);
  if ( dB < dTmp )
    dB = dTmp;

  double dMax = 10.;
 
  double dDelta = (dB - dA) / nNbPoints;
  
  double dDeltaMin = dDelta / pMeshParams->GetSpaceAccumulation();
  double dStretch = pMeshParams->GetSpaceStretch();
  double dDeltaMax = 0.5;

  // Use an upper bound for the normal delta, otherwise the result can be very
  // bad
  if ( dDelta > 0.15 )
    dDelta = 0.15;

  m_nNbS = GenRafSize(dB - dA, dDeltaMin, dDelta, dStretch);

  int iNbS2 = 1;
  
  if ( dB < dMax )
    iNbS2 = GenRafSize( dMax - dB, dDelta, dDeltaMax, dStretch);
  
  m_pdLogS = Array<double>(m_nNbS);
  
  int iNbS;

  // Generate the mesh from dA to dB and refine around the center which for now
  // is just the barrier
  GenRaf(dA, dB, dDeltaMin, dDelta, dStretch, 1, m_pdLogS.Get(), iNbS);

  // Generate the mesh from dB to a large number, beginning with previous 
  // mesh spacing, and growing to a max delta.
  //dDelta = m_pdLogS[iNbS-1] - m_pdLogS[iNbS-2];
  if ( dB < dMax )
    GenRaf(dB, dMax, dDelta, dDeltaMax, dStretch, 1, 
           m_pdLogS.Get() + iNbS - 1, iNbS2);

  m_nNbS = iNbS + iNbS2 - 1;
  m_pdS = Array<double>(m_nNbS);
 
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);

  // Revert the mesh accordingly if the barrier is down
  if (    m_oneTouchParams.GetOneTouch().GetBarrierType()
       == finance::Barrier_UpAndOut )
  {
    Array<double> pdLogSTmp(m_nNbS);
    Array<double> pdSTmp(m_nNbS);
    for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    {
      pdLogSTmp[nIdx] = 2. * m_pdLogS[0] - m_pdLogS[m_nNbS - 1 - nIdx];
      pdSTmp[nIdx] = exp(pdLogSTmp[nIdx]);
    }

    swap(m_pdLogS, pdLogSTmp);
    swap(m_pdS, pdSTmp);
  }
}

void OneTouchMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  OneTouch& oneTouch = m_oneTouchParams.GetOneTouch();

  m_pdBoundaryValues = Array<double>(m_nNbTimes);

  // Precompute the boundary value
  if ( oneTouch.IsImmediateRebate() )
    for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
      m_pdBoundaryValues[nIdxT] = 1.;
  else
  {
    Array<double> pdCF(m_nNbTimes);
    m_params.GetYieldCurve()->GetCompoundFactor
                              ( m_pdTimes.Get(), pdCF.Get(), m_nNbTimes );

    double dTmp = 1. / pdCF[m_nNbTimes - 1];
    for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
      m_pdBoundaryValues[nIdxT] = pdCF[nIdxT] * dTmp;
  }

  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);
  if ( oneTouch.GetBarrierType() == finance::Barrier_UpAndOut )
    for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
      m_pdRecoveryValues[nIdxT] = 0.;
  else
    for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
       m_pdRecoveryValues[nIdxT] = m_pdBoundaryValues[nIdxT];
}


} // namespace pricing

} // namespace ito33
