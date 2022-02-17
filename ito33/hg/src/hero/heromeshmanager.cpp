/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/heromeshmanager.cpp
// Purpose:     option mesh manager for backward PDE problems
// Created:     2004/02/11
// RCS-ID:      $Id: heromeshmanager.cpp,v 1.3 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/optionliketype.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "hg/model.h"
#include "hg/heroparams.h"
#include "hg/heromeshmanager.h"

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

namespace ito33
{

namespace hg
{

HeroMeshManager::HeroMeshManager(HeroParams& params, hg::Model& model)
                               : BackwardMeshManagerFix(params, model),
                                 m_heroParams(params), m_HGmodel(model) 
{ 
}

void HeroMeshManager::ConstructSpaceMesh()
{
  if (m_params.GetMeshParams()->GetUniformSpaceGrid() == true)
  {
    ConstructUniformSpaceMesh();
    return;
  }

  // TODO: Make a better mesh?  How?  Include special points from
  // all hedge contracts?
  // For lack of anything better at the moment, use an option mesh
  // centered on the spot (instead of the strike).
  OptionSpaceMesh osgGrid;

  osgGrid.SetOptionType( pricing::Option_Call );

  double dMaturityTime = m_heroParams.GetHedgingData().GetMaturityTime();
  double dSpot = m_params.GetSpotSharePrice();

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = m_model.GetSquaredTotalVol(dMaturityTime, dSpot);

  osgGrid.SetDiffusionSize
          ( m_params.GetDiffusionSize(dSquaredTotalVol) );

  // Compute the convection size for this problem
  double dConvection = m_model.GetConvection( dMaturityTime, dSpot );

  osgGrid.SetConvectionSize
          ( m_params.GetBackwardConvectionSize(dConvection) );

  
  // Generate a space mesh centered at zero
  std::vector<double> vecGrid;
  
  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_params.GetNumParams()->GetNbSpaceSteps();

  osgGrid.Build(nNbPoints, 0.0, vecGrid);

  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);

  // Center at the spot price
  double dLogSpot = log( dSpot );

  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx] + dLogSpot;
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }

}


void HeroMeshManager::ConstructUniformSpaceMesh()
{

  // Force the points 0.0 and the spot.  Continue until approximately
  // three times the spot.  This means about 1/3 of the points will be
  // between 0 and the strike.

  m_nNbS = m_params.GetNumParams()->GetNbSpaceSteps();
  size_t nNbTmp = m_nNbS / 3;
  double dDeltaT = m_heroParams.GetSpotSharePrice() / nNbTmp;
  
  // Construct the actual mesh
  m_pdS = Array<double>(m_nNbS);
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS[0] = 0.0;
  m_pdLogS[0] = log(dDeltaT / 100.0);
  for (size_t nIdx = 1; nIdx < m_nNbS; nIdx++)
  {    
    m_pdS[nIdx] = nIdx * dDeltaT;
    m_pdLogS[nIdx] = log(m_pdS[nIdx]);
  }

}


void HeroMeshManager::ComputeRecoveryValues()
{
  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
    m_pdRecoveryValues[nIdxT] = 0.;
}

void HeroMeshManager::ComputeHValues()
{
  // Pre-compute the 'h' values: h(t) = exp( \int (2 r - SR^2) ds )
  // Re-write as h(t) = exp( 2 \int r ds  -  \int SR^2 ds)
  //                  = exp(\int r ds)^2 * 1 / exp( \int SR^2 ds)
  //                  = exp(\int r ds)^2 * 1/ exp( (T-t) * SR^2 )
  m_pdHValues = Array<double>(m_nNbTimes);

  Array<double> pdCF(m_nNbTimes);
     
  m_heroParams.GetYieldCurve()->GetCompoundFactor
    (m_pdTimes.Get(), pdCF.Get(), m_nNbTimes);

  double dSR = m_HGmodel.GetSharpeRatio();
  double dSR2 = dSR * dSR;
  
  for (size_t nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
  {    
    double dTmp = pdCF[m_nNbTimes-1] / pdCF[nIdxT];
    double dFact1 = dTmp * dTmp;

    double dLength = m_pdTimes[m_nNbTimes-1] - m_pdTimes[nIdxT];
    double dFact2 = exp(dLength * dSR2);

    m_pdHValues[nIdxT] = dFact1 / dFact2;
  }


}

void HeroMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  ComputeRecoveryValues();

  ComputeHValues();
}

} // namespace hg

} // namespace ito33
