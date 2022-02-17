/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/heroinstdata.cpp
// Purpose:     Implementation of HG HeroInstData class
// Created:     2005/09/26
// RCS-ID:      $Id: heroinstdata.cpp,v 1.4 2006/01/22 22:02:14 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/src/hero/heroinstdata.cpp
   @brief Implementation of HG HeroInstData class
 */

#include "ito33/array.h"

#include "ito33/finance/payoff.h"

#include "hg/heroparams.h"
#include "hg/heromeshmanager.h"
#include "hg/heroinstdata.h"
#include "hg/herostepper.h"

namespace ito33
{

namespace hg
{

HeroInstData::HeroInstData(HeroParams& params,
                           Model& model, 
                           HeroMeshManager& meshes)
                         : BackwardInstData(params, model, meshes), 
                           m_heroParams(params), 
                           m_heroMeshes(meshes)
{
}

void HeroInstData::Init()
{
  // Get the space mesh
  m_pdS = m_heroMeshes.GetS();

  m_nNbS = m_heroMeshes.GetNbS();

  m_pdLogS = m_heroMeshes.GetLogS();

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);
 
}

void HeroInstData::SetInitialValue()
{
  double* pdPrices = m_pdPrices.Get();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    m_heroParams.GetPayoff()->Get(m_pdS, pdPrices + nIdxR * m_nNbS, m_nNbS);
 
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];
  }

  m_bHasEvent = false;
}

void HeroInstData::UpdateBeforeStep()
{
  BackwardInstData::UpdateBeforeStep();

  m_dRecoveryValue = m_heroMeshes.GetRecoveryValue();

  // Calculate the hero term: h * (V - V^T V(x,x) V)
  double dTime = m_heroMeshes.GetTime();

  const HedgingData& hedgingData = m_heroParams.GetHedgingData();

  m_pdHeroTerms = hedgingData.ComputePDETermsAt(dTime, m_pdS, m_nNbS);

  // Original equation is Dg^x + h * (V - V^T V(x,x) V) = 0
  double dH = m_heroMeshes.GetHValue();
  
  for (size_t nIdx = 0; nIdx < m_nNbRegimes * m_nNbS; nIdx++)
    m_pdHeroTerms[nIdx] *= dH;
}

} // namespace hg

} // namespace ito33
