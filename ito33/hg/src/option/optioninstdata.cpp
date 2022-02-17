/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/option/optioninstdata.cpp
// Purpose:     Implementation of HG OptionInstData class
// Created:     2005/01/13
// RCS-ID:      $Id: optioninstdata.cpp,v 1.14 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/src/option/optioninstdata.cpp
   @brief Implementation of HG OptionInstData class
 */

#include "ito33/sharedptr.h"

#include "ito33/finance/payoff.h"

#include "ito33/pricing/minconstraint.h"
#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/optionmeshmanager.h"

#include "hg/optioninstdata.h"

namespace ito33
{

namespace hg
{

  class Model;

OptionInstData::OptionInstData(pricing::OptionParams& params,
                               Model& model, 
                               pricing::OptionMeshManager& meshes)
                             : InstDataWithConstraints(params, model, meshes), 
                               m_optionParams(params), 
                               m_optionMeshes(meshes)
{
}

void OptionInstData::Init()
{
  // Get the space mesh
  m_pdS = m_optionMeshes.GetS();

  m_nNbS = m_optionMeshes.GetNbS();

  m_pdLogS = m_optionMeshes.GetLogS();

  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  m_optionMeshes.ComputeRecoveryValues(); 
  
  if (m_optionParams.GetOption().GetExerciseType() 
      == finance::ExerciseType_American)
  {
    double* pdPrices = m_pdPrices.Get();

    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_optionParams.GetPayoff()->Get(m_pdS, pdPrices + nIdxR * m_nNbS, m_nNbS);

    m_constraint.Update(pdPrices, m_nNbX);
    
    m_pConstraints = &m_constraint;
  }
  else
    m_pConstraints = 0;
}

void OptionInstData::SetInitialValue()
{
  double* pdPrices = m_pdPrices.Get();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    m_optionParams.GetPayoff()->Get(m_pdS, pdPrices + nIdxR * m_nNbS, m_nNbS);
 
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];
    m_piFrozenFlags[nIdx] = 0;
  }

  SetInitialSensitivityValue();

  m_bHasEvent = false;
}

void OptionInstData::UpdateBeforeStep()
{
  InstDataWithConstraints::UpdateBeforeStep();

  m_dRecoveryValue = m_optionMeshes.GetRecoveryValue();
}


} // namespace hg

} // namespace ito33
