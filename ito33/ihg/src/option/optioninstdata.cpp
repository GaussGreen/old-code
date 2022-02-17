/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/option/optioninstdata.cpp
// Purpose:     Implementation of OptionInstData class
// Author:      David Pooley, WANG Xuewen
// Created:     2003/12/17
// RCS-ID:      $Id: optioninstdata.cpp,v 1.25 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/option/optioninstdata.cpp
    @brief Implementation of OptionInstData class
*/

#include "ito33/sharedptr.h"

#include "ito33/finance/payoff.h"

#include "ito33/pricing/minconstraint.h"
#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/optionmeshmanager.h"

#include "ihg/optioninstdata.h"

namespace ito33
{

namespace ihg
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

  Alloc(m_nNbS);

  m_optionMeshes.ComputeRecoveryValues(); 
  
  if ( m_optionParams.GetOption().GetExerciseType() 
      == finance::ExerciseType_American)
  {
    m_optionParams.GetPayoff()->Get(m_pdS, m_pdPrices.Get(), m_nNbS);

    m_constraint.Update(m_pdPrices.Get(), m_nNbS);
    
    m_pConstraints = &m_constraint;
  }
  else
    m_pConstraints = 0;
}

void OptionInstData::SetInitialValue()
{

  m_optionParams.GetPayoff()->Get(m_pdS, m_pdPrices.Get(), m_nNbS);
 
  if (m_bComputeVega)
    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdVegas[nIdxS] = 0.;

  m_bHasEvent = false;
}

void OptionInstData::UpdateBeforeStep()
{
  InstDataWithConstraints::UpdateBeforeStep();

  m_dRecoveryValue = m_optionMeshes.GetRecoveryValue();
}


} // namespace ihg

} // namespace ito33
