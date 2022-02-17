/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/parbond/parbondpricer.cpp
// Purpose:     ParBond pricer class
// Created:     2005/06/09
// RCS-ID:      $Id: parbondpricer.cpp,v 1.3 2006/03/20 14:54:10 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/list.h"
#include "ito33/vector.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

#include "hg/model.h"
#include "hg/numoutputtimeonly.h"
#include "hg/steppertimeonly.h"
#include "hg/parbondpricer.h"

namespace ito33
{

namespace hg
{


typedef pricing::Engine
                 <
                   pricing::ParBondParams,
                   pricing::ParBondMeshManager,
                   ParBondInstData,
                   StepperTimeOnly, 
                   NumOutputTimeOnly
                 > ParBondEngine; 

ParBondPricer::ParBondPricer
    (pricing::ParBondParams& params, Model& model, 
    const finance::ComputationalFlags& flags)
   : m_params(params), 
     m_model(model),
     m_flags(flags),
     m_meshes(m_params, m_model),
     m_instdata(m_params, model, m_meshes),
     m_stepper(m_instdata),
     m_pNumOutput( new NumOutputTimeOnly(params) )
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  if ( flags.AreAllSensitivitiesActivated() )
    m_instdata.m_pbComputeSensitivities.resize
               (m_model.GetNbParameters(), true);
  else
    m_instdata.m_pbComputeSensitivities = flags.GetSensitivityFlags();

  m_pNumOutput->GetComputationalFlags().ActivateAllSensitivities
                                        ( flags.AreAllSensitivitiesActivated() );

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );
}

AutoPtr<NumOutputTimeOnly> ParBondPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMeTimeOnly();

  ParBondEngine
    engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);

  engine.Run();

  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33
