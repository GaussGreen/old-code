/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cds/cdspricer.cpp
// Purpose:     CDS pricer class
// Created:     2005/02/16
// RCS-ID:      $Id: cdspricer.cpp,v 1.15 2006/08/21 15:29:20 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/cdsparams.h"

#include "hg/model.h"
#include "hg/numoutputtimeonly.h"
#include "hg/steppertimeonly.h"
#include "hg/cdspricer.h"

namespace ito33
{

namespace hg
{


typedef pricing::Engine
                 <
                   pricing::CDSParams,
                   pricing::CDSMeshManager,
                   CDSInstData,
                   StepperTimeOnly, 
                   NumOutputTimeOnly
                 > CDSEngine; 

CDSPricer::CDSPricer(pricing::CDSParams& params, 
                     Model& model, 
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

AutoPtr<NumOutputTimeOnly> CDSPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMeTimeOnly();

  CDSEngine
    engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);

  engine.Run();

  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33
