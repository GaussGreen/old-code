/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/eds/edspricer.cpp
// Purpose:     EDS pricer class
// Created:     2005/01/26
// RCS-ID:      $Id: edspricer.cpp,v 1.5 2006/06/13 15:43:40 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/edsparams.h"

#include "ihg/model.h"
#include "ihg/edsnumoutput.h"
#include "ihg/edspricer.h"

namespace ito33
{

namespace ihg
{

typedef pricing::Engine
                 <
                   pricing::EDSParams,
                   pricing::EDSMeshManager,
                   EDSInstData, 
                   EDSStepper, 
                   EDSNumOutput
                 > EDSEngine; 

EDSPricer::EDSPricer(pricing::EDSParams& params, 
                     Model& model, 
                     const finance::ComputationalFlags& flags)
                   : m_params(params), 
                     m_model(model),
                     m_flags(flags),
                     m_meshes(m_params, m_model),
                     m_instdata(m_params, model, m_meshes),
                     m_stepper(m_instdata, m_flags),
                     m_pNumOutput(new EDSNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  m_instdata.m_bComputeVega = flags.GetComputeVega();

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                          (flags.GetComputeSurface());
}

AutoPtr<EDSNumOutput> EDSPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMe();

  EDSEngine
    engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);

  engine.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
