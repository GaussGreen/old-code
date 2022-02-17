/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/onetouch/onetouchpricer.cpp
// Purpose:     OneTouch pricer class
// Created:     2006/08/11
// RCS-ID:      $Id: onetouchpricer.cpp,v 1.1 2006/08/10 23:12:02 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/onetouchparams.h"

#include "ihg/model.h"
#include "ihg/onetouchnumoutput.h"
#include "ihg/onetouchpricer.h"

namespace ito33
{

namespace ihg
{

typedef pricing::Engine
                 <
                   pricing::OneTouchParams,
                   pricing::OneTouchMeshManager,
                   OneTouchInstData, 
                   OneTouchStepper, 
                   OneTouchNumOutput
                 > OneTouchEngine; 

OneTouchPricer::OneTouchPricer(pricing::OneTouchParams& params, 
                               Model& model, 
                               const finance::ComputationalFlags& flags)
                             : m_params(params), 
                               m_model(model),
                               m_flags(flags),
                               m_meshes(m_params, m_model),
                               m_instdata(m_params, model, m_meshes),
                               m_stepper(m_instdata, m_flags),
                               m_pNumOutput(new OneTouchNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  m_instdata.m_bComputeVega = flags.GetComputeVega();

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                          (flags.GetComputeSurface());
}

AutoPtr<OneTouchNumOutput> OneTouchPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMe();

  OneTouchEngine
    engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);

  engine.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
