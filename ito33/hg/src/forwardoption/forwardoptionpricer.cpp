/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/forward/forwardoptionpricer.h
// Purpose:     HG forward option pricer class
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionpricer.cpp,v 1.12 2006/06/13 16:25:07 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/forwardoptionparams.h"

#include "ito33/finance/computationalflags.h"

#include "hg/model.h"
#include "hg/forwardoptionnumoutput.h"
#include "hg/forwardoptionpricer.h"
#include "hg/forwardoptionstepper.h"

namespace ito33
{

namespace hg
{

ForwardOptionPricer::ForwardOptionPricer(pricing::ForwardOptionParams& params, 
                           Model& model,
                           const finance::ComputationalFlags& flags)
                         : m_params(params), 
                           m_model(model),
                           m_flags(flags),
                           m_meshes(m_params, m_model),
                           m_instdata(m_params, model, m_meshes),
                           m_pNumOutput(new ForwardOptionNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  // Let instdata setup its own flags
  m_instdata.SetupFlags(flags);

  // Adjoint method not yet implemented for finite difference

  int iDiscretizationMethod = flags.GetDiscretizationMethod();

  if ( !iDiscretizationMethod )
  {
    // Should check also dividend!
    if ( m_instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
      m_instdata.m_sensitivityMethod = SensitivityMethod_PDE;
  }

  // For forward option, surface can be obtained for price and Greeks
  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );

  m_pNumOutput->GetComputationalFlags().SetDiscretizationMethod
                                        ( iDiscretizationMethod );
}

AutoPtr<ForwardOptionNumOutput> ForwardOptionPricer::Price()
{ 
  m_params.Init();

  m_meshes.SetupMe();
  
  if ( m_flags.GetDiscretizationMethod() )
  {
    typedef ForwardOptionStepperFE MyStepper;

    MyStepper stepper(m_instdata, m_flags);

    pricing::Engine
        <
          pricing::ForwardOptionParams,
          pricing::ForwardOptionMeshManager,
          ForwardOptionInstData, 
          MyStepper, 
          ForwardOptionNumOutput
        > 
      engine(m_params, m_meshes, m_instdata, stepper, *m_pNumOutput);

    engine.Run();
  }
  else
  {
    typedef ForwardOptionStepperFD MyStepper;

    MyStepper stepper(m_instdata, m_flags);

    pricing::Engine
        <
          pricing::ForwardOptionParams,
          pricing::ForwardOptionMeshManager,
          ForwardOptionInstData, 
          MyStepper, 
          ForwardOptionNumOutput
        > 
      engine(m_params, m_meshes, m_instdata, stepper, *m_pNumOutput);

    engine.Run();
  }

  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33

