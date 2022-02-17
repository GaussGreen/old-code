/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/eds/edspricer.cpp
// Purpose:     EDS pricer class
// Created:     2005/01/31
// RCS-ID:      $Id: edspricer.cpp,v 1.14 2006/06/13 15:49:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/edsparams.h"
#include "ito33/pricing/engine.h"

#include "hg/model.h"
#include "hg/backwardnumoutput.h"
#include "hg/edspricer.h"

#include "hg/stepper_fix.h"
#include "hg/stepper_fe_fix.h"

namespace ito33
{

namespace hg
{

EDSPricer::EDSPricer(pricing::EDSParams& params, 
                     Model& model, 
                     const finance::ComputationalFlags& flags)
                   : m_params(params), 
                     m_model(model),
                     m_flags(flags),
                     m_meshes(m_params, m_model),
                     m_instdata(m_params, model, m_meshes),
                     m_pNumOutput( new BackwardNumOutput(params) )
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  // Let instdata setup its own flags
  m_instdata.SetupFlags(flags);

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );

  m_pNumOutput->GetComputationalFlags().SetDiscretizationMethod
                                        ( flags.GetDiscretizationMethod() );
}

AutoPtr<BackwardNumOutput> EDSPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMe();
 
  if ( m_flags.GetDiscretizationMethod() )
  {
    typedef StepperFEFix EDSStepper;

    EDSStepper stepper(m_instdata, m_flags);

    pricing::Engine
        <
          pricing::Params,
          pricing::BackwardMeshManagerFix,
          BackwardInstData, 
          EDSStepper, 
          BackwardNumOutput
        > 
      engine(m_params, m_meshes, m_instdata, stepper, *m_pNumOutput);

    engine.Run();
  }
  else
  {
    typedef StepperFix EDSStepper;

    EDSStepper stepper(m_instdata, m_flags);

    pricing::Engine
        <
          pricing::Params,
          pricing::BackwardMeshManagerFix,
          BackwardInstData, 
          EDSStepper, 
          BackwardNumOutput
        > 
      engine(m_params, m_meshes, m_instdata, stepper, *m_pNumOutput);

    engine.Run();
  }

  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33
