/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/option/optionpricer.cpp
// Purpose:     HG option pricer class
// Created:     2005/01/13
// RCS-ID:      $Id: optionpricer.cpp,v 1.16 2006/06/13 15:49:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/optionparams.h"

#include "ito33/finance/computationalflags.h"

#include "hg/model.h"
#include "hg/optionnumoutput.h"
#include "hg/optionpricer.h"
#include "hg/optionstepper.h"

namespace ito33
{

namespace hg
{

OptionPricer::OptionPricer(pricing::OptionParams& params, 
                           Model& model,
                           const finance::ComputationalFlags& flags)
                         : m_params(params), 
                           m_model(model),
                           m_flags(flags),
                           m_meshes(m_params, m_model),
                           m_instdata(m_params, model, m_meshes),
                           m_pNumOutput(new OptionNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  // Let instdata setup its own flags
  m_instdata.SetupFlags(flags);

  // For option, Surface can be obtained for price and Greeks
  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );

  m_pNumOutput->GetComputationalFlags().SetDiscretizationMethod
                                        ( flags.GetDiscretizationMethod() );
}

AutoPtr<OptionNumOutput> OptionPricer::Price()
{ 
  m_params.Init();

  m_meshes.SetupMe();
  
  if ( m_flags.GetDiscretizationMethod() )
  {
    typedef OptionStepper<StepperFEFix> MyStepper;

    MyStepper stepper(m_instdata, m_flags);

    pricing::Engine
        <
          pricing::Params,
          pricing::BackwardMeshManagerFix,
          OptionInstData, 
          MyStepper, 
          OptionNumOutput
        > 
      engine(m_params, m_meshes, m_instdata, stepper, *m_pNumOutput);

    engine.Run();
  }
  else
  {
    typedef OptionStepper<StepperFix> MyStepper;

    MyStepper stepper(m_instdata, m_flags);

    pricing::Engine
        <
          pricing::Params,
          pricing::BackwardMeshManagerFix,
          OptionInstData, 
          MyStepper, 
          OptionNumOutput
        > 
      engine(m_params, m_meshes, m_instdata, stepper, *m_pNumOutput);

    engine.Run();
  }

  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33
