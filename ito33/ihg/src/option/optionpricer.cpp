/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/option/optionpricer.h
// Purpose:     option pricer class
// Author:      David
// Created:     2003/12/19
// RCS-ID:      $Id: optionpricer.cpp,v 1.25 2006/06/13 15:42:54 wang Exp $
// Copyright:   (c) 2003 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/optionparams.h"

#include "ihg/model.h"
#include "ihg/optionnumoutput.h"
#include "ihg/optionpricer.h"

namespace ito33
{

namespace ihg
{
  
typedef pricing::Engine< pricing::OptionParams, pricing::OptionMeshManager,
                         OptionInstData, OptionStepper, OptionNumOutput
                       > OptionEngine; 


OptionPricer::OptionPricer(pricing::OptionParams& params, 
                           Model& model,
                           const finance::ComputationalFlags& flags)
                         : m_params(params), 
                           m_model(model),
                           m_flags(flags),
                           m_meshes(m_params, m_model),
                           m_instdata(m_params, model, m_meshes),
                           m_stepper(m_instdata, flags),
                           m_pNumOutput(new OptionNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  // For option, Vega can be computed by PDE. 
  m_instdata.m_bComputeVega = flags.GetComputeVega();

  // For option, Surface can be obtained for price and Greeks
  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );
}

ito33::AutoPtr<OptionNumOutput> OptionPricer::Price()
{ 
  m_params.Init();

  m_meshes.SetupMe();
  
  OptionEngine
    runner(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);
  
  runner.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33

