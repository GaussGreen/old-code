/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/heropricer.cpp
// Purpose:     HG HERO pricer class
// Created:     2005/09/26
// RCS-ID:      $Id: heropricer.cpp,v 1.5 2006/06/13 15:50:22 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/engine.h"

#include "ito33/finance/computationalflags.h"

#include "hg/model.h"
#include "hg/heronumoutput.h"
#include "hg/heropricer.h"
#include "hg/heroparams.h"

namespace ito33
{

namespace hg
{
  
typedef pricing::Engine< HeroParams, HeroMeshManager,
                         HeroInstData, HeroStepper, HeroNumOutput
                       > HeroEngine; 


HeroPricer::HeroPricer(HeroParams& params, 
                       Model& model,
                       const finance::ComputationalFlags& flags)
                     : m_params(params), 
                       m_model(model),
                       m_flags(flags),
                       m_meshes(m_params, m_model),
                       m_instdata(m_params, model, m_meshes),
                       m_stepper(m_instdata, m_flags),
                       m_pNumOutput(new HeroNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  // Let instdata setup its own flags
  m_instdata.SetupFlags(flags);

  // For option, Surface can be obtained for price and Greeks
  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );

  // save so we can get the hero
  //m_pNumOutput->SetFinalSave(true);
}

AutoPtr<HeroNumOutput> HeroPricer::Price()
{ 
  m_params.Init();

  m_meshes.SetupMe();
  
  HeroEngine
    runner(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);
  
  runner.Run();

  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33
