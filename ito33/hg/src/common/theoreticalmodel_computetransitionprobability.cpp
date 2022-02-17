/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_computetransitionprobability.cpp
// Purpose:     transition probability computation
// Created:     2006/03/28
// RCS-ID:      $Id: theoreticalmodel_computetransitionprobability.cpp,v 1.5 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/constants.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/transitionprobability.h"
#include "ito33/pricing/transitionprobabilityparams.h"
#include "ito33/pricing/transitionprobabilitymeshmanager.h"
#include "ito33/pricing/engine.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/transitionprobabilityoutput.h"

#include "hg/model.h"
#include "hg/transitionprobabilitynumoutput.h"
#include "hg/transitionprobabilitypricer.h"

namespace ito33
{

namespace hg
{
  
typedef pricing::Engine< pricing::TransitionProbabilityParams, 
                         pricing::TransitionProbabilityMeshManager,
                         TransitionProbabilityInstData, 
                         TransitionProbabilityStepper, 
                         TransitionProbabilityNumOutput
                       > TransitionProbabilityEngine; 

TransitionProbabilityPricer::TransitionProbabilityPricer
    (pricing::TransitionProbabilityParams& params, 
     Model& model,
     const finance::ComputationalFlags& flags)
   : m_params(params), 
     m_model(model),
     m_flags(flags),
     m_meshes(m_params, m_model),
     m_instdata(m_params, model, m_meshes),
     m_stepper(m_instdata, m_flags),
     m_pNumOutput(new TransitionProbabilityNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  m_instdata.SetupFlags(flags);

  // For option, Surface can be obtained for price and Greeks
  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );
}

shared_ptr<TransitionProbabilityOutput>
TransitionProbabilityPricer::Price()
{ 
  m_params.Init();

  m_meshes.SetupMe();

  TransitionProbabilityEngine
    runner(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);
  
  runner.Run();

  return m_pNumOutput->Solve(m_instdata);
}

shared_ptr<TransitionProbabilityOutput>
TheoreticalModel::ComputeTransitionProbability(double dPeriod)
{
  // use a faked spot/strike
  double dStrike = 100;

  pricing::TransitionProbability tp(dStrike, dPeriod);

  shared_ptr<numeric::NumParams> 
    pNumParams( new numeric::NumParams(*m_pQualityControl, dPeriod) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::TransitionProbabilityParams params(tp, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcess);

  // use default flags for transition probability
  finance::ComputationalFlags flags;

  TransitionProbabilityPricer pricer(params, modelParams, flags);
  
  return pricer.Price();
}

} // namespace hg

} // namespace ito33
