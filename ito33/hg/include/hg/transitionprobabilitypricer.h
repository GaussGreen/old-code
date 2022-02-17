/////////////////////////////////////////////////////////////////////////////
// Name:        hg/transitionprobabilitypricer.h
// Purpose:     implementation of pricer class for transition probability
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobabilitypricer.h,v 1.2 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/transitionprobabilitypricer.h
    @brief implementation of pricer class transition probability
 */

#ifndef _HG_TRANSITIONPROBABILITYPRICER_H_ 
#define _HG_TRANSITIONPROBABILITYPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/optionmeshmanager.h"

#include "hg/transitionprobabilityinstdata.h"
#include "hg/stepper_fix.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace pricing
{
  class TransitionPrpbabilityParams;
}

namespace hg
{
  class Model;
  class TransitionProbabilityNumOutput;
  class TransitionProbabilityOutput;

  typedef StepperFix TransitionProbabilityStepper;

/// Transition Probability Pricer class.
class TransitionProbabilityPricer
{
public:

  /// Construct from transition probability parameters
  TransitionProbabilityPricer(pricing::TransitionProbabilityParams& params, 
                              Model& model,
                              const finance::ComputationalFlags& flags);
  
  /// Computes the transition probability
  shared_ptr<TransitionProbabilityOutput> Price();

private:

  pricing::TransitionProbabilityParams& m_params;

  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::TransitionProbabilityMeshManager m_meshes;
  
  TransitionProbabilityInstData m_instdata;
  
  TransitionProbabilityStepper m_stepper;
  
  AutoPtr<TransitionProbabilityNumOutput> m_pNumOutput;

  NO_COPY_CLASS(TransitionProbabilityPricer);

}; // class TransitionProbabilityPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_TRANSITIONPROBABILITYPRICER_H_
