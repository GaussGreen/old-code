/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/transitionprobabilityparams.h
// Purpose:     params class for transition probability
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobabilityparams.h,v 1.2 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/transitionprobabilityparams.h
    @brief Implementation of the params class for transition probability.
 */

#ifndef _ITO33_PRICING_TRANSITIONPROBABILITYPARAMS_H_
#define _ITO33_PRICING_TRANSITIONPROBABILITYPARAMS_H_

#include "ito33/pricing/transitionprobability.h"
#include "ito33/pricing/params.h"

namespace ito33 
{
 
namespace pricing
{

/// Transition probability params class
class TransitionProbabilityParams : public Params
{

public: 

  /**
      Ctor initializes a param suitable to compute the time homogeneous
      transition probability: no yieldcurve, dividends, borrow curve.
   */
  TransitionProbabilityParams
  (TransitionProbability& tp,    
   const shared_ptr<numeric::NumParams>& pNumParams,
   const shared_ptr<numeric::MeshParams>& pMeshParams);

  // Default dtor is ok

  void Init();

  /// Gets a reference to the transition probability contract
  TransitionProbability& GetTransitionProbability() const { return m_tp; }


private:
  
  TransitionProbability& m_tp; 

  NO_COPY_CLASS(TransitionProbabilityParams);

}; // class TransitionProbabilityParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_TRANSITIONPROBABILITYPARAMS_H_
