/////////////////////////////////////////////////////////////////////////////
// Name:        hg/transitionprobabilityinstdata.h
// Purpose:     instdata class for HG transition probability computation
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobabilityinstdata.h,v 1.1 2006/03/31 17:43:59 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/transitionprobabilityinstdata.h
   @brief instdata class for HG transition probability computation
 */

#ifndef _HG_TRANSITIONPROBABILITYINSTDATA_H_
#define _HG_TRANSITIONPROBABILITYINSTDATA_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class TransitionProbabilityParams;
  class TransitionProbabilityMeshManager;
}

namespace hg
{


class TransitionProbabilityInstData : public BackwardInstData
{

public:

  TransitionProbabilityInstData
  (pricing::TransitionProbabilityParams& params,
   Model& model,
   pricing::TransitionProbabilityMeshManager& meshes);

  virtual void Init();

  virtual void UpdateBeforeStep();

  virtual void SetInitialValue();

  size_t GetIndexSpot() const { return m_nIdxSpot; }

protected:

  /// The params related to transition probability
  pricing::TransitionProbabilityParams& m_tpParams;
  
  /// The transition probability mesh manager
  pricing::TransitionProbabilityMeshManager& m_tpMeshes;

private:

  size_t m_nIdxReferenceRegime;

  size_t m_nIdxSpot;

  NO_COPY_CLASS(TransitionProbabilityInstData);

}; // class TransitionProbabilityInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_TRANSITIONPROBABILITYINSTDATA_H_
