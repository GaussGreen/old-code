/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/transitionprobabilitymeshmanager.h
// Purpose:     mesh manager class for transition probability
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobabilitymeshmanager.h,v 1.1 2006/03/31 17:43:48 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/transitionprobabilitymeshmanager.h
    @brief mesh manager class for transition probability
 */

#ifndef _ITO33_PRICING_TRANSITIONPROBABILITYMESHMANAGER_H_
#define _ITO33_PRICING_TRANSITIONPROBABILITYMESHMANAGER_H_

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class TransitionProbabilityParams;

/// Class for managing the space and time meshes for transition probability.
class TransitionProbabilityMeshManager: public BackwardMeshManagerFix
{
public:

  /// Constructor
  TransitionProbabilityMeshManager
  (TransitionProbabilityParams& params, Model& model);

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Setup this option meshmanager
  virtual void SetupMe();
  
  /// Pre-computes recovery values
  virtual void ComputeRecoveryValues();

protected:
  
  /// Construct a non-uniform space mesh 
  virtual void ConstructSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// The params for option
  TransitionProbabilityParams& m_tpParams;

private:

  NO_COPY_CLASS(TransitionProbabilityMeshManager);

}; // class TransitionProbabilityMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_TRANSITIONPROBABILITYMESHMANAGER_H_
