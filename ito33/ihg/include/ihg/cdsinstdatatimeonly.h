/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cdsinstdatatimeonly.h
// Purpose:     time only cds instdata class
// Author:      Wang
// Created:     2004/03/18
// RCS-ID:      $Id: cdsinstdatatimeonly.h,v 1.7 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_CDSINSTDATATIMEONLY_H_
#define _IHG_CDSINSTDATATIMEONLY_H_

#include "ito33/pricing/instdatatimeonly.h"

namespace ito33
{

namespace pricing
{
  class Event;
  class CDSParams;
  class CDSMeshManager;
}

namespace ihg
{
  class Model;

/// Inst data class for cds in the ihg model
class CDSInstDataTimeOnly : public pricing::InstDataTimeOnly
{
public:

  CDSInstDataTimeOnly(pricing::CDSParams& params,
                      Model& model,
                      pricing::CDSMeshManager& meshes);

  // Default dtor is ok
  
  void Init();

  void UpdateBeforeStep();

  // DoEvents is the same as in the base class InstDataTimeOnly

  void SetInitialValue();

  /** 
     Normally it's implemented in Stepper, but it's just more convient
     to do it here
   */
  void Run();

  void TurnAllFlagsOff() {  }

  // the hazard rate at this time
  double m_dHazardRate;

  /// the price at this time
  double m_dPrice;

  /// the price at the old time step
  double m_dOldPrice;

  /// the price used for BDF time stepping
  double m_dOldOldPrice;

  /// Current default value
  double m_dRecoveryValue;

  /// The decay coefficient
  double m_dCoeZero;

  /// The source term
  double m_dCoeConst;

  /// The old decay coefficient, used by Crank-Nicolson
  double m_dOldCoeZero;

  /// The old source term
  double m_dOldCoeConst;


protected:

  void ApplyEvent(const pricing::Event* pEvent);

  pricing::CDSParams& m_cdsParams;
  
  Model& m_model;
  
  pricing::CDSMeshManager& m_cdsMeshes;


private:

  NO_COPY_CLASS(CDSInstDataTimeOnly);

}; // class CDSInstDataTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CDSINSTDATATIMEONLY_H_


