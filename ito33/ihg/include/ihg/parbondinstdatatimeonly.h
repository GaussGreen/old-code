/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parbondinstdatatimeonly.h
// Purpose:     time only parbond instdata class
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondinstdatatimeonly.h,v 1.1 2005/06/08 15:53:35 zhang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_ParBondINSTDATATIMEONLY_H_
#define _IHG_ParBondINSTDATATIMEONLY_H_

#include "ito33/pricing/instdatatimeonly.h"

namespace ito33
{

namespace pricing
{
  class Event;
  class ParBondParams;
  class ParBondMeshManager;
}

namespace ihg
{
  class Model;

/// Inst data class for parbond in the ihg model
class ParBondInstDataTimeOnly : public pricing::InstDataTimeOnly
{
public:

  ParBondInstDataTimeOnly(pricing::ParBondParams& params,
                      Model& model,
                      pricing::ParBondMeshManager& meshes);

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

  pricing::ParBondParams& m_parbondParams;
  
  Model& m_model;
  
  pricing::ParBondMeshManager& m_parbondMeshes;


private:

  NO_COPY_CLASS(ParBondInstDataTimeOnly);

}; // class ParBondInstDataTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ParBondINSTDATATIMEONLY_H_


