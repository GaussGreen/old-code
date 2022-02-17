/////////////////////////////////////////////////////////////////////////////
// Name:        hg/stepper_fe_fix.h
// Purpose:     Stepper class for problem with fixed space mesh
// Created:     2005/05/06
// RCS-ID:      $Id: stepper_fe_fix.h,v 1.8 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/stepper_fe_fix.h
    @brief Stepper class for problem with fixed space mesh using finite element.
 */

#ifndef _HG_STEPPER_FE_FIX_H_
#define _HG_STEPPER_FE_FIX_H_

#include "hg/stepper_fe.h"
#include "hg/backwardinstdata.h"

namespace ito33
{

namespace numeric
{
  class TriSparseLinearSolver;
}

namespace hg
{


/// Stepper class for problem with fixed space mesh
class StepperFEFix : public StepperFE
{
public:
  
  StepperFEFix(BackwardInstData& instdata,
               const finance::ComputationalFlags& flags) 
             : StepperFE(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~StepperFEFix()
  {
  }

  /// Initialize class data
  virtual void Init();

  /// Make a timestep
  virtual void Run();


protected:
  
  /// Construct the constant discrete equation coefficients
  virtual void MakeCoefficients();

  /// Construct the sensitivity equation coefficients
  virtual void MakeSensitivityCoefficients( const SensitivityData& data);

  /// Build the main system for price solving
  void BuildMainSystem();

  virtual void BuildDerivedJumpSystem();

  /// Store the matrix and coefficients if dual system is required
  void SetupDualSystemData();
 
  /// Store datas required for adjoint method, including dual system data
  void SetupSensivityByAdjointData();

  /// The solver for linear equations (no constraints)
  AutoPtr<numeric::TriSparseLinearSolver> m_pLinearSolver;

  /// The instdata
  BackwardInstData& m_instdata;

  
private:
  
  NO_COPY_CLASS(StepperFEFix);

}; // class StepperFEFix


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_STEPPER_FE_FIX_H_
