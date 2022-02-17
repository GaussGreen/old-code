/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parbondstepper.h
// Purpose:     parbond stepper class
// Author:      Nabil
// Created:     2005/05/20
// RCS-ID:      $Id: parbondstepper.h,v 1.2 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parbondstepper.h
    @brief parbond stepper class

    Implementation of the system solving (part of timestepping) class for
    parbond.
 */

#ifndef _IHG_ParBondSTEPPER_H_
#define _IHG_ParBondSTEPPER_H_

#include "ito33/autoptr.h"

#include "ihg/stepper.h"
#include "ihg/parbondinstdata.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalSolver;
}

namespace ihg
{


/**
    Parbond stepper class.

    Class to make a single timestep for regular parbond contracts. Each step
    must also calculate Greek data, if requested by the user.
 */
class ParBondStepper : public Stepper
{
public:
  
  ParBondStepper(ParBondInstData& instdata,
                 const finance::ComputationalFlags& flags) 
               : Stepper(instdata, flags),
                 m_instdata(instdata)
  {
  }

  ~ParBondStepper()
  {
  }

  /**
      Initialization from m_meshes and m_params.

      All Stepper class must have this member function.

      REQUIRE : m_meshes must have been done
      
      @todo add an assert instrument to check this condition
   */
  void Init();

  /**
      Makes a single timestep.
      
      Called by engine to make a single (possibly iterative) timestep.
      This function is pure virtual in the base classes.
   */
  void Run();


protected:
  
  /// Construct the constant discrete equation coefficients
  void MakeCoefficients();

  /// Construct the vega equation coefficients
  void MakeVegaCoefficients();

  /// the space mesh
  const double* m_pdX;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;

  ParBondInstData& m_instdata;


private:
  
  NO_COPY_CLASS(ParBondStepper);

}; // class ParBondStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_ParBondSTEPPER_H_
