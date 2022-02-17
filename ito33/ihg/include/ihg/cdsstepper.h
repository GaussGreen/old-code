/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cdsstepper.h
// Purpose:     cds stepper class
// Author:      Nabil
// Created:     2003/10/31
// RCS-ID:      $Id: cdsstepper.h,v 1.8 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cdsstepper.h
    @brief cds stepper class

    Implementation of the system solving (part of timestepping) class for cds.
 */

#ifndef _IHG_CDSSTEPPER_H_
#define _IHG_CDSSTEPPER_H_

#include "ito33/autoptr.h"

#include "ihg/stepper.h"
#include "ihg/cdsinstdata.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalSolver;
}

namespace ihg
{


/**
    cds stepper class.

    Class to make a single timestep for regular cds contracts. Each step
    must also calculate Greek(by PDE) data, if requested by the user.
 */
class CDSStepper : public Stepper
{
public:
  
  CDSStepper(CDSInstData& instdata, const finance::ComputationalFlags& flags) 
           : Stepper(instdata, flags),
             m_instdata(instdata)
  {
  }

  ~CDSStepper()
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

  // the space mesh
  const double* m_pdX;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;

  CDSInstData& m_instdata;


private:
  
  NO_COPY_CLASS(CDSStepper);

}; // class CDSStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_CDSSTEPPER_H_
