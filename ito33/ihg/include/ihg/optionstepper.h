/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/optionstepper.h
// Purpose:     option stepper including Greek computations
// Author:      David, ZHANG Yunzhi
// Created:     2003/12/17
// RCS-ID:      $Id: optionstepper.h,v 1.28 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ihg/optionstepper.h
  @brief option stepper class for price and Greeks

  Implementation of the system solving (part of timestepping) class for 
  options for the price and the Greeks.
*/

#ifndef _IHG_OPTIONSTEPPER_H_
#define _IHG_OPTIONSTEPPER_H_


#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"

#include "ihg/optioninstdata.h"
#include "ihg/stepper.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/tridiagonalpenaltysolver.h"
#include "ito33/numeric/tridiagonalfrozensolver.h"



namespace ito33
{

namespace ihg
{

/**
    Option stepper class.

    Class to make a single timestep for regular option contracts.  Each step
    must also calculate Greek data (by PDE), if requested by the user.
 */
class OptionStepper : public Stepper
{

public:

  OptionStepper(OptionInstData& instdata, 
                const finance::ComputationalFlags& flags) 
              : Stepper(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~OptionStepper() { }

  /**
      Make a single timestep.
      
      Called by engine to make a single (possibly iterative) timestep.
      This function is pure virtual in the base classes.
   */
  void Run();

  /**
      Initialization the class data.

      Most data is initialized by the base class Init() function, but the
      computational grid must still be defined (and delta, gamma).
   */
  void Init()
  {  
    m_nNbX = m_instdata.m_nNbS;

    Alloc(m_nNbX);

    CalculateAreaArrays(m_instdata.m_pdS, m_nNbX);

    m_pdX = m_instdata.m_pdS;

    if(m_instdata.m_pConstraints == 0)
      m_pLinearSolver = AutoPtr<numeric::TridiagonalSolver>
                        ( new numeric::TridiagonalSolver(m_nNbX) );
    else   
    {
      if ( m_iSolverType )      
        m_pIterativeSolver = AutoPtr<numeric::TridiagonalConstraintSolver>
                             ( new numeric::TridiagonalFrozenSolver(m_nNbX) );
      else
        m_pIterativeSolver = AutoPtr<numeric::TridiagonalConstraintSolver>
                             ( new numeric::TridiagonalPenaltySolver(m_nNbX) );
    }
  }

protected:

  /// Construct the constant discrete equation coefficients
  void MakeCoefficients();

  /// Construct the vega equation coefficients
  void MakeVegaCoefficients();
  
  /** 
      The space mesh.

      We are using S grid here.
   */
  const double* m_pdX;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;
  
  AutoPtr<numeric::TridiagonalConstraintSolver> m_pIterativeSolver;

  /// The type info of instdata is needed (e.g. for Greeks)
  OptionInstData& m_instdata;


private:

  NO_COPY_CLASS(OptionStepper);

}; // class OptionStepper;


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_OPTIONSTEPPER_H_

