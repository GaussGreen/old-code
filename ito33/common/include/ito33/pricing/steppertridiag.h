/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/steppertridiag.h
// Purpose:     base tridiagonal stepper class
// Author:      ZHANG Yunzhi
// Created:     2003/12/29
// RCS-ID:      $Id: steppertridiag.h,v 1.8 2006/06/13 15:28:51 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/steppertridiag.h
    @brief base tridiagonal stepper class

    Implementation of the system solving (part of timestepping) 
    using tridiagonal matrix.
 */

#ifndef _ITO33_PRICING_STEPPERTRIDIAG_H_
#define _ITO33_PRICING_STEPPERTRIDIAG_H_

#include "ito33/common.h"
#include "ito33/array.h"

#include "ito33/numeric/tridiagonalmatrix.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace numeric
{
  class TridiagonalMatrix;
  class TridiagonalSolver;
}

namespace pricing
{

class Constraints;

/**
    Timestepper base class for tridiagonal problems.

    Base class of timestepping code for problems with tridiagonal
    matrices (such as in the IHG framework, or the HG framework with
    iterative solving for the jump terms). The derived classes
    determine what type of step to make (implicit, Crank-Nicolson,
    etc.). This class simply provides the tridiagonal matrix. 
 */
class StepperTriDiag
{
public:
  StepperTriDiag(const finance::ComputationalFlags& flags);

  virtual ~StepperTriDiag() { }

  /**
      Initialization from m_meshes and m_params.

      All Stepper class must have this member function.

      REQUIRE : m_meshes must have been done
   */
  virtual void Init() = 0;

  /// Run the stepper for a single timestep. Pure virtual.
  virtual void Run() = 0;

protected:

  /**
      Memory allocation for all members.

      It is to be called by Init() function of my derived classes.

      @param nNb max number of the unknown that the stpper will handle
   */
  void Alloc(size_t nNb)
  {
    m_nNbX = nNb;

    m_pdRHS = Array<double>(m_nNbX);

    m_tridiagonalMatrix.SetDimension(m_nNbX);
  }

  /// The size of the grid
  size_t m_nNbX;

  /// The (tridiagonal) matrix
  numeric::TridiagonalMatrix m_tridiagonalMatrix;

  /// The second member (right hand side) 
  Array<double> m_pdRHS; 

  /// The type of the solver for non linear system, 0 = Penalty, 1 = Frozen
  int m_iSolverType;

private:

  NO_COPY_CLASS(StepperTriDiag);
};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_STEPPERTRIDIAG_H_
