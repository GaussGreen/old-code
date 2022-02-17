/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/stepper.h
// Purpose:     option stepper class (finite differences)
// Author:      David
// Created:     2003/12/11
// RCS-ID:      $Id: stepper.h,v 1.5 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/stepper.h
    @brief IHG Base stepper class

    Implementation of the system building/solving common to all IHG pricing.
 */

#ifndef _IHG_STEPPER_H_
#define _IHG_STEPPER_H_

#include "ito33/pricing/steppertridiag.h"

#include "ihg/instdata.h"

namespace ito33
{

namespace ihg
{


class Stepper : public pricing::StepperTriDiag
{
public:

  Stepper(InstData& instdata, const finance::ComputationalFlags& flags);

  virtual ~Stepper() { }

  /**
      Initialization from m_meshes and m_params.

      All Stepper class must have this member function.
 
      REQUIRE : m_meshes must have been done.
      
      @todo add an assert instrument to check this condition
   */
  virtual void Init() = 0;

  /**
      Runs the system solver.

      REQUIRE : the instanous parameter part in m_instdata must have
                been updated.
      PROMISE : the result on the current time is calculated and m_instdata
                is thus modified.
   */
  virtual void Run() = 0;


protected:

  /// Precalculate grid spacings for alpha/beta construction
  void CalculateAreaArrays(const double *pdX, size_t nNbX);

  /// Modify the linear boundary condition coefficient for space grid pdX
  void ModifyLinearBoundaryConditionCoef(const double *pdX, size_t nNbX);

  /// Construct the discrete equation coefficients at a timestep
  virtual void MakeCoefficients() = 0;

  /// Build the the right hand side using previously calculated coefficients
  void BuildRHS(const double* pdPrice1, const double* pdPriceM);

  /// Build the matrix using previously calculated coefficients
  void BuildMatrix();

  /// Build the discrete system using previously calculated coefficients
  void BuildAlphaBeta();

  /**
      Memory allocation for all members.

      It is to be called by Init() function of my derived classes

      @param nNb max number of the unknown that the stpper will handle
   */
  void Alloc(size_t nNb);

  /// Linear boundary condition coefficient on the first space point.
  double m_dLBC_CoefLeft;

  /// Linear boundary condition coefficient for the last unkown point.
  double m_dLBC_CoefRight;

  /// The grid spacings
  Array<double> m_pdDeltaX;

  /// The grid finite volumes/areas
  Array<double> m_pdAreaX;

  /// The grid spacings
  Array<double> m_pdInverseDeltaX;

  /// The grid finite volumes/areas
  Array<double> m_pdInverseAreaX;

  /// Diffusion coefficient (Uss)
  Array<double> m_pdCoe2nd;

  /// Convection coefficient (Us)
  Array<double> m_pdCoe1st;

  /// Decay coefficient (U)
  Array<double> m_pdCoeZero;

  /// Constant coefficient (added to RHS)
  Array<double> m_pdCoeConst;

  /// Discretization values (alpha), constant for the price and Greek PDEs
  Array<double> m_pdAlpha;

  /// Discretization values (beta), constant for the price and Greek PDEs
  Array<double> m_pdBeta;

  /// The old alpha values (for Crank-Nicolson)
  Array<double> m_pdAlphaOld;

  /// The old alpha values (for Crank-Nicolson)
  Array<double> m_pdBetaOld;

  /// The old decay PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdCoeZeroOld;

  /// The old constant PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdCoeConstOld;

  /// The old constant vega PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdVegaCoeConstOld;


protected:

  InstData& m_instdata;

private:
  
  typedef pricing::StepperTriDiag BaseClass;

  NO_COPY_CLASS(Stepper);
};

} // namespace ihg

} // namespace ito33 


#endif // #ifndef _IHG_STEPPER_H_
