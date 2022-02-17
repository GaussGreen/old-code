////////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cbstepper.h
// Purpose:     cb stepper including greek, straight bond and fugit computations
// Author:      Nabil
// Created:     2004/04/09
// RCS-ID:      $Id: cbstepper.h,v 1.19 2006/07/20 15:14:00 nabil Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
////////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cbstepper.h
    @brief cb stepper class for price, greeks, straight bond and fugit

    Implementation of the system solving (part of timestepping) class for 
    cb for the price, greeks, and fugit.
 */

#ifndef _IHG_CBSTEPPER_H_
#define _IHG_CBSTEPPER_H_

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalfrozensolver.h"
#include "ito33/numeric/tridiagonalsolver.h"

#include "ihg/cbinstdata.h"
#include "ihg/stepper.h"

namespace ito33
{

namespace ihg
{

/**
  @brief cb stepper class

  Class to make a single timestep for cb contracts.  Each step calculate
  price and straight bond. It must also calculate greeks and fugit 
  if requested by the user.
*/
class CBStepper : public Stepper
{

public:

  CBStepper(CBInstData& cbinstdata, const finance::ComputationalFlags& flags) 
          : Stepper(cbinstdata, flags),
            m_cbinstdata(cbinstdata)
  {
  }

  virtual ~CBStepper() { }

  /**
      Makes a single timestep.
      
      Called by engine to make a single (possibly iterative) timestep.
      This function is pure virtual in the base classes.
   */
  void Run();

  /**
      Initialization the class data.

      This function is required by engine and is pure virtual in the base class.
   */
  void Init();

protected:
  
  /// Construct the constant discrete equation coefficients
  void MakeCoefficients();
  
  /// Construct the new share equation coefficients
  void MakeNewShareCoefficients();
  
  /// Make a single time step for new share pricing
  void RunNewShareStep();

  /** 
      Constructs the vega equation coefficients.

      @important{This function only compute the m_pdCoeConst then,
                 to compute m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero,
                 MakeCoefficients() function must be called before}
   */
  void MakeVegaCoefficients();

  /** 
      Constructs the fugit equation coefficients.

      @important{This function only computes the m_pdCoeZero and m_pdCoeConst 
                 then, to compute m_pdCoe2nd and m_pdCoe1st, MakeCoefficients() 
                 function must be called before}
    */
  void MakeFugitCoefficients();

  /** 
      Changes the diagonal of the matrix after saving it in the member 
      m_pdDiagonalOfMatrix.
    
      This is due to the fact that, for some PDEs as the one of the fugit, 
      not only m_pdCoeConst changes compared with the price case but also 
      m_pdCoeZero which intervene in the diagonal of the matrix.
    
      @important{If the PDE of the computation of XXX have a m_pdCoeZero different
                 from the one of the price, the function MakeXXXCoefficients() 
                 must be called before and at the end of the computation of the 
                 XXX, the function RestoreDiagonalOfMatrix() must be called.}
    
   */
  void ChangeDiagonalOfMatrix();

  /** 
      Restores the diagonal of the matrix.
    
      @important{This function must be called at the end of the computation
                 where the function ChangeDiagonalOfMatrix() has been used.}

   */
  void RestoreDiagonalOfMatrix();
  
  
  /** 
      The space mesh.

      We are using S grid here.
   */
  const double* m_pdX;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;
  
  AutoPtr<numeric::TridiagonalConstraintSolver> m_pIterativeSolver;  

  /// The old zero fugit PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdFugitCoeZeroOld;

  /// The old constant fugit PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdFugitCoeConstOld; 

  /// The old zero new share PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdNewShareCoeZeroOld;

  /// The old constant new share PDE coefficient (for Crank-Nicolson)
  Array<double> m_pdNewShareCoeConstOld;

  /// The diagonal of the matrix (for must be saved and retored in some case)
  Array<double> m_pdDiagonalOfMatrix;

  /// The type info of instdata is needed (for greeks, straight bond and fugit)
  CBInstData& m_cbinstdata;

private:

  NO_COPY_CLASS(CBStepper);

  friend class CBOptionStepper;

};

} // namespace ihg

} // namespace ito33 


#endif // #ifndef _IHG_CBSTEPPER_H_
