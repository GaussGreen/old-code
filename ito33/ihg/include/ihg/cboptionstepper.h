////////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cboptionstepper.h
// Purpose:     cb option stepper
// Author:      Nabil
// Created:     2005/10/17
// RCS-ID:      $Id: cboptionstepper.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
////////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cboptionstepper.h
    @brief cb option stepper class

    Implementation of the system solving (part of timestepping) class for 
    cb option for the price and greeks.
 */

#ifndef _IHG_CBOPTIONSTEPPER_H_
#define _IHG_CBOPTIONSTEPPER_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/tridiagonalconstraintsolver.h"
#include "ito33/numeric/tridiagonalsolver.h"

#include "ihg/cboptioninstdata.h"
#include "ihg/cbstepper.h"
#include "ihg/stepper.h"

namespace ito33
{

namespace ihg
{

/**
    CB option stepper class.

    Class to make a single timestep for cb option contracts.  Each step
    calculates the price. It must also calculate greeks(by PDE, if possible)
    if requested by the user.
 */
class CBOptionStepper : public Stepper
{

public:

  CBOptionStepper(CBOptionInstData& cboptioninstdata,
                  const finance::ComputationalFlags& flags);

  virtual ~CBOptionStepper() {}
  
  /**
      Initialization the class data.

      This function is required by engine and is pure virtual in the base class.
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

  /** 
      Constructs the cb option vega equation coefficients. 

      @important{This function only compute the m_pdCoeConst because the other
                 coefficients are the same as those of the cb price then,
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

  /// The diagonal of the matrix (for must be saved and retored in some case)
  Array<double> m_pdDiagonalOfMatrix;  

  /// Working vector
  std::vector<double> m_pdTmp1;

  /// The type info of instdata is needed
  CBOptionInstData& m_cboptioninstdata;

  /// cbstepper for the underlying cb
  shared_ptr<CBStepper> m_pCBStepper;


private:

  NO_COPY_CLASS(CBOptionStepper);
};

} // namespace ihg

} // namespace ito33 


#endif // #ifndef _IHG_CBOPTIONSTEPPER_H_
