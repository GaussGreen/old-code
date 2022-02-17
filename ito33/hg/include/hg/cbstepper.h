////////////////////////////////////////////////////////////////////////////////
// Name:        hg/cbstepper.h
// Purpose:     cb stepper including greek, and fugit computations
// Created:     2005/04/11
// RCS-ID:      $Id: cbstepper.h,v 1.7 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
////////////////////////////////////////////////////////////////////////////////

/**
   @file hg/cbstepper.h
   @brief cb stepper class for price, greeks and fugit

   Implementation of the system solving (part of timestepping) class for 
   cb for the price, greeks, and fugit.
 */

#ifndef _HG_CBSTEPPER_H_
#define _HG_CBSTEPPER_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/trisparselinearsolver.h"
#include "ito33/numeric/trisparseconstraintsolver.h"

#include "hg/cbinstdata.h"
#include "hg/stepper.h"

namespace ito33
{

namespace hg
{

/**
   cb stepper class

   Class to make a single timestep for cb contracts.  Each step calculate
   price and straight bond. It must also calculate greeks and fugit 
   if requested by the user.
 */
class CBStepper : public Stepper
{

public:

  CBStepper(CBInstData& cbinstdata,
            const finance::ComputationalFlags& flags) 
          : Stepper(cbinstdata, flags), 
            m_cbinstdata(cbinstdata)
  {
  }

  virtual ~CBStepper() { }
  
  /**
      Initialization the class data.

      This function is required by engine and is pure virtual in the base classes.
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
  
  /// Construct the new share equation coefficients
  void MakeNewShareCoefficients();

  /// Make a single time step for new share pricing
  void RunNewShareStep();

  /** 
     Construct the fugit equation coefficients.

     @important{This function only compute the m_pdCoeZero and m_pdCoeConst then,
                to compute m_pdCoe2nd and m_pdCoe1st, MakeCoefficients() 
                function must be called before}
   */
  void MakeFugitCoefficients();

  /** 
     Changes the matrix to compute the fugit.
    
     This is due to the fact that, for the fugit, not only m_pdCoeConst 
     changes compared with the price case but also m_pdCoeZero which 
     intervene in the diagonal of the matrix.
    
     @important{This function must be called for the fugit computation.
                The function MakeFugitCoefficients() must be called before
                and at the end of the computation of the fugit, the function
                RestoreMatrixAfterFugit() must be called.}
    
   */
  void ChangeMatrixForFugit();

  /** 
     Restore the matrix after the computation of the fugit.
    
     @important{This function must be called at the end of the computation of
               the fugit to restore the initial matrix.}
   */
  void RestoreMatrixAfterFugit();

  AutoPtr<numeric::TriSparseLinearSolver> m_pLinearSolver;
  
  AutoPtr<numeric::TriSparseConstraintSolver> m_pIterativeSolver;


protected:

  /// The type info of instdata is needed (for greeks, straight bond and fugit)
  CBInstData& m_cbinstdata;

private:

  NO_COPY_CLASS(CBStepper);
};


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_CBSTEPPER_H_
