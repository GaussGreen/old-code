////////////////////////////////////////////////////////////////////////////////
// Name:        hg/cboptionstepper.h
// Purpose:     cb option stepper
// Created:     2006/01/19
// RCS-ID:      $Id: cboptionstepper.h,v 1.2 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
////////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cboptionstepper.h
    @brief cb option stepper class

    Implementation of the system solving (part of timestepping) class for 
    cb option for the price and greeks.
 */

#ifndef _HG_CBOPTIONSTEPPER_H_
#define _HG_CBOPTIONSTEPPER_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/trisparselinearsolver.h"
#include "ito33/numeric/trisparseconstraintsolver.h"

#include "hg/cboptioninstdata.h"
#include "hg/cbstepper.h"

namespace ito33
{

namespace hg
{

/**
    Class to make a single timestep for cb option contracts.  Each step 
    calculates the price, and greeks if requested by the user.
 */
class CBOptionStepper : public Stepper
{

public:

  CBOptionStepper(CBOptionInstData& cboptioninstdata,
                  const finance::ComputationalFlags& flags);

  virtual ~CBOptionStepper() { }

  void Init();
  
  void Run();

protected:
  
  /// Constructs the constant discrete equation coefficients
  void MakeCoefficients();

  AutoPtr<numeric::TriSparseLinearSolver> m_pLinearSolver;
  
  AutoPtr<numeric::TriSparseConstraintSolver> m_pIterativeSolver; 

  /// The instdata that is needed
  CBOptionInstData& m_cboptioninstdata;

  /// cbstepper for the underlying cb
  CBStepper m_cbstepper;

private:

  NO_COPY_CLASS(CBOptionStepper);
};

} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_CBOPTIONSTEPPER_H_
