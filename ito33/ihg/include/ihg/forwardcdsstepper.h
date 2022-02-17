/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/forwardcdsstepper.h
// Purpose:     CDS stepper class for forward PDE
// Author:      David
// Created:     2004/03/29
// RCS-ID:      $Id: forwardcdsstepper.h,v 1.6 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/forwardcdsstepper.h
   @brief CDS stepper class for forward PDE
 */

#ifndef _IHG_FORWARDCDSSTEPPER_H_
#define _IHG_FORWARDCDSSTEPPER_H_

#include "ito33/array.h"
#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalpreconditioner.h"
#include "ito33/numeric/gmressolver.h"

#include "ihg/forwardcdsinstdata.h"
#include "ihg/stepper.h"

namespace ito33
{

namespace ihg
{

class ForwardCDSStepper : public Stepper
{
public:

  ForwardCDSStepper(ForwardCDSInstData& instdata, 
                    const finance::ComputationalFlags& flags) 
                  : Stepper(instdata, flags), m_instdata(instdata)
  {
  }

  // Default dtor is ok

  void Init();

  void Run();

  CVecteurDouble operator*(CVecteurDouble &X) const;

protected:

  /// Construct the constant discrete equation coefficients
  void MakeCoefficients();
  
  /// Compute the integral term in the forward PDE
  void ComputeIntegralTerm(double* pdPrices);

  /// Calculate the terms required for the CDS pricing
  void CalculateCDSTerms();

  /// The space grid
  const double *m_pdX;

  /// Value of the integral term at each sace point
  //Array<double> m_pdIntegralValues;

  /// The solver
  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;

  AutoPtr<numeric::GmresSolver> m_pGmresSolver;

  AutoPtr<numeric::TridiagonalPreconditioner> m_pPreconditioner;

  /// The class storing all the timestepping data  
  ForwardCDSInstData& m_instdata;


private:

  NO_COPY_CLASS(ForwardCDSStepper);

}; // class ForwardCDSStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_FORWARDCDSSTEPPER_H_
