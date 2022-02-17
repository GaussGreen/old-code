/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/forwardoptionstepper.h
// Purpose:     option stepper class for forward PDE
// Author:      Wang
// Created:     2004/03/10
// RCS-ID:      $Id: forwardoptionstepper.h,v 1.7 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/forwardoptionstepper.h
   @brief option stepper class for forward PDE
 */

#ifndef _IHG_FORWARDOPTIONSTEPPER_H_
#define _IHG_FORWARDOPTIONSTEPPER_H_

#include "ito33/array.h"
#include "ito33/autoptr.h"

#include "ito33/mv/mvvd.h"

#include "ihg/forwardoptioninstdata.h"
#include "ihg/stepper.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalSolver;
  class TridiagonalPreconditioner;
  class GmresSolver;
}

namespace ihg
{

class ForwardOptionStepper : public Stepper
{
public:

  ForwardOptionStepper(ForwardOptionInstData& instdata,
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

  /// The space grid
  const double *m_pdX;

  /// Value of the integral term at each space point
  // No longer needed. Store directly in m_pdCoeConst
  //Array<double> m_pdIntegralValues;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;

  AutoPtr<numeric::GmresSolver> m_pGmresSolver;

  AutoPtr<numeric::TridiagonalPreconditioner> m_pPreconditioner;


protected:
  
  ForwardOptionInstData& m_instdata;


private:

  NO_COPY_CLASS(ForwardOptionStepper);

}; // class ForwardOptionStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_FORWARDOPTIONSTEPPER_H_
