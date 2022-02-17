/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/varinceswapstepper.h
// Purpose:     variance swap stepper class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapstepper.h,v 1.3 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/varianceswapstepper.h
    @brief variance swap stepper class

    @todo This is nearly the same as XXXStepper. We should find a way to share
          the implementation by using the base class BackWardInstData.
 */

#ifndef _IHG_VARIANCESWAPSTEPPER_H_
#define _IHG_VARIANCESWAPSTEPPER_H_

#include "ito33/autoptr.h"

#include "ihg/stepper.h"
#include "ihg/varianceswapinstdata.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalSolver;
}

namespace ihg
{


/// Variance swap stepper class
class VarianceSwapStepper : public Stepper
{
public:
  
  VarianceSwapStepper(VarianceSwapInstData& instdata, 
                      const finance::ComputationalFlags& flags) 
                    : Stepper(instdata, flags), m_instdata(instdata)
  {
  }

  ~VarianceSwapStepper()
  {
  }

  // implement fucntions in base class
  void Init();

  void Run();


protected:
  
  /// Construct the constant discrete equation coefficients
  void MakeCoefficients();

  /// the space mesh
  const double* m_pdX;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;

  VarianceSwapInstData& m_instdata;


private:
  
  NO_COPY_CLASS(VarianceSwapStepper);

}; // class VarianceSwapStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_VARIANCESWAPSTEPPER_H_
