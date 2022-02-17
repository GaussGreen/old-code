/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/onetouchstepper.h
// Purpose:     OneTouch stepper class
// Created:     2006/08/10
// RCS-ID:      $Id: onetouchstepper.h,v 1.1 2006/08/10 23:10:42 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/onetouchstepper.h
    @brief OneTouch stepper class

    @todo This is nearly the same as EDSStepper. We should find a way to share
          the implementation by using the base class BackWardInstData. See HG.
 */

#ifndef _IHG_ONETOUCHSTEPPER_H_
#define _IHG_ONETOUCHSTEPPER_H_

#include "ito33/autoptr.h"

#include "ihg/stepper.h"
#include "ihg/onetouchinstdata.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalSolver;
}

namespace ihg
{


/// OneTouch stepper class.
class OneTouchStepper : public Stepper
{
public:
  
  OneTouchStepper(OneTouchInstData& instdata,
                  const finance::ComputationalFlags& flags) 
                : Stepper(instdata, flags), m_instdata(instdata)
  {
  }

  ~OneTouchStepper()
  {
  }

  // implement functions in base class
  void Init();

  void Run();


protected:
  
  /// Construct the constant discrete equation coefficients
  void MakeCoefficients();

  /// Construct the vega equation coefficients
  void MakeVegaCoefficients();

  // the space mesh
  const double *m_pdX;

  AutoPtr<numeric::TridiagonalSolver> m_pLinearSolver;

  OneTouchInstData& m_instdata;


private:
  
  NO_COPY_CLASS(OneTouchStepper);

}; // class EDSStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_ONETOUCHSTEPPER_H_

