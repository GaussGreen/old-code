/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/edsstepper.h
// Purpose:     EDS stepper class
// Created:     2005/01/26
// RCS-ID:      $Id: edsstepper.h,v 1.2 2006/06/13 15:34:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/edsstepper.h
    @brief EDS stepper class

    @todo This is nearly the same as CDSStepper. We should find a way to share
          the implementation by using the base class BackWardInstData.
 */

#ifndef _IHG_EDSSTEPPER_H_
#define _IHG_EDSSTEPPER_H_

#include "ito33/autoptr.h"

#include "ihg/stepper.h"
#include "ihg/edsinstdata.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalSolver;
}

namespace ihg
{


/**
    EDS stepper class.
 */
class EDSStepper : public Stepper
{
public:
  
  EDSStepper(EDSInstData& instdata, const finance::ComputationalFlags& flags) 
           : Stepper(instdata, flags),
             m_instdata(instdata)
  {
  }

  ~EDSStepper()
  {
  }

  // implement fucntions in base class
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

  EDSInstData& m_instdata;


private:
  
  NO_COPY_CLASS(EDSStepper);

}; // class EDSStepper


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_EDSSTEPPER_H_

