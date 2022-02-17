/////////////////////////////////////////////////////////////////////////////
// Name:        hg/herostepper.h
// Purpose:     HG stepper class for HERO
// Created:     2005/09/26
// RCS-ID:      $Id: herostepper.h,v 1.4 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/herostepper.h
    @brief HG stepper class for HERO
 */

#ifndef _HG_HEROSTEPPER_H_
#define _HG_HEROSTEPPER_H_

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparselinearsolver.h"

#include "hg/heroinstdata.h"

#include "hg/stepper_fix.h"

namespace ito33
{

namespace hg
{

typedef StepperFix HeroStepperBase;

/// HG stepper class for HERO
class HeroStepper : public HeroStepperBase
{

public:

  HeroStepper(HeroInstData& instdata, const finance::ComputationalFlags& flags) 
            : HeroStepperBase(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~HeroStepper() { }

  /**
      Makes a single timestep.
      
      Called by engine to make a single (possibly iterative) timestep.
      This function is pure virtual in the base classes.
   */
  void Run();

  /**
      Initialization of the class data.

      Most data is initialized by the base class Init() function, but the
      computational grid must still be defined (and delta, gamma).
   */
  void Init();

protected:

  /// Construct the discrete equation coefficients at a timestep
  virtual void MakeCoefficients();

  /// The type info of instdata is needed (e.g. for Greeks)
  HeroInstData& m_instdata;


private:

  NO_COPY_CLASS(HeroStepper);

}; // class HeroStepper;


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_HEROSTEPPER_H_
