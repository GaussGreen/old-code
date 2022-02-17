/////////////////////////////////////////////////////////////////////////////
// Name:        hg/varianceswapstepper.h
// Purpose:     HG stepper class for variance swaps
// Created:     2006/03/05
// RCS-ID:      $Id: varianceswapstepper.h,v 1.4 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/varianceswapstepper.h
    @brief HG stepper class for variance swaps
 */

#ifndef _HG_VARIANCESWAPSTEPPER_H_
#define _HG_VARIANCESWAPSTEPPER_H_

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparseconstraintsolver.h"
#include "ito33/numeric/trisparselinearsolver.h"

#include "hg/stepper_fix.h"
#include "hg/varianceswapinstdata.h"

namespace ito33
{

namespace hg
{


typedef StepperFix VarianceSwapStepperBase;

/// HG stepper class for variance swaps
class VarianceSwapStepper : public VarianceSwapStepperBase
{

public:

  VarianceSwapStepper(VarianceSwapInstData& instdata, 
                      const finance::ComputationalFlags& flags) 
      : VarianceSwapStepperBase(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~VarianceSwapStepper() { }

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
  
  /// The type info of instdata is needed (e.g. for Greeks)
  VarianceSwapInstData& m_instdata;


private:

  NO_COPY_CLASS(VarianceSwapStepper);

}; // class VarianceSwapStepper;


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_VARIANCESWAPSTEPPER_H_
