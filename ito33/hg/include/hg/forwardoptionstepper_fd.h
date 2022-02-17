/////////////////////////////////////////////////////////////////////////////
// Name:        hg/forwardoptionstepper_fd.h
// Purpose:     HG stepper class for forward options using finite difference
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionstepper_fd.h,v 1.3 2006/06/13 16:25:24 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/forwardoptionstepper.h
    @brief HG finite difference stepper class for forward options
 */

#ifndef _HG_FORWARDOPTIONSTEPPER_FD_H_
#define _HG_FORWARDOPTIONSTEPPER_FD_H_

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparseconstraintsolver.h"
#include "ito33/numeric/trisparselinearsolver.h"

#include "hg/forwardoptioninstdata.h"
#include "hg/stepper_fix.h"

namespace ito33
{

namespace hg
{


/// HG stepper class for forward options
class ForwardOptionStepperFD : public StepperFix
{

public:

  ForwardOptionStepperFD(ForwardOptionInstData& instdata,
                         const finance::ComputationalFlags& flags) 
                       : StepperFix(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~ForwardOptionStepperFD() { }

 
protected:

  /**
      Constructs the discrete equation coefficients.

      Redefinition of implementation in StepperFix class
   */
  void MakeCoefficients();

  /**
      Constructs the discrete sensitivity equation coefficients.

      Redefinition of implementation in StepperFix class
   */
  void MakeSensitivityCoefficients( const SensitivityData& data);

  /**
      Builds the constant sparse matrix due to the jump terms.

      Redefinition of implementation in Stepper class
   */
  void BuildJumpSystem();

  /**
      Builds the sparse matrix for the sensitivity of a jump .

      Redefinition of implementation in Stepper class
   */
  AutoPtr<numeric::MorseMatrix>
  BuildJumpSensitivitySystem(size_t nIdxR1, size_t nIdxR2,
                             double dAmplitude, double dIntensity, 
                             bool bIsAmplitude);


  /// The type info of instdata is needed (e.g. for Greeks)
  ForwardOptionInstData& m_instdata;


private:

  NO_COPY_CLASS(ForwardOptionStepperFD);

}; // class ForwardOptionStepperFD;


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_FORWARDOPTIONSTEPPER_FD_H_
