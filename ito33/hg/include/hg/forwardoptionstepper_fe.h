/////////////////////////////////////////////////////////////////////////////
// Name:        hg/forwardoptionstepper_fe.h
// Purpose:     HG stepper class for forward options using finite elements
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionstepper_fe.h,v 1.8 2006/06/13 16:25:24 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/forwardoptionstepper_fe.h
    @brief HG finite element stepper class for forward options
 */

#ifndef _HG_FORWARDOPTIONSTEPPER_FE_H_
#define _HG_FORWARDOPTIONSTEPPER_FE_H_

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparseconstraintsolver.h"
#include "ito33/numeric/trisparselinearsolver.h"

#include "hg/forwardoptioninstdata.h"
#include "hg/stepper_fe_fix.h"

namespace ito33
{

namespace hg
{


/// HG stepper class for forward options
class ForwardOptionStepperFE : public StepperFEFix
{

public:

  ForwardOptionStepperFE(ForwardOptionInstData& instdata,
                         const finance::ComputationalFlags& flags) 
                       : StepperFEFix(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~ForwardOptionStepperFE() { }

 
protected:

  /**
      Constructs the discrete equation coefficients.

      Redefinition of implementation in StepperFEFix class
   */
  void MakeCoefficients();

  /**
      Constructs the discrete sensitivity equation coefficients.

      Redefinition of implementation in StepperFEFix class
   */
  virtual void MakeSensitivityCoefficients( const SensitivityData& data);

  /**
      Builds the constant sparse matrix due to the jump terms.

      Redefinition of implementation in StepperFE class
   */
  void BuildJumpSystem();
  
  /**
      Builds the constant sparse matrix for the intensity sensitivity.

      Redefinition of implementation in StepperFE class
   */
  AutoPtr<numeric::MorseMatrix>
  BuildIntensitySensitivityJumpSystem(size_t nIdxR1, size_t nIdxR2,
                                      double dAmplitude, double dIntensity);

  /**
      Builds the constant sparse matrix for the amplitude sensitivity.

      Redefinition of implementation in StepperFE class
   */
  AutoPtr<numeric::MorseMatrix>
  BuildAmplitudeSensitivityJumpSystem(size_t nIdxR1, size_t nIdxR2,
                                      double dAmplitude, double dIntensity);
  
  void BuildDerivedJumpSystem();

  /// The type info of instdata is needed (e.g. for Greeks)
  ForwardOptionInstData& m_instdata;


private:

  void
  BuildSensitivityJumpSystemStructure
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude,
   std::list<size_t>* ppnColumnLists);

  void
  BuildPartialJumpSystem
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dMultiplier,
   numeric::MorseMatrix& sparseMatrix);

  void
  BuildAmplitudeSensitivityJumpSystem
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity,
   numeric::MorseMatrix& sparseMatrix);

  NO_COPY_CLASS(ForwardOptionStepperFE);

}; // class ForwardOptionStepperFE;


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_FORWARDOPTIONSTEPPER_FE_H_
