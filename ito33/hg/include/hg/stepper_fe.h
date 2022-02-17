/////////////////////////////////////////////////////////////////////////////
// Name:        hg/stepper_fe.h
// Purpose:     HG Base stepper class (finite element)
// Created:     2005/01/18
// RCS-ID:      $Id: stepper_fe.h,v 1.12 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/stepper_fe.h
    @brief HG Base stepper class (finite element)
 */

#ifndef _HG_STEPPER_FE_H_
#define _HG_STEPPER_FE_H_

#include "ito33/vector.h"

#include "ito33/pricing/steppertridiag.h"

#include "ito33/numeric/morsematrix.h"

#include "hg/instdata.h"

namespace ito33
{

namespace hg
{

  class Model;

 
/// HG Base stepper class using finite element space discretization. 
class StepperFE : public pricing::StepperTriDiag
{
public:

  StepperFE(InstData& instdata, const finance::ComputationalFlags& flags)
    : pricing::StepperTriDiag(flags), m_instdata(instdata)
  {
  }

  virtual ~StepperFE() { }

  /**
      Initialization from params and meshes.

      Required by Engine.
   */
  virtual void Init();

  /**
      Runs the stepper to do one step solving on time.
   */
  virtual void Run() = 0;


protected:

  /// Precalculate grid spacings for alpha/beta construction
  void BuildMassMatrix();

  /// Construct the discrete equation coefficients at a timestep
  virtual void MakeCoefficients() = 0;

  /// Build the the right hand side using previously calculated coefficients
  void BuildRHS(const double* pdPrice1, const double* pdPriceM);
  
  /// Build the matrix using previously calculated coefficients
  void BuildMatrix();

  /// Build the sparse matrix due to jumps
  virtual void BuildJumpSystem();

  /// Compute some helper arrays for sensitivity
  void ComputeHelperArrays();
 
  void
  BuildSensitivityJumpSystemStructure
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude,
   std::list<size_t>* ppnColumnLists); 

  virtual AutoPtr<numeric::MorseMatrix>
  BuildIntensitySensitivityJumpSystem
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity);

  void
  BuildPartialSensitivityJumpSystem
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dMultiplier,
   numeric::MorseMatrix& sparseMatrix);

  void
  BuildAmplitudeSensitivityJumpSystem
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity,
   numeric::MorseMatrix& sparseMatrix);

  virtual AutoPtr<numeric::MorseMatrix>
  BuildAmplitudeSensitivityJumpSystem
  (size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity);

  numeric::TridiagonalMatrix m_massMatrix;
  
  numeric::MorseMatrix m_sparseMatrix;

  const double* m_pdX;

  size_t m_nNbRegimes;

  size_t m_nNbS;

  std::vector<double> m_pdVols;

  std::vector<double> m_pdDefaultIntensities;

  const Model* m_pModel;

  /**
      Memory allocation for all members.
 
      It is to be called by Init() function of my derived classes
 
      @param nNb max number of the unknown that the stpper will handle
   */
  void Alloc(size_t nNb);

  /// Linear boundary condition coefficient on the first space point.
  double m_dLBC_CoefLeft;

  /// Linear boundary condition coefficient for the last unkown point.
  double m_dLBC_CoefRight;

  /// The grid spacings
  Array<double> m_pdDeltaX;

  /// The grid spacings
  Array<double> m_pdInverseDeltaX;

  /// Constant part of diffusion coefficient (Uss)
  Array<double> m_pdCoe2nd0;

  /// Constant part of convection coefficient (Us)
  Array<double> m_pdCoe1st0;

  /// Constant part of decay coefficient (U)
  Array<double> m_pdCoeZero0;

  /// Diffusion coefficient (Uss)
  Array<double> m_pdCoe2nd;

  /// Convection coefficient (Us)
  Array<double> m_pdCoe1st;

  /// Decay coefficient (U)
  Array<double> m_pdCoeZero;

  /// Constant coefficient (added to RHS)
  Array<double> m_pdCoeConst;

  /// Boolean indicates if we are building sensitivity right hand side
  bool m_bSensitivityRHS;

  /// Constant coefficient (added to RHS)
  Array<double> m_pdCoeConstSensitivity;


protected:

  InstData& m_instdata;

  NO_COPY_CLASS(StepperFE);
};


} // namespace hg

} // namespace ito33 


#endif // #ifndef _HG_STEPPER_FE_H_
