/////////////////////////////////////////////////////////////////////////////
// Name:        hg/stepper.h
// Purpose:     HG Base stepper class (finite differences)
// Created:     2005/01/13
// RCS-ID:      $Id: stepper.h,v 1.8 2006/06/13 15:51:43 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/stepper.h
   @brief HG Base stepper class (finite differences)

   @remark Crank-Nicolson won't be supproted.
 */

#ifndef _HG_STEPPER_H_
#define _HG_STEPPER_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"

#include "ito33/pricing/steppertridiag.h"

#include "ito33/numeric/morsematrix.h"

#include "hg/instdata.h"

namespace ito33
{

namespace hg
{

  class Model;

class Stepper : public pricing::StepperTriDiag
{
public:

  Stepper(InstData& instdata, const finance::ComputationalFlags& flags) 
        : pricing::StepperTriDiag(flags), m_instdata(instdata)
  {
  }

  virtual ~Stepper() { }

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

  /// Construct the discrete equation coefficients at a timestep
  virtual void MakeCoefficients() = 0;

  /// Precalculate grid spacings for alpha/beta construction
  void CalculateAreaArrays(const double *pdX, size_t nNbX);

  /// Modify the linear boundary condition coefficient for space grid pdX
  void ModifyLinearBoundaryConditionCoef(const double *pdX, size_t nNbX);

  /// Build the the right hand side using previously calculated coefficients
  void BuildRHS(const double* pdPrice1, const double* pdPriceM);
 
  /// Build the discrete system using previously calculated coefficients
  void BuildAlphaBeta();
  
  /// Build the matrix using previously calculated coefficients
  void BuildMatrix();

  /// Build the sparse matrix due to jumps
  virtual void BuildJumpSystem();

  /// Compute some helper arrays for sensitivity
  void ComputeHelperArrays();

  /// Build the sparse matrix for the sensitivity of a jump 
  virtual AutoPtr<numeric::MorseMatrix>
  BuildJumpSensitivitySystem(size_t nIdxR1, size_t nIdxR2,
                             double dAmplitude, double dIntensity, 
                             bool bIsAmplitude);

  /**
      Memory allocation for all members.
 
      It is to be called by Init() function of my derived classes
 
      @param nNb max number of the unknown that the stpper will handle
   */
  void Alloc(size_t nNb);

  numeric::MorseMatrix m_sparseMatrix;

  const double* m_pdX;

  size_t m_nNbRegimes;

  size_t m_nNbS;

  Array<int> m_piFD;

  Array<double> m_pdS2;

  std::vector<double> m_pdVols;

  std::vector<double> m_pdDefaultIntensities;

  const Model* m_pModel;

  /// Linear boundary condition coefficient on the first space point.
  double m_dLBC_CoefLeft;

  /// Linear boundary condition coefficient for the last unkown point.
  double m_dLBC_CoefRight;

  /// The grid spacings
  Array<double> m_pdDeltaX;

  /// The grid finite volumes/areas
  Array<double> m_pdAreaX;

  /// The grid spacings
  Array<double> m_pdInverseDeltaX;

  /// The grid finite volumes/areas
  Array<double> m_pdInverseAreaX;

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

  /// Discretization values (alpha), constant for the price and Greek PDEs
  Array<double> m_pdAlpha;

  /// Discretization values (beta), constant for the price and Greek PDEs
  Array<double> m_pdBeta;


protected:

  InstData& m_instdata;

  NO_COPY_CLASS(Stepper);
};

} // namespace hg

} // namespace ito33 


#endif // #ifndef _HG_STEPPER_H_
