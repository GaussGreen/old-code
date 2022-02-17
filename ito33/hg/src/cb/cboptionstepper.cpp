/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cboptionstepper.cpp
// Purpose:     cb option stepper
// Created:     2006/01/19
// RCS-ID:      $Id: cboptionstepper.cpp,v 1.4 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/src/cb/cboptionstepper.cpp
    @brief implement CBOptionStepper class
 */

#include "ito33/numeric/boundary1d.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/schemetype.h"

#include "ito33/numeric/trisparselinearsolver_gmres.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparseconstraintsolver_penalty.h"
#include "ito33/numeric/trisparseconstraintsolver_frozen.h"

#include "hg/cbinstdata.h"
#include "hg/cbstepper.h"
#include "hg/cboptionstepper.h"

namespace ito33
{

namespace hg
{

CBOptionStepper::CBOptionStepper
(CBOptionInstData& cboptioninstdata, const finance::ComputationalFlags& flags) 
    : Stepper(cboptioninstdata, flags),
      m_cboptioninstdata(cboptioninstdata),
      m_cbstepper(cboptioninstdata.GetCBInstData(), flags)
{
}

void CBOptionStepper::Init()
{
  // Model specific initialization at once since the model is homogeneous
  Stepper::Init();

  // Init for the underlying cb
  m_cbstepper.Init();

  // Some allocations.

  // Use the max size of the space mesh to allocate helper arrays
  m_nNbS = m_cboptioninstdata.GetNbSpotsMax();
  m_nNbX = m_nNbS * m_nNbRegimes;
 
  Alloc(m_nNbX);

  m_pLinearSolver = AutoPtr<numeric::TriSparseLinearSolver>
      ( new numeric::TriSparseLinearSolverFixedPoint(m_nNbX) );

  if ( m_iSolverType )      
    m_pIterativeSolver = AutoPtr<numeric::TriSparseConstraintSolver>
      ( new numeric::TriSparseConstraintSolverFrozen(m_nNbRegimes, m_nNbX) );
  else
    m_pIterativeSolver = AutoPtr<numeric::TriSparseConstraintSolver>
      ( new numeric::TriSparseConstraintSolverPenalty(m_nNbRegimes, m_nNbX) );
}

void CBOptionStepper::Run()
{
  // Run for the underlying cb
  m_cbstepper.Run();
  
  if ( !m_cboptioninstdata.IsCBOptionWindow() )
    return;

  // Update some data
  m_nNbS = m_cboptioninstdata.m_nNbS;
  m_nNbX = m_nNbS * m_nNbRegimes;
  m_pdX = m_cboptioninstdata.m_pdS;
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdS2[nIdx] = m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];

  m_tridiagonalMatrix.SetDimension(m_nNbX);

  CalculateAreaArrays(m_pdX, m_nNbS);

  // Build the sparse matrix with entries corresponding to jump terms
  BuildJumpSystem();

  // Compute the CB option price =============================================

  // Make the PDE coefficients for the price
  MakeCoefficients();

  // Make the alpha/beta coefficients
  BuildAlphaBeta();
  
  // Build the matrix
  BuildMatrix();
  
  // Build the right hand side for the price equation
  BuildRHS(m_cboptioninstdata.m_pdOldPrices.Get(),
           m_cboptioninstdata.m_pdOldOldPrices.Get());
  
  // Compute the cb option constraints for the cb option price
  m_cboptioninstdata.UpdateCBOptionConstraint();

  m_pIterativeSolver->Solve(m_tridiagonalMatrix,
                            m_sparseMatrix,
                            m_pdRHS.Get(),
                            m_cboptioninstdata.GetCBOptionConstraint(),
                            m_cboptioninstdata.m_pdPrices.Get(),
                            m_cboptioninstdata.m_piFrozenFlags.Get());
}

void CBOptionStepper::MakeCoefficients()
{
  // Get the various interest rates. These usually change at each timestep
  double dRate = m_cboptioninstdata.m_dRate;
  double dForeignRate = m_cboptioninstdata.m_dForeignRate; 
  double dDerivativeRate = m_cboptioninstdata.m_dDerivativeRate;

  // Recovery assumption: no recovery at all in case of default

  /*
  // we have only two cases now. so we can use if{}else{}. if, in the future,
  // we get more condition check to do, we'd better use m_pdHROfDerivative
  // array in cbinstdata to avoid "if" check. In this way, we should 
  // initialize the array according to the value of m_bDerivativeHasOwnHR
  // in UpdateBeforeStep().
  if ( m_cbinstdata.m_bDerivativeHasOwnHR )
  {
  }
  else */
  {
    double dSpeedCorrection = m_cboptioninstdata.GetSpeedCorrection();

    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
      m_pdCoe1st[nIdxR] = m_pdCoe1st0[nIdxR] + dRate - dForeignRate  
                        - dSpeedCorrection;
      m_pdCoeZero[nIdxR] = m_pdCoeZero0[nIdxR] + dDerivativeRate;
    }
    
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeConst[nIdx] = 0.;
  }
}

} // namespace hg

} // namespace ito33
