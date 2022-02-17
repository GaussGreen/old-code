/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cbstepper.cpp
// Purpose:     cb stepper including greek, and fugit computations
// Created:     2005/04/11
// RCS-ID:      $Id: cbstepper.cpp,v 1.15 2006/06/13 15:47:09 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "hg/cbstepper.h"

#include "ito33/numeric/trisparselinearsolver_gmres.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparseconstraintsolver_penalty.h"
#include "ito33/numeric/trisparseconstraintsolver_frozen.h"


namespace ito33
{

namespace hg
{

void CBStepper::Init()
{
  // Model specific initialization at once since the model is homogeneous
  Stepper::Init();

  // Use the max size of the space mesh to allocate helper arrays
  m_nNbS = m_cbinstdata.GetNbSpotsMax();
  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  // Solver initialization, mainly memory allocation at once
  m_pLinearSolver = AutoPtr<numeric::TriSparseLinearSolver>
      ( new numeric::TriSparseLinearSolverFixedPoint(m_nNbX) );
  
  if ( m_iSolverType )      
    m_pIterativeSolver = AutoPtr<numeric::TriSparseConstraintSolver>
      ( new numeric::TriSparseConstraintSolverFrozen(m_nNbRegimes, m_nNbX) );
  else
    m_pIterativeSolver = AutoPtr<numeric::TriSparseConstraintSolver>
      ( new numeric::TriSparseConstraintSolverPenalty(m_nNbRegimes, m_nNbX) );


}

void CBStepper::Run()
{
  // Update some data
  m_nNbS = m_cbinstdata.m_nNbS;
  m_nNbX = m_cbinstdata.m_nNbS * m_nNbRegimes;
  m_pdX = m_cbinstdata.m_pdS;
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdS2[nIdx] = m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];

  m_tridiagonalMatrix.SetDimension(m_nNbX);

  CalculateAreaArrays(m_pdX, m_nNbS);
 
  // Build the sparse matrix with entries corresponding to jump terms
  BuildJumpSystem();

  // Compute the new share price first (if required) ========================
  if ( m_cbinstdata.m_bHasNewShare )
    RunNewShareStep();

  // Make the PDE coefficients for the price
  MakeCoefficients();

  // Make the alpha/beta coefficients
  BuildAlphaBeta();
  
  // Build the matrix
  BuildMatrix();
  
  // Build the right hand side for the price equation
  BuildRHS(m_cbinstdata.m_pdOldPrices.Get(),
           m_cbinstdata.m_pdOldOldPrices.Get());

  // Check if constraints are active, and choose the solver accordingly
  if ( m_cbinstdata.m_pConstraints )
    m_pIterativeSolver->Solve(m_tridiagonalMatrix,
                              m_sparseMatrix,
                              m_pdRHS.Get(),
                              m_cbinstdata.m_pConstraints,
                              m_cbinstdata.m_pdPrices.Get(),
                              m_cbinstdata.m_piFrozenFlags.Get());     
  else
    m_pLinearSolver->Solve(m_tridiagonalMatrix,
                           m_sparseMatrix,
                           m_pdRHS.Get(),
                           m_cbinstdata.m_pdPrices.Get());
}

void CBStepper::MakeCoefficients()
{
  // Get the various interest rates. These usually change at each timestep
  double dRate = m_cbinstdata.m_dRate;
  double dForeignRate = m_cbinstdata.m_dForeignRate; 
  double dDerivativeRate = m_cbinstdata.m_dDerivativeRate;

  /*
  // we have only two cases now. so we can use if{}else{}. if, in the future,
  // we get more condition check to do, we'd better use m_pdHROfDerivative
  // array in cbinstdata to avoid "if" check. In this way, we should 
  // initialize the array according to the value of m_bDerivativeHasOwnHR
  // in UpdateBeforeStep().
  if(m_cbinstdata.m_bDerivativeHasOwnHR)
  {
    double
      dSpeedCorrection = m_cbinstdata.GetSpeedCorrection();

    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = 0.5 * m_cbinstdata.m_pdVolsSquared[nIdx];
      m_pdCoe1st[nIdx] = dRate - dForeignRate + pdHazardRates[nIdx] 
        - dSpeedCorrection;

      m_pdCoeZero[nIdx] = dDerivativeRate + m_cbinstdata.m_dHROfDerivative;
      m_pdCoeConst[nIdx] = m_cbinstdata.m_dHROfDerivative *
                          m_cbinstdata.m_pdRecoveryValues[nIdx];
    }
  }
  else */
  {
    double dSpeedCorrection = m_cbinstdata.GetSpeedCorrection();

    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
      m_pdCoe1st[nIdxR] = m_pdCoe1st0[nIdxR] + dRate - dForeignRate  
                        - dSpeedCorrection;
      m_pdCoeZero[nIdxR] = m_pdCoeZero0[nIdxR] + dDerivativeRate;

      double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;
      double dIntensity = m_pdDefaultIntensities[nIdxR];
      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] = dIntensity * m_cbinstdata.m_pdRecoveryValues[nIdx];
    }
  }
}

void CBStepper::MakeNewShareCoefficients()
{
  double dRate = m_cbinstdata.m_dRate;
  double dForeignRate = m_cbinstdata.m_dForeignRate; //borrow rate

  double dSpeedCorrection = m_cbinstdata.GetSpeedCorrection();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
    m_pdCoe1st[nIdxR] = m_pdCoe1st0[nIdxR] + dRate - dForeignRate  
                      - dSpeedCorrection;
    m_pdCoeZero[nIdxR] = m_pdCoeZero0[nIdxR] + dRate;
  }

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = 0.;
}

void CBStepper::RunNewShareStep()
{
  MakeNewShareCoefficients();    

  // Make the alpha/beta coefficients 
  BuildAlphaBeta();

  // Build the matrix
  BuildMatrix(); 

  BuildRHS(m_cbinstdata.m_pdOldNewSharePrices.Get(), 
           m_cbinstdata.m_pdOldOldNewSharePrices.Get()); 

  m_pLinearSolver->Solve(m_tridiagonalMatrix,
                         m_sparseMatrix, 
                         m_pdRHS.Get(),
                         m_cbinstdata.m_pdNewSharePrices.Get());
  
  if ( m_cbinstdata.m_pConstraints )
    m_cbinstdata.UpdateConstraints(); 
}


} // namespace hg

} // namespace ito33
