/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/cbstepper.cpp
// Purpose:     cb stepper including greek, straight bond and fugit computations
// Author:      Nabil
// Created:     2004/04/09
// RCS-ID:      $Id: cbstepper.cpp,v 1.51 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ihg/src/cb/cbstepper.cpp
  @brief implement CBStepper class
*/

#include "ihg/cbstepper.h"

#include "ito33/numeric/boundary1d.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/schemetype.h"

#include "ito33/numeric/tridiagonalpenaltysolver.h"
#include "ito33/numeric/tridiagonalfrozensolver.h"

using namespace ito33::numeric;

namespace ito33
{

namespace ihg
{

void CBStepper::Init()
{
  // Some allocations.

  m_nNbX = m_cbinstdata.GetNbSpotsMax();

  Alloc(m_nNbX);

  // Also need to save some fugit values for Crank-Nicolson
  m_pdFugitCoeZeroOld = Array<double>(m_nNbX);
  m_pdFugitCoeConstOld = Array<double>(m_nNbX);

  if( m_cbinstdata.HasNewShare() )
  {
    // Also need to save some new share values for Crank-Nicolson
    m_pdNewShareCoeZeroOld = Array<double>(m_nNbX);
    m_pdNewShareCoeConstOld = Array<double>(m_nNbX);
  }

  m_pdDiagonalOfMatrix = Array<double>(m_nNbX);

  m_pLinearSolver = AutoPtr<TridiagonalSolver>
                    ( new TridiagonalSolver(m_nNbX) );

  if ( m_iSolverType )
    m_pIterativeSolver = AutoPtr<TridiagonalConstraintSolver>
                         ( new TridiagonalFrozenSolver(m_nNbX, m_nNbX) );
  else     
    m_pIterativeSolver = AutoPtr<TridiagonalConstraintSolver>
                         ( new TridiagonalPenaltySolver(m_nNbX, m_nNbX) );

  // Some initializations.
  m_pdX = m_cbinstdata.m_pdS;
}

void CBStepper::Run()
{
  // Update some data
  m_nNbX = m_cbinstdata.m_nNbS;
  m_tridiagonalMatrix.SetDimension(m_nNbX);
  CalculateAreaArrays(m_pdX, m_nNbX);

  // Compute the new share price first (if required) ========================
  if( m_cbinstdata.HasNewShare() )
    RunNewShareStep();

  // Compute the CB price ===================================================

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
  if( m_cbinstdata.m_pConstraints == 0 )
    m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                           m_pdRHS.Get(),
                           m_cbinstdata.m_pdPrices.Get());
  else
  {
    // Copy m_pdOldprice into m_pdPrices.
    // This is done in order to have consistent valid
    // initial values for the matrix solve.
    // Note that for greeks the copy is not necessary since
    // we are doing a direct solve
    memcpy( m_cbinstdata.m_pdPrices.Get(), m_cbinstdata.m_pdOldPrices.Get(), m_nNbX);

    m_pIterativeSolver->Solve(m_tridiagonalMatrix,
                              m_pdRHS.Get(),
                              m_cbinstdata.m_pConstraints,
                              m_cbinstdata.m_pdPrices.Get(),
                              m_cbinstdata.m_piFrozenFlags.Get());
  }
  // Calculate vega if required. Note that gamma is needed ==================

  // Pricing data comes from the solution computed above
  if( m_cbinstdata.m_bComputeVega )
  {
    // Remark: Only CoeConst changes compare to the price

    // Some swap to prepare the vega computation.
    swap(m_pdCoeConst, m_pdCoeConstOld);
    swap(m_pdCoeConstOld, m_pdVegaCoeConstOld);

    MakeVegaCoefficients();

    BuildRHS(m_cbinstdata.m_pdOldVegas.Get(), 
             m_cbinstdata.m_pdOldOldVegas.Get());

    if( m_cbinstdata.m_pConstraints == 0 )
      m_pLinearSolver->Solve(m_tridiagonalMatrix,
                             m_pdRHS.Get(),
                             m_cbinstdata.m_pdVegas.Get());
    else
      m_pIterativeSolver->SolveGreek(m_tridiagonalMatrix,
                                m_pdRHS.Get(),
                                m_cbinstdata.m_piFrozenFlags.Get(),
                                m_cbinstdata.m_pdVegas.Get());

    // Restore m_pdCoeConst and initialize m_pdVegaCoeConstOld to VegaCoeConst
    swap( m_pdCoeConst, m_pdVegaCoeConstOld );
  }

  // Calculate fugit (life) if required ============================================
 
  if( m_cbinstdata.m_bComputeFugit )
  {    
    // Remark: CoeZero and CoeConst changes compare to the price

    // Some swap to prepare the vega computation.
    swap(m_pdCoeConst, m_pdCoeConstOld);
    swap(m_pdCoeConstOld, m_pdFugitCoeConstOld);
    
    swap(m_pdCoeZero, m_pdCoeZeroOld);
    swap(m_pdCoeZeroOld, m_pdFugitCoeZeroOld);
    
    MakeFugitCoefficients();
    
    ChangeDiagonalOfMatrix();
    
    BuildRHS(m_cbinstdata.m_pdOldFugits.Get(), 
             m_cbinstdata.m_pdOldOldFugits.Get());
   
    if( m_cbinstdata.m_pConstraints == 0 )
      m_pLinearSolver->Solve(m_tridiagonalMatrix,
                             m_pdRHS.Get(),
                             m_cbinstdata.m_pdFugits.Get());
    else
      m_pIterativeSolver->SolveGreek(m_tridiagonalMatrix,
                                m_pdRHS.Get(),
                                m_cbinstdata.m_piFrozenFlags.Get(),
                                m_cbinstdata.m_pdFugits.Get());

    // Restore m_pdCoeConst and m_pdCoeZero and initialize m_pdFugitCoeConstOld
    //to FugitCoeConst and m_pdFugitCoeZeroOld to FugitCoeZero.
    swap( m_pdCoeConst, m_pdFugitCoeConstOld );
    swap( m_pdCoeZero, m_pdFugitCoeZeroOld );
    
    RestoreDiagonalOfMatrix();
  }

  // Now that solving is done, save values in case Crank-Nicolson is used
  swap(m_pdCoeConst, m_pdCoeConstOld);
  swap(m_pdCoeZero, m_pdCoeZeroOld);
  swap(m_pdAlpha, m_pdAlphaOld);
  swap(m_pdBeta, m_pdBetaOld);
  
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
                          m_pdRHS.Get(),
                          m_cbinstdata.m_pdNewSharePrices.Get());
  
  if( m_cbinstdata.m_pConstraints != 0 )
    m_cbinstdata.UpdateConstraints(); 

  // TODO: Treatment for Crank-Nicolson:
  // Initialize m_pdNewShareCoeConstOld to NewShareCoeConst and 
  // m_pdNewShareCoeZeroOld to NewShareCoeZero...
  /*swap( m_pdCoeConst, m_pdNewShareCoeConstOld );
  swap( m_pdCoeZero, m_pdNewShareCoeZeroOld );
  swap( m_pdAlpha, m_pdNewShareAlphaOld );
  swap( m_pdBeta, m_pdNewShareBetaOld );*/ 
}

void CBStepper::MakeCoefficients()
{
  size_t nIdx;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_cbinstdata.m_dRate;
  double dForeignRate = m_cbinstdata.m_dForeignRate; 
  double dDerivativeRate = m_cbinstdata.m_dDerivativeRate;
  
  // Get the hazard rates
  double* pdHazardRates = m_cbinstdata.m_pdHazardRates.Get();

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
  else
  {
    double
      dSpeedCorrection = m_cbinstdata.GetSpeedCorrection();

    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = 0.5 * m_cbinstdata.m_pdVolsSquared[nIdx];
      m_pdCoe1st[nIdx] = dRate - dForeignRate + pdHazardRates[nIdx] 
        - dSpeedCorrection;

      m_pdCoeZero[nIdx] = dDerivativeRate + pdHazardRates[nIdx];
      m_pdCoeConst[nIdx] = pdHazardRates[nIdx] *
                          m_cbinstdata.m_pdRecoveryValues[nIdx];
    }
  }

  if ( m_cbinstdata.IsFixedQuanto() )
  {
    double dFXVol = m_cbinstdata.GetFXRateVolatility();
    
    double 
      dCorrelation = m_cbinstdata.GetCorrelationBetweenUnderlyingAndFXRate();

    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoe1st[nIdx] -= dCorrelation * dFXVol * 
                          sqrt( m_cbinstdata.m_pdVolsSquared[nIdx] );
  }

  for(nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdCoe2nd[nIdx] *= m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
    m_pdCoe1st[nIdx] *= m_instdata.m_pdS[nIdx];
  }

}

void CBStepper::MakeNewShareCoefficients()
{
  double dRate = m_cbinstdata.m_dRate;
  double dForeignRate = m_cbinstdata.m_dForeignRate; //borrow rate

  // Get the hazard rates (new share dosn't work with exchangeable)
  double* pdHazardRates = m_cbinstdata.m_pdHazardRates.Get();
    
  double dSpeedCorrection = m_cbinstdata.GetSpeedCorrection();

  for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdCoe2nd[nIdx] = 0.5 * m_cbinstdata.m_pdVolsSquared[nIdx];
    m_pdCoe1st[nIdx] = dRate - dForeignRate + pdHazardRates[nIdx] 
                     - dSpeedCorrection;
    m_pdCoeZero[nIdx] = dRate + pdHazardRates[nIdx];
    m_pdCoeConst[nIdx] = 0.;
    
    m_pdCoe2nd[nIdx] *= m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
    m_pdCoe1st[nIdx] *= m_instdata.m_pdS[nIdx];
  }
}

void CBStepper::MakeVegaCoefficients()
{
  ComputeGammaFD(m_pdX, m_cbinstdata.m_pdPrices.Get(), 
                 m_nNbX, m_cbinstdata.m_pdGammas.Get());

  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for pricing
  for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = m_cbinstdata.m_pdVols[nIdx]
                       * m_cbinstdata.m_pdGammas[nIdx]
                       * m_cbinstdata.m_pdS[nIdx] * m_cbinstdata.m_pdS[nIdx];                       

  if ( m_cbinstdata.IsFixedQuanto() )
  {
    /*
    @todo Use ComputeDelta(...) which takes the finite difference scheme flag
          to compute the delta according to the finite difference scheme used 
          during discretization. 
    */
    ComputeDelta(m_pdX, m_cbinstdata.m_pdPrices.Get(), 
                 m_nNbX, m_cbinstdata.m_pdDeltas.Get());
  
    for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeConst[nIdx] -= 
        m_cbinstdata.GetCorrelationBetweenUnderlyingAndFXRate() *
        m_cbinstdata.GetFXRateVolatility() *
        m_cbinstdata.m_pdS[nIdx] *
        m_cbinstdata.m_pdDeltas[nIdx];

  }
}

void CBStepper::MakeFugitCoefficients()
{
  size_t nIdx;
  
  // Get the hazard rates
  double* pdHazardRates = m_cbinstdata.m_pdHazardRates.Get();

  // m_pdCoe2nd and m_pdCoe1st are the same as for the price
  
  // m_pdCoeZero
  if(m_cbinstdata.m_bDerivativeHasOwnHR)
    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeZero[nIdx] = m_cbinstdata.m_dHROfDerivative;
  else
    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeZero[nIdx] = pdHazardRates[nIdx];
  
  // m_pdCoeConst
  for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = 1.;

}

void CBStepper::ChangeDiagonalOfMatrix()
{
  const Boundary1D 
    &boundaryCondition = m_instdata.GetBoundaryCondition();
  
  double
    dOneMinusTheta = 1.0;

  if(m_instdata.m_schemeType == SchemeType_CrankNicolson)
    dOneMinusTheta = 0.5;

  // Get the diagonal of the tridiagonal matrix for easy access.
  // Only the diagonal change.
  double* pdB = m_tridiagonalMatrix.GetB();

  size_t nIdx;
  // Save the diagonal of the matrix
  for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdDiagonalOfMatrix[nIdx] = pdB[nIdx];

  // Construct the new diagonal of the matrix
  
  // Handle left boundary
  if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
    pdB[0] = 1.0;

  else if (boundaryCondition.GetLeftType() == BCType_Gamma)
    pdB[0] = m_instdata.m_dTimeWeight + dOneMinusTheta * (m_pdBeta[0] + 
      m_pdCoeZero[0]);
   
  else
  {
    ASSERT_MSG(false, "Unknown boundary type");
  }

  // handle all interior points. Alpha and beta should already be constructed.
  // m_pdCoeZero entries change. The matrix should be an M-matrix
  for (nIdx=1; nIdx < m_nNbX-1; nIdx++)
    pdB[nIdx] = m_instdata.m_dTimeWeight + 
                dOneMinusTheta
                  * (m_pdAlpha[nIdx] + m_pdBeta[nIdx] + m_pdCoeZero[nIdx]);

  // handle right boundary
  if (boundaryCondition.GetRightType() == BCType_Dirichlet)
    pdB[m_nNbX-1] = 1.0;

  else if (boundaryCondition.GetRightType() == BCType_Gamma)
    pdB[m_nNbX-1] = m_instdata.m_dTimeWeight + 
                    dOneMinusTheta 
                      * (m_pdAlpha[m_nNbX-1] + m_pdCoeZero[m_nNbX-1]);

  else
  {
    ASSERT_MSG(false, "Unknown boundary type");
  }
    
}

void CBStepper::RestoreDiagonalOfMatrix()
{
  double* pdB = m_tridiagonalMatrix.GetB();

  // Restore the diagonal of the matrix
  for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdB[nIdx] = m_pdDiagonalOfMatrix[nIdx];
}

} // namespace ihg

} // namespace ito33
