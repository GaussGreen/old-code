/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/cboptionstepper.cpp
// Purpose:     cb option stepper
// Author:      Nabil
// Created:     2005/10/17
// RCS-ID:      $Id: cboptionstepper.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ihg/cboptionstepper.h"

#include "ito33/numeric/boundary1d.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/schemetype.h"

#include "ito33/numeric/tridiagonalfrozensolver.h"
#include "ito33/numeric/tridiagonalpenaltysolver.h"

using namespace ito33::numeric;

namespace ito33
{

namespace ihg
{

CBOptionStepper::CBOptionStepper
(CBOptionInstData& cboptioninstdata, const finance::ComputationalFlags& flags) 
                : Stepper(cboptioninstdata, flags), 
                  m_cboptioninstdata(cboptioninstdata)
{
  m_pCBStepper = shared_ptr<CBStepper>
    ( new CBStepper(*cboptioninstdata.GetCBInstData(), flags) );
}

void CBOptionStepper::Init()
{
  // Init for the underlying cb
  m_pCBStepper->Init();

  // Some allocations.

  m_nNbX = m_cboptioninstdata.GetNbSpotsMax();

  Alloc(m_nNbX);

  // Also need to save some fugit values for Crank-Nicolson
  m_pdFugitCoeZeroOld = Array<double>(m_nNbX);
  m_pdFugitCoeConstOld = Array<double>(m_nNbX);

  m_pdDiagonalOfMatrix = Array<double>(m_nNbX);
  
  // Allocate the working vector
  m_pdTmp1.resize(m_nNbX);

  m_pLinearSolver = AutoPtr<numeric::TridiagonalSolver>
                    (
                      new numeric::TridiagonalSolver(m_nNbX)
                    );
  
  if ( m_iSolverType )
  {
    m_pIterativeSolver = AutoPtr<TridiagonalConstraintSolver>
                         ( new TridiagonalFrozenSolver(m_nNbX) );
  }
  else
  {
    m_pIterativeSolver = AutoPtr<TridiagonalConstraintSolver>
                         ( new TridiagonalPenaltySolver(m_nNbX) );
  }

  // Some initializations.
  m_pdX = m_cboptioninstdata.m_pdS;
}

void CBOptionStepper::Run()
{
  // Run for the underlying cb
  m_pCBStepper->Run();
  
  if( !m_cboptioninstdata.IsCBOptionWindow() )
    return;

  // Update some data
  m_nNbX = m_cboptioninstdata.m_nNbS;
  m_tridiagonalMatrix.SetDimension(m_nNbX);
  CalculateAreaArrays(m_pdX, m_nNbX);

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
  
  //Compute the cb option constraints to be able to compute the cb 
  //option price
  m_cboptioninstdata.UpdateCBOptionConstraint();
   
  // Copy m_pdOldprice into m_pdPrices.
  // This is done in order to have consistent valid
  // initial values for the matrix solve.
  // Note that for greeks the copy is not necessary since
  // we are doing a direct solve
  memcpy( m_cboptioninstdata.m_pdPrices.Get(), 
    m_cboptioninstdata.m_pdOldPrices.Get(), m_nNbX);

  m_pIterativeSolver->Solve(m_tridiagonalMatrix,
                            m_pdRHS.Get(),
                            m_cboptioninstdata.GetCBOptionConstraint(),
                            m_cboptioninstdata.m_pdPrices.Get(),
                            m_cboptioninstdata.m_piFrozenFlags.Get());

  // Save the cb option price coefficient constant for the next step.
  swap(m_pdCoeConst, m_pdCoeConstOld);
  
  // Price data comes from the solution computed above
  if( m_cboptioninstdata.m_bComputeVega )
  {
    // Remark: Only CoeConst changes compare to the cb price

    // Some swap to prepare the vega computation.
    swap(m_pdCoeConst, m_pdCoeConstOld);
    swap(m_pdCoeConstOld, m_pdVegaCoeConstOld);
    
    MakeVegaCoefficients();

    BuildRHS(m_cboptioninstdata.m_pdOldVegas.Get(), 
             m_cboptioninstdata.m_pdOldOldVegas.Get());

    m_pIterativeSolver->SolveGreekWithSpecialConstraint
                        ( m_tridiagonalMatrix,
                          m_pdRHS.Get(),
                          m_cboptioninstdata.m_piFrozenFlags.Get(),
                          //vega of the cb
                          m_pCBStepper->m_cbinstdata.m_pdVegas.Get(),
                          m_cboptioninstdata.m_pdVegas.Get() );
    
    // Save the cb option vega coefficient constant for the next step.
    swap(m_pdCoeConst, m_pdVegaCoeConstOld);
  }

  // Calculate fugit (life) if required ============================================
 
  if( m_cboptioninstdata.m_bComputeFugit )
  {    
    // Remark: CoeZero and CoeConst changes compare to the price

    // Some swap to prepare the vega computation.
    swap(m_pdCoeConst, m_pdCoeConstOld);
    swap(m_pdCoeConstOld, m_pdFugitCoeConstOld);
    
    swap(m_pdCoeZero, m_pdCoeZeroOld);
    swap(m_pdCoeZeroOld, m_pdFugitCoeZeroOld);
    
    MakeFugitCoefficients();
    
    ChangeDiagonalOfMatrix();
    
    BuildRHS(m_cboptioninstdata.m_pdOldFugits.Get(), 
             m_cboptioninstdata.m_pdOldOldFugits.Get());
   
    if( m_cboptioninstdata.m_pConstraints == 0 )
      m_pLinearSolver->Solve(m_tridiagonalMatrix,
                             m_pdRHS.Get(),
                             m_cboptioninstdata.m_pdFugits.Get());
    else
      m_pIterativeSolver->SolveGreek(m_tridiagonalMatrix,
                                m_pdRHS.Get(),
                                m_cboptioninstdata.m_piFrozenFlags.Get(),
                                m_cboptioninstdata.m_pdFugits.Get());

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

void CBOptionStepper::MakeCoefficients()
{
  size_t nIdx;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_cboptioninstdata.m_dRate;
  double dForeignRate = m_cboptioninstdata.m_dForeignRate; 
  double dDerivativeRate = m_cboptioninstdata.m_dDerivativeRate;
  
  // Get the hazard rates
  double* pdHazardRates = m_cboptioninstdata.m_pdHazardRates.Get();

  // we have only two cases now. so we can use if{}else{}. if, in the future,
  // we get more condition check to do, we'd better use m_pdHROfDerivative
  // array in cbinstdata to avoid "if" check. In this way, we should 
  // initialize the array according to the value of m_bDerivativeHasOwnHR
  // in UpdateBeforeStep().
  if(m_cboptioninstdata.m_bDerivativeHasOwnHR)
  {
    double
      dSpeedCorrection = m_cboptioninstdata.GetSpeedCorrection();

    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = 0.5 * m_cboptioninstdata.m_pdVolsSquared[nIdx];
      m_pdCoe1st[nIdx] = dRate - dForeignRate + pdHazardRates[nIdx] 
        - dSpeedCorrection;

      m_pdCoeZero[nIdx] = dDerivativeRate 
        + m_cboptioninstdata.m_dHROfDerivative;
      m_pdCoeConst[nIdx] = 0.;
    }
  }
  else
  {
    double
      dSpeedCorrection = m_cboptioninstdata.GetSpeedCorrection();

    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = 0.5 * m_cboptioninstdata.m_pdVolsSquared[nIdx];
      m_pdCoe1st[nIdx] = dRate - dForeignRate + pdHazardRates[nIdx] 
        - dSpeedCorrection;

      m_pdCoeZero[nIdx] = dDerivativeRate + pdHazardRates[nIdx];
      m_pdCoeConst[nIdx] = 0.;
    }
  }
  
  CBInstData& 
    cbinstdata = *m_cboptioninstdata.GetCBInstData();
  
  if ( cbinstdata.IsFixedQuanto() )
  {
    double dFXVol = cbinstdata.GetFXRateVolatility();
    
    double 
      dCorrelation = cbinstdata.GetCorrelationBetweenUnderlyingAndFXRate();

    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoe1st[nIdx] -= dCorrelation * dFXVol * 
                          sqrt( m_cboptioninstdata.m_pdVolsSquared[nIdx] );
  }

  for(nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdCoe2nd[nIdx] *= m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
    m_pdCoe1st[nIdx] *= m_instdata.m_pdS[nIdx];
  }

}

void CBOptionStepper::MakeVegaCoefficients()
{
  ito33::numeric::ComputeGammaFD(m_pdX, m_cboptioninstdata.m_pdPrices.Get(), 
               m_nNbX, &m_pdTmp1[0]);

  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for the price
  for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = m_cboptioninstdata.m_pdVols[nIdx]
                       * m_pdTmp1[nIdx]
                       * m_cboptioninstdata.m_pdS[nIdx] 
                       * m_cboptioninstdata.m_pdS[nIdx];                      

  CBInstData& 
    cbinstdata = *m_cboptioninstdata.GetCBInstData();

  if ( cbinstdata.IsFixedQuanto() )
  {
    /*
    @todo Use ComputeDelta(...) which takes the finite difference scheme flag
          to compute the delta according to the finite difference scheme used 
          during discretization. 
    */
    ComputeDelta(m_pdX, cbinstdata.m_pdPrices.Get(), 
                 m_nNbX, cbinstdata.m_pdDeltas.Get());
  
    for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeConst[nIdx] -= 
        cbinstdata.GetCorrelationBetweenUnderlyingAndFXRate() *
        cbinstdata.GetFXRateVolatility() *
        cbinstdata.m_pdS[nIdx] *
        cbinstdata.m_pdDeltas[nIdx];

  }
}

void CBOptionStepper::MakeFugitCoefficients()
{
  size_t nIdx;
  
  // Get the hazard rates
  double* pdHazardRates = m_cboptioninstdata.m_pdHazardRates.Get();

  // m_pdCoe2nd and m_pdCoe1st are the same as for the price
  
  // m_pdCoeZero
  if(m_cboptioninstdata.m_bDerivativeHasOwnHR)
    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeZero[nIdx] = m_cboptioninstdata.m_dHROfDerivative;
  else
    for(nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdCoeZero[nIdx] = pdHazardRates[nIdx];
  
  // m_pdCoeConst
  for(nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = 1.;

}

void CBOptionStepper::ChangeDiagonalOfMatrix()
{
  const Boundary1D 
    &boundaryCondition = m_cboptioninstdata.GetBoundaryCondition();
  
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
    pdB[0] = m_cboptioninstdata.m_dTimeWeight + dOneMinusTheta * (m_pdBeta[0] + 
      m_pdCoeZero[0]);
   
  else
  {
    ASSERT_MSG(false, "Unknown boundary type");
  }

  // handle all interior points. Alpha and beta should already be constructed.
  // m_pdCoeZero entries change. The matrix should be an M-matrix
  for (nIdx=1; nIdx < m_nNbX-1; nIdx++)
    pdB[nIdx] = m_cboptioninstdata.m_dTimeWeight + 
                dOneMinusTheta
                  * (m_pdAlpha[nIdx] + m_pdBeta[nIdx] + m_pdCoeZero[nIdx]);

  // handle right boundary
  if (boundaryCondition.GetRightType() == BCType_Dirichlet)
    pdB[m_nNbX-1] = 1.0;

  else if (boundaryCondition.GetRightType() == BCType_Gamma)
    pdB[m_nNbX-1] = m_cboptioninstdata.m_dTimeWeight + 
                    dOneMinusTheta 
                      * (m_pdAlpha[m_nNbX-1] + m_pdCoeZero[m_nNbX-1]);

  else
  {
    ASSERT_MSG(false, "Unknown boundary type");
  }
    
}

void CBOptionStepper::RestoreDiagonalOfMatrix()
{
  double* pdB = m_tridiagonalMatrix.GetB();

  // Restore the diagonal of the matrix
  for(size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdB[nIdx] = m_pdDiagonalOfMatrix[nIdx];
}

} // namespace ihg

} // namespace ito33
