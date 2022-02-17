/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/varianceswap/varianceswapstepper.cpp
// Purpose:     variance swap stepper class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapstepper.cpp,v 1.2 2006/04/10 11:38:56 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/tridiagonalsolver.h"

#include "ihg/varianceswapstepper.h"

using ito33::numeric::ComputeGammaFD;

using namespace ito33::ihg;

void VarianceSwapStepper::Init(void)
{
  m_nNbX = m_instdata.m_nNbS;

  Alloc(m_nNbX);

  CalculateAreaArrays(m_instdata.m_pdS, m_nNbX);

  m_pdX = m_instdata.m_pdS;

  m_pLinearSolver = AutoPtr<numeric::TridiagonalSolver>
                    (
                      new numeric::TridiagonalSolver(m_nNbX)
                    );
}

void VarianceSwapStepper::Run()
{
  // Make the PDE coefficients
  MakeCoefficients();

  // Make the alpha/beta coefficients
  BuildAlphaBeta();
  
  BuildMatrix();

  // Build the discrete system
  BuildRHS(m_instdata.m_pdOldPrices.Get(),
           m_instdata.m_pdOldOldPrices.Get());

  m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                         m_pdRHS.Get(),
                         m_instdata.m_pdPrices.Get());
  
  // Now that solving is done, save values in case Crank-Nicolson is used
  swap(m_pdCoeConst, m_pdCoeConstOld);
  swap(m_pdAlpha, m_pdAlphaOld);
  swap(m_pdBeta, m_pdBetaOld);
  swap(m_pdCoeZero, m_pdCoeZeroOld);

  // Variance swaps are path dependent. Greeks are computed by perturbation.

}

void VarianceSwapStepper::MakeCoefficients()
{
  size_t nIdx;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  // Get the hazard rates
  double* pdHazardRates = m_instdata.m_pdHazardRates.Get();

  for (nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdCoe2nd[nIdx] = 0.5 * m_instdata.m_pdVolsSquared[nIdx];
    m_pdCoe1st[nIdx] = dRate - dForeignRate + pdHazardRates[nIdx];
    m_pdCoeZero[nIdx] = dRate + pdHazardRates[nIdx];
    m_pdCoeConst[nIdx] = pdHazardRates[nIdx] * m_instdata.m_dRecoveryValue;
  }

  for (nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdCoe2nd[nIdx] *= m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
    m_pdCoe1st[nIdx] *= m_instdata.m_pdS[nIdx];
  }
}
