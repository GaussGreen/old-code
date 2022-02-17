/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/eds/edsstepper.cpp
// Purpose:     EDS stepper class
// Created:     2005/01/26
// RCS-ID:      $Id: edsstepper.cpp,v 1.2 2006/03/07 12:19:53 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/tridiagonalsolver.h"

#include "ihg/edsstepper.h"

using ito33::numeric::ComputeGamma;
using ito33::numeric::ComputeGammaFD;

using namespace ito33::ihg;

void EDSStepper::Init(void)
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

void EDSStepper::Run()
{
  // Compute the option price first
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

  // Save the constant PDed coefficients in case Crank-Nicolson is used.
  // Do not swap alpha and beta yet, since they are used by the Greek PDEs.
  swap(m_pdCoeConst, m_pdCoeConstOld);

  // Determine which Greeks we need to compute
  bool bComputeVega = m_instdata.m_bComputeVega;

  // Now calculate the Greeks
  if (bComputeVega)
  {
    MakeVegaCoefficients();

    // Temporarily swap the old const vega values before building the RHS.
    // After building the RHS, swap back, and save the vega values
    swap( m_pdCoeConstOld, m_pdVegaCoeConstOld);

    BuildRHS(m_instdata.m_pdOldVegas.Get(), m_instdata.m_pdOldOldVegas.Get());

    swap( m_pdCoeConstOld, m_pdVegaCoeConstOld );
    swap( m_pdCoeConst, m_pdVegaCoeConstOld );

    m_pLinearSolver->Solve(m_tridiagonalMatrix,
                           m_pdRHS.Get(),
                           m_instdata.m_pdVegas.Get());
  }

  // Now that solving is done, save values in case Crank-Nicolson is used
  swap(m_pdAlpha, m_pdAlphaOld);
  swap(m_pdBeta, m_pdBetaOld);
  swap(m_pdCoeZero, m_pdCoeZeroOld);
}

void EDSStepper::MakeCoefficients()
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

void EDSStepper::MakeVegaCoefficients()
{
  ComputeGammaFD(m_pdX, m_instdata.m_pdPrices.Get(), 
                 m_nNbX, m_instdata.m_pdGammas.Get());

  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for pricing
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = m_instdata.m_pdVols[nIdx] * m_instdata.m_pdGammas[nIdx]
                       * m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
}
