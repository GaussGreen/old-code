/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/option/optionstepper.cpp
// Purpose:     option stepper including Greek computations
// Author:      David, ZHANG Yunzhi
// Created:     2003/12/17
// RCS-ID:      $Id: optionstepper.cpp,v 1.29 2006/06/13 15:42:34 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ihg/optionstepper.h"

#include "ito33/numeric/deltagamma.h"

using ito33::numeric::ComputeGamma;
using ito33::numeric::ComputeGammaFD;

using ito33::ihg::OptionStepper;

void OptionStepper::Run()
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

  // Check if constraints are active, and choose the solver accordingly
  if (m_instdata.m_pConstraints == 0)
    m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                           m_pdRHS.Get(),
                           m_instdata.m_pdPrices.Get());
  else
  {
    // Copy m_pdOldprice into m_pdPrices.
    // This is done in order to have consistent valid
    // initial values for the matrix solve.
    // Note that for greeks the copy is not necessary since
    // we are doing a direct solve
    memcpy( m_instdata.m_pdPrices.Get(), m_instdata.m_pdOldPrices.Get(),
            m_instdata.m_nNbS);

    m_pIterativeSolver->Solve(m_tridiagonalMatrix,
                              m_pdRHS.Get(),
                              m_instdata.m_pConstraints,
                              m_instdata.m_pdPrices.Get(),
                              m_instdata.m_piFrozenFlags.Get());
  }


  // Save the constant PDed coefficients in case Crank-Nicolson is used.
  // Do not swap alpha, beta and CoeZero, since they are used by the vega PDE.
  swap(m_pdCoeConst, m_pdCoeConstOld);

  // Calculate vega if required. Note that gamma is needed.
  // Pricing data comes from the solution computed above
  if ( m_instdata.m_bComputeVega )
  {
    MakeVegaCoefficients();

    // Temporarily swap the old const vega values before building the RHS.
    // After building the RHS, swap back, and save the vega values
    swap( m_pdCoeConstOld, m_pdVegaCoeConstOld);
    BuildRHS(m_instdata.m_pdOldVegas.Get(), m_instdata.m_pdOldOldVegas.Get());
    swap( m_pdCoeConstOld, m_pdVegaCoeConstOld );
    swap( m_pdCoeConst, m_pdVegaCoeConstOld );

    if (m_instdata.m_pConstraints == 0)
      m_pLinearSolver->Solve(m_tridiagonalMatrix,
                             m_pdRHS.Get(),
                             m_instdata.m_pdVegas.Get());
    else
      m_pIterativeSolver->SolveGreek(m_tridiagonalMatrix,
                                m_pdRHS.Get(),
                                m_instdata.m_piFrozenFlags.Get(),
                                m_instdata.m_pdVegas.Get());

  }

  // Now that solving is done, save values in case Crank-Nicolson is used
  swap(m_pdAlpha, m_pdAlphaOld);
  swap(m_pdBeta, m_pdBetaOld);
  swap(m_pdCoeZero, m_pdCoeZeroOld);
}


void OptionStepper::MakeCoefficients()
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

void OptionStepper::MakeVegaCoefficients()
{
  ComputeGammaFD(m_pdX, m_instdata.m_pdPrices.Get(), 
                 m_nNbX, m_instdata.m_pdGammas.Get());

  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for pricing
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdCoeConst[nIdx] = m_instdata.m_pdVols[nIdx] * m_instdata.m_pdGammas[nIdx]
                       * m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
}
