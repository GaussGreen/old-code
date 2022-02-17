/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/stepper_fix.cpp
// Purpose:     Stepper class for problem with fixed space mesh
// Created:     2005/01/31
// RCS-ID:      $Id: stepper_fix.cpp,v 1.15 2006/03/31 17:43:59 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/trisparselinearsolver_gmres.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"

#include "hg/model.h"
#include "hg/stepper_fix.h"
#include "hg/sensitivitydata.h"

#include "ito33/numeric/boundary1d.h"
#include "ito33/numeric/extrapolationmode.h"
#include "ito33/numeric/interpolation.h"

using namespace ito33::hg;
using namespace ito33::numeric;

void StepperFix::Init()
{
  Stepper::Init();

  m_nNbS = m_instdata.m_nNbS;
  
  m_pdX = m_instdata.m_pdS;
  
  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  m_pdS2 = Array<double>(m_nNbS);
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdS2[nIdx] = m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];

  CalculateAreaArrays(m_pdX, m_nNbS);

  // Build the sparse matrix with entries corresponding to jump terms
  BuildJumpSystem();

  m_pLinearSolver = AutoPtr<numeric::TriSparseLinearSolver>
                    ( new numeric::TriSparseLinearSolverGmres(m_nNbX) );
                    //( new numeric::TriSparseLinearSolverFixedPoint(m_nNbX) );

  BuildDerivedJumpSystem();

  m_instdata.m_pSparseMatrix = &m_sparseMatrix;
}

void StepperFix::BuildDerivedJumpSystem()
{
  for (size_t n = 0; n < m_instdata.m_pSensitivityData.size(); n++)
  {
    SensitivityData& data = m_instdata.m_pSensitivityData[n];
    
    if ( data.m_sensitivityType == SensitivityType_JumpIntensity )
      data.m_pSparseMatrix = BuildJumpSensitivitySystem
                             (data.m_nRegime1, data.m_nRegime2, 
                              data.m_dAmplitude, data.m_dIntensity, false);

    if ( data.m_sensitivityType == SensitivityType_JumpAmplitude )
      data.m_pSparseMatrix = BuildJumpSensitivitySystem
                             (data.m_nRegime1, data.m_nRegime2, 
                              data.m_dAmplitude, data.m_dIntensity, true);
  }
}

void StepperFix::Run()
{
  BuildMainSystem();

  m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                         m_sparseMatrix,
                         m_pdRHS.Get(),
                         m_instdata.m_pdPrices.Get());

  // Calculate sensitivities if required. Note that gamma is needed.
  // Pricing data comes from the solution computed above
  size_t nNbSensitivity = m_instdata.m_pSensitivityData.size();

  if ( nNbSensitivity 
   && m_instdata.m_sensitivityMethod == SensitivityMethod_PDE )
  {
    // Reset the boundary condition, if applicable (eg. Dircihlet type)
    m_instdata.SetSensitivityBoundary();

    // Compute helper arrays 
    ComputeHelperArrays();

    for (size_t nIdxD = 0; nIdxD < nNbSensitivity; nIdxD++)
    {      
      MakeSensitivityCoefficients(m_instdata.m_pSensitivityData[nIdxD]);

      BuildRHS(m_instdata.m_ppdOldSensitivities[nIdxD].Get(), 
               m_instdata.m_ppdOldOldSensitivities[nIdxD].Get());

      m_pLinearSolver->Solve(m_tridiagonalMatrix,
                             m_sparseMatrix,
                             m_pdRHS.Get(),
                             m_instdata.m_ppdSensitivities[nIdxD].Get());

    } // loop over the model params
  } // if computing sensitivities
  
  if ( m_instdata.m_bDualSystemRequired )
    SetupDualSystemData();

  if (  m_instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
    SetupSensivityByAdjointData();
}

void StepperFix::SetupDualSystemData()
{
  SensitivityByAdjointData& data = m_instdata.m_aData;
  data.m_bIsValid = true;
  data.m_dRecoveryValue = m_instdata.m_dRecoveryValue;
  data.m_dOldOldTimeWeight = m_instdata.m_dOldOldTimeWeight;
  data.m_dOldTimeWeight = m_instdata.m_dOldTimeWeight;
  data.m_pMatrix = AutoPtr<TridiagonalMatrix>
                   ( new TridiagonalMatrix(m_tridiagonalMatrix) );
}

void StepperFix::SetupSensivityByAdjointData()
{
  SensitivityByAdjointData& data = m_instdata.m_aData;

  Array<double> pdPrices(m_nNbX);
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdPrices[nIdx] = m_instdata.m_pdPrices[nIdx];
  data.m_pdPrices = pdPrices;

  Array<int> piFDFlags(m_nNbX);
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    piFDFlags[nIdx] = m_piFD[nIdx];

  data.m_piFDFlags = piFDFlags;
}

void StepperFix::MakeCoefficients()
{
  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
    m_pdCoe1st[nIdxR] = m_pdCoe1st0[nIdxR] + dRate - dForeignRate;
    m_pdCoeZero[nIdxR] = m_pdCoeZero0[nIdxR] + dRate;
 
    double dTmp = m_pdDefaultIntensities[nIdxR] * m_instdata.m_dRecoveryValue;

    double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;
    for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
      pdCoeConst[nIdx] = dTmp;
  }
}

void StepperFix::MakeSensitivityCoefficients(const SensitivityData& data)
{
  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for pricing.
  // Looking at the actual HG PDE should make the code below easy to
  // understand. Just derive the equation with respet to the appropriate
  // variable
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++) 
  { 
    // get the constant coefficients for this regime
    double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;

    // Typically, the sensitivity parameter only applies to one regime,
    // so all other regimes have zero entries for the constant coeff's.
    if (data.m_nRegime1 != nIdxR)
    {
      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] = 0.0;

      continue;
    }

    // Determine what type of sensitivity is being computed.
    switch (data.m_sensitivityType)
    {
    case SensitivityType_Volatility:
      {
      // volatility coefficient
      double* pdGammas = m_instdata.m_pdGammas.Get() + nIdxR * m_nNbS;

      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] = m_pdVols[nIdxR] * pdGammas[nIdx] * m_pdS2[nIdx];

      break;
      }
    case SensitivityType_DefaultIntensity:
    {
      // default intensity coefficient
      double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;
      double* pdPrices = m_instdata.m_pdPrices.Get() + nIdxR * m_nNbS;

      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] = m_instdata.m_dRecoveryValue - pdPrices[nIdx]
                         + m_pdX[nIdx] * pdDeltas[nIdx];

      break;
    }
    case SensitivityType_JumpIntensity:
      {
      // intensity
      data.m_pSparseMatrix->ProductMatrixVector
                            (m_instdata.m_pdPrices.Get(), m_pdCoeConst.Get());

      double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;
      double* pdPrices = m_instdata.m_pdPrices.Get() + nIdxR * m_nNbS;

      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] += - pdPrices[nIdx]
                            - data.m_dAmplitude * m_pdX[nIdx] * pdDeltas[nIdx];

      break;
      }
    case SensitivityType_JumpAmplitude:
      {
      // amplitude
      double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;

      data.m_pSparseMatrix->ProductMatrixVector
                            (m_instdata.m_pdPrices.Get(), m_pdCoeConst.Get());

      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] -= data.m_dIntensity * m_pdX[nIdx] * pdDeltas[nIdx];

      break;
      }
    } // end the switch

  } // loop over regimes
}

void StepperFix::BuildMainSystem()
{
  // Make the PDE coefficients
  MakeCoefficients();

  // Make the alpha/beta coefficients
  BuildAlphaBeta();
  
  BuildMatrix();

  // Build the discrete system
  BuildRHS(m_instdata.m_pdOldPrices.Get(),
           m_instdata.m_pdOldOldPrices.Get());
}
