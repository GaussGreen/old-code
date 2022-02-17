/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/stepper_fe_fix.cpp
// Purpose:     Stepper class for problem with fixed space mesh
// Created:     2005/01/31
// RCS-ID:      $Id: stepper_fe_fix.cpp,v 1.21 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/trisparselinearsolver_gmres.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"

#include "hg/stepper_fe_fix.h"
#include "hg/sensitivitydata.h"

using namespace ito33::hg;
using namespace ito33::numeric;

void StepperFEFix::Init()
{
  StepperFE::Init();

  m_nNbS = m_instdata.m_nNbS;
  
  // We are using logS as variable
  m_pdX = m_instdata.m_pdLogS;
  
  m_nNbX = m_nNbS * m_nNbRegimes;

  Alloc(m_nNbX);

  for (size_t nIdx = 0; nIdx < m_nNbS - 1; nIdx++)
  {
    m_pdDeltaX[nIdx] = m_pdX[nIdx + 1] - m_pdX[nIdx];
    m_pdInverseDeltaX[nIdx] = 1. / m_pdDeltaX[nIdx];
  }

  BuildMassMatrix();

  // Build the sparse matrix with entries corresponding to jump terms
  BuildJumpSystem();

  m_pLinearSolver = AutoPtr<TriSparseLinearSolver>
                    ( new TriSparseLinearSolverGmres(m_nNbX) );
                    //( new TriSparseLinearSolverFixedPoint(m_nNbX) );

  BuildDerivedJumpSystem();

  m_instdata.m_pSparseMatrix = &m_sparseMatrix;
  m_instdata.m_pMassMatrix = &m_massMatrix;
}

void StepperFEFix::BuildDerivedJumpSystem()
{
  for (size_t n = 0; n < m_instdata.m_pSensitivityData.size(); n++)
  {
    SensitivityData& data = m_instdata.m_pSensitivityData[n];
    
    const size_t nIdxR1 = data.m_nRegime1;
    const size_t nIdxR2 = data.m_nRegime2;
    const double dIntensity = data.m_dIntensity;
    const double dAmplitude = data.m_dAmplitude;

    const size_t nOffSetRow = nIdxR1 * m_nNbS;
    Array< std::list<size_t> > ppnColumnLists(m_nNbX);
    
    for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
    {  
      const size_t nIdxRow = nOffSetRow + nIdxS;
      ppnColumnLists[nIdxRow].push_back(nIdxRow - 1);
      ppnColumnLists[nIdxRow].push_back(nIdxRow);
      ppnColumnLists[nIdxRow].push_back(nIdxRow + 1);
    }

    data.m_pSparseMatrix = AutoPtr<MorseMatrix>(new MorseMatrix);
    MorseMatrix& sparseMatrix = *(data.m_pSparseMatrix);

    // Determine what type of sensitivity is being computed.
    switch (data.m_sensitivityType)
    {
    case SensitivityType_Volatility:
      {
      shared_ptr<MorseStruct> 
        pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );      
      sparseMatrix.Init(pMS);

      const double dVol = m_pdVols[nIdxR1];

      for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
      {  
        const size_t nIdxRow = nOffSetRow + nIdxS;

        double dInvDeltaX0 = m_pdInverseDeltaX[nIdxS - 1];
        double dInvDeltaX = m_pdInverseDeltaX[nIdxS];   
        
        sparseMatrix(nIdxRow, nIdxRow - 1) = dVol * dInvDeltaX0 + 0.5 * dVol;      
        sparseMatrix(nIdxRow, nIdxRow) = - dVol * (dInvDeltaX0 + dInvDeltaX); 
        sparseMatrix(nIdxRow, nIdxRow + 1) = dVol * dInvDeltaX - 0.5 * dVol; 
      }

      break;
      }
    case SensitivityType_DefaultIntensity:
      { 

      shared_ptr<MorseStruct> 
        pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );      
      sparseMatrix.Init(pMS);

      for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
      {  
        const size_t nIdxRow = nOffSetRow + nIdxS;  
        
        double dDeltaX0 = m_pdDeltaX[nIdxS - 1];  
        double dDeltaX = m_pdDeltaX[nIdxS];  

        sparseMatrix(nIdxRow, nIdxRow - 1) = - 0.5 - dDeltaX0 / 6.;      
        sparseMatrix(nIdxRow, nIdxRow) = - (dDeltaX0 + dDeltaX) / 3.; 
        sparseMatrix(nIdxRow, nIdxRow + 1) = 0.5 - dDeltaX / 6.;  
      }

      break;
      }
    case SensitivityType_JumpIntensity:
      {
      BuildSensitivityJumpSystemStructure
      ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

      shared_ptr<MorseStruct> 
        pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );      
      sparseMatrix.Init(pMS, 0);
      
      for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
      {  
        const size_t nIdxRow = nOffSetRow + nIdxS;  
        
        double dDeltaX0 = m_pdDeltaX[nIdxS - 1];  
        double dDeltaX = m_pdDeltaX[nIdxS];  

        sparseMatrix(nIdxRow, nIdxRow - 1) = - 0.5 * (- dAmplitude) - dDeltaX0 / 6.;      
        sparseMatrix(nIdxRow, nIdxRow) = - (dDeltaX0 + dDeltaX) / 3.; 
        sparseMatrix(nIdxRow, nIdxRow + 1) = 0.5 * (- dAmplitude) - dDeltaX / 6.;  
      }
      
      BuildPartialSensitivityJumpSystem
      ( nIdxR1, nIdxR2, dAmplitude, 1, sparseMatrix );

      break;
      }
    case SensitivityType_JumpAmplitude:
      {
      BuildSensitivityJumpSystemStructure
      ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

      shared_ptr<MorseStruct> 
        pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );      
      sparseMatrix.Init(pMS, 0);
      
      for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
      {  
        const size_t nIdxRow = nOffSetRow + nIdxS;  
        
        sparseMatrix(nIdxRow, nIdxRow - 1) = - 0.5 * (- dIntensity);      
        sparseMatrix(nIdxRow, nIdxRow) = 0.;
        sparseMatrix(nIdxRow, nIdxRow + 1) = 0.5 * (- dIntensity);  
      }
      
      BuildAmplitudeSensitivityJumpSystem
      ( nIdxR1, nIdxR2, dAmplitude, dIntensity, sparseMatrix);

      break;
      }
   } // end the switch
  }
}

void StepperFEFix::BuildMainSystem()
{
  // Make the PDE coefficients
  MakeCoefficients();

  BuildMatrix();

  // Build the discrete system
  BuildRHS(m_instdata.m_pdOldPrices.Get(),
           m_instdata.m_pdOldOldPrices.Get());
}

void StepperFEFix::Run()
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
    &&  m_instdata.m_sensitivityMethod == SensitivityMethod_PDE )
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

void StepperFEFix::SetupDualSystemData()
{
  SensitivityByAdjointData& data = m_instdata.m_aData;
  data.m_bIsValid = true;
  data.m_dRecoveryValue = m_instdata.m_dRecoveryValue;
  data.m_dOldOldTimeWeight = m_instdata.m_dOldOldTimeWeight;
  data.m_dOldTimeWeight = m_instdata.m_dOldTimeWeight;
  data.m_pMatrix = AutoPtr<TridiagonalMatrix>
                   ( new TridiagonalMatrix(m_tridiagonalMatrix) );
}

void StepperFEFix::SetupSensivityByAdjointData()
{
  SensitivityByAdjointData& data = m_instdata.m_aData;

  Array<double> pdPrices(m_nNbX);
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdPrices[nIdx] = m_instdata.m_pdPrices[nIdx];

  data.m_pdPrices = pdPrices;
}

void StepperFEFix::MakeCoefficients()
{
  m_bSensitivityRHS = false;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
    m_pdCoe1st[nIdxR] = m_pdCoe2nd0[nIdxR] + m_pdCoe1st0[nIdxR] 
                      + dForeignRate - dRate;
    m_pdCoeZero[nIdxR] = m_pdCoeZero0[nIdxR] + dRate + m_instdata.m_dTimeWeight;
 
    double dTmp = m_pdDefaultIntensities[nIdxR] * m_instdata.m_dRecoveryValue;

    double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;
    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      pdCoeConst[nIdxS] = dTmp;
  }
}

void StepperFEFix::MakeSensitivityCoefficients(const SensitivityData& data)
{
  m_bSensitivityRHS = true;

  for (size_t nIdxX = 0; nIdxX < m_nNbX; nIdxX++)
  {
    m_pdCoeConst[nIdxX] = 0.;
  }

  data.m_pSparseMatrix->ProductMatrixVector
                        (m_instdata.m_pdPrices.Get(), 
                         m_pdCoeConstSensitivity.Get());
  
  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for pricing.
  // Looking at the actual HG PDE should make the code below easy to
  // understand. Just derive the equation with respect to the appropriate
  // variable

  // Needs only to take care of the values corresponding to iRegime1
  const size_t nOffset = data.m_nRegime1 * m_nNbS;
  
  double* pdCoeConst = m_pdCoeConst.Get() + nOffset;

  // Determine what type of sensitivity is being computed.
  if (data.m_sensitivityType == SensitivityType_DefaultIntensity)
  {
    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      pdCoeConst[nIdxS] = m_instdata.m_dRecoveryValue;
  }
}
