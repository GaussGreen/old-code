/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/sensitivitybyadjoint.cpp
// Purpose:     implementation of BackwardNumOutput class 
// Created:     2005/05/30
// RCS-ID:      $Id: sensitivitybyadjoint.cpp,v 1.18 2006/07/30 18:21:21 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
  @todo: American option might need to pay more attention. Strange enough that 
         when fd is used, PDE solving gives worse result when there is 
         dividend.

         We are supposing that there is at least a few points(2 or 3?) between
         two events which is reasonable but still might not be true in extreme
         case.

         Also, if there is event at observation date, the implementation needs
         to be refined and backward and forward will be different.
*/

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/morsematrix.h"
#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"

#include "hg/backwardnumoutput.h"
#include "hg/backwardinstdata.h"
#include "hg/model.h"

namespace ito33
{

namespace hg
{
  using namespace numeric;

void BackwardNumOutput::SolveTransposeSystem
      (BackwardInstData& instdata, const double* pdRHS, double* pdValues)
{
  if ( m_computationalFlags.GetDiscretizationMethod() )
    SolveTransposeSystemFE(instdata, pdRHS, pdValues);
  else
    SolveTransposeSystemFD(instdata, pdRHS, pdValues);
}

void BackwardNumOutput::SensitivityByAdjoint(BackwardInstData& instdata)
{
  if ( m_computationalFlags.GetDiscretizationMethod() )
    SensitivityByAdjointFE(instdata);
  else
    SensitivityByAdjointFD(instdata);
}

void BackwardNumOutput::SolveTransposeSystemFD
     (BackwardInstData& instdata, const double* pdPrices0, double* pdPrices)
{
  const MorseMatrix* pSparseMatrix = instdata.m_pSparseMatrix;

  Array<double> pdAPrices(m_nNbX);
  Array<double> pdOldAPrices(m_nNbX);
  Array<double> pdOldOldAPrices(m_nNbX);
  
  Array<double> pdRHS(m_nNbX);
  Array<double> pdTmp(m_nNbX);
  Array<double> pdS2(m_nNbS);
  for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    pdS2[nIdxS] = m_pdS[nIdxS] * m_pdS[nIdxS];

  TriSparseLinearSolverGmres solver(m_nNbX);
  solver.EnableTransposeMatrixSolving();

  // setup the initial guess
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdAPrices[nIdx] = pdOldAPrices[nIdx] = pdOldOldAPrices[nIdx] = pdPrices0[nIdx];

  size_t nNbSystems = m_pAdjointDatas.size();

  // UpdateMe updates for initial values which is not associated to a system
  for (size_t n = nNbSystems - 1; n < nNbSystems; n--)
  {    
    swap(pdOldOldAPrices, pdAPrices);
    swap(pdOldOldAPrices, pdOldAPrices);

    const TridiagonalMatrix* pMatrix = m_pAdjointDatas[n].m_pMatrix.get();

    if ( n == nNbSystems - 1 )
    {
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        pdRHS[nIdx] = pdAPrices[nIdx];
    }
    else
    {
      double dOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if (  n + 2 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 2].m_dOldOldTimeWeight;

      if ( n + 1 < nNbSystems && m_pAdjointDatas[n + 1].m_pInterpMatrix )
        dOldOldTimeWeight = 0.;

      const InterpolationMatrix* 
        pInterpMatrix = m_pAdjointDatas[n].m_pInterpMatrix.get();

      if ( pInterpMatrix )
      {
        pInterpMatrix->ProductTransposeMatrixVector
                       ( pdOldAPrices.Get(), pdRHS.Get() );

        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight; 

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdRHS[nIdx] += pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
      }
      else
      {
        if ( dOldOldTimeWeight == 0 )
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdRHS[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
        else
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdRHS[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                        + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
      }

      instdata.ApplyBoundaryConditionToSensitivityRHS( pdRHS.Get() );
    }

    if ( m_pAdjointDatas[n].m_pRMatrix )
      m_pAdjointDatas[n].m_pRMatrix->ProductTransposeMatrixVector
          (m_pAdjointDatas[n].m_pdResiduals.Get(), pdRHS.Get(), true);

    solver.Solve( *pMatrix, *pSparseMatrix, pdRHS.Get(), pdAPrices.Get() );

    // Dividend caused interpolation
    if ( n > 0 && m_pAdjointDatas[n - 1].m_pInterpMatrix )
    {
      swap(pdOldOldAPrices, pdAPrices);
      swap(pdOldOldAPrices, pdOldAPrices);

      double dOldTimeWeight = m_pAdjointDatas[n].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if ( n + 1 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight;

      if ( dOldOldTimeWeight == 0 )
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdAPrices[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
      else
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdAPrices[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                          + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
    }
  }

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdPrices[nIdx] = m_pAdjointDatas[0].m_dOldTimeWeight * pdAPrices[nIdx];
}

void BackwardNumOutput::SensitivityByAdjointFD(BackwardInstData& instdata)
{  
  const size_t nNbSensitivities = instdata.m_pSensitivityData.size();
  
  if ( !nNbSensitivities )
    return;

  m_pdSensitivities.resize(nNbSensitivities);
  for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
    m_pdSensitivities[nIdxD] = 0.;

  const MorseMatrix* pSparseMatrix = instdata.m_pSparseMatrix;

  Array<double> pdAPrices(m_nNbX);
  Array<double> pdOldAPrices(m_nNbX);
  Array<double> pdOldOldAPrices(m_nNbX);
  
  Array<double> pdRHS(m_nNbX);
  Array<double> pdTmp(m_nNbX);
  Array<double> pdS2(m_nNbS);
  for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    pdS2[nIdxS] = m_pdS[nIdxS] * m_pdS[nIdxS];

  TriSparseLinearSolverGmres solver(m_nNbX);
  solver.EnableTransposeMatrixSolving();

  // setup the initial guess
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdAPrices[nIdx] = pdOldAPrices[nIdx] = pdOldOldAPrices[nIdx] = 0.;

  size_t nNbSystems = m_pAdjointDatas.size();

  // UpdateMe updates for initial values which is not associated to a system
  for (size_t n = nNbSystems - 1; n < nNbSystems; n--)
  {    
    swap(pdOldOldAPrices, pdAPrices);
    swap(pdOldOldAPrices, pdOldAPrices);

    const TridiagonalMatrix* pMatrix = m_pAdjointDatas[n].m_pMatrix.get();

    if ( n == nNbSystems - 1 )
    {
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        pdRHS[nIdx] = 0.;
    }
    else
    {
      double dOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if (  n + 2 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 2].m_dOldOldTimeWeight;

      if ( n + 1 < nNbSystems && m_pAdjointDatas[n + 1].m_pInterpMatrix )
        dOldOldTimeWeight = 0.;

      const InterpolationMatrix* 
        pInterpMatrix = m_pAdjointDatas[n].m_pInterpMatrix.get();

      if ( pInterpMatrix )
      {
        pInterpMatrix->ProductTransposeMatrixVector
                       ( pdOldAPrices.Get(), pdRHS.Get() );

        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight; 

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdRHS[nIdx] += pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
      }
      else
      {
        if ( dOldOldTimeWeight == 0 )
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdRHS[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
        else
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdRHS[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                        + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
      }

      instdata.ApplyBoundaryConditionToSensitivityRHS( pdRHS.Get() );
    }

    if ( m_pAdjointDatas[n].m_pRMatrix )
      m_pAdjointDatas[n].m_pRMatrix->ProductTransposeMatrixVector
          (m_pAdjointDatas[n].m_pdResiduals.Get(), pdRHS.Get(), true);

    solver.Solve( *pMatrix, *pSparseMatrix, pdRHS.Get(), pdAPrices.Get() );
    
    double* pdPrices = m_pAdjointDatas[n].m_pdPrices.Get();
   
    // Code duplicated from MakeSensitivityCoefficient. Should probably
    // find a better way.
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      size_t nOffset = nIdxR * m_nNbS;

      double* pdPricesTmp = pdPrices + nOffset;
      double* pdDeltasTmp = instdata.m_pdDeltas.Get() + nOffset;
      double* pdGammasTmp = instdata.m_pdGammas.Get() + nOffset;
      int* piFDTmp = m_pAdjointDatas[n].m_piFDFlags.Get() + nOffset;

      ComputeGammaFD(m_pdS, pdPricesTmp, m_nNbS, pdGammasTmp);

      if ( instdata.GetBoundaryCondition().GetLeftType() == BCType_Gamma )
        pdGammasTmp[0] = 0.;
      if ( instdata.GetBoundaryCondition().GetRightType() == BCType_Gamma )
        pdGammasTmp[m_nNbS - 1] = 0.;

      ComputeDelta(m_pdS, piFDTmp, pdPricesTmp, m_nNbS, pdDeltasTmp);
    }

    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
    {
      SensitivityData& data = instdata.m_pSensitivityData[nIdxD];

      const size_t nIdxR = data.m_nRegime1;
      const size_t nOffset = nIdxR * m_nNbS;
      
      double* pdGammasTmp = instdata.m_pdGammas.Get() + nOffset;
      double* pdPricesTmp = pdPrices + nOffset;
      double* pdDeltasTmp = instdata.m_pdDeltas.Get() + nOffset;
      double* pdTmpTmp = pdTmp.Get() + nOffset;
      double* pdAPricesTmp = pdAPrices.Get() + nOffset;

      // Determine what type of sensitivity is being computed.
      switch (data.m_sensitivityType)
      {
      case SensitivityType_Volatility:
        {
        // volatility coefficient
        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          pdTmpTmp[nIdxS] = instdata.GetModel().GetVolatilities()[nIdxR]
                          * pdGammasTmp[nIdxS] * pdS2[nIdxS];

        instdata.ApplyBoundaryConditionToSensitivityRHS( pdTmp.Get() );

        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          m_pdSensitivities[nIdxD] += pdAPricesTmp[nIdxS] * pdTmpTmp[nIdxS];

        break;
        } 
      case SensitivityType_DefaultIntensity:
      {
        // default intensity coefficient
        double dRecoveryValue = m_pAdjointDatas[n].m_dRecoveryValue;

        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          pdTmpTmp[nIdxS] = dRecoveryValue - pdPricesTmp[nIdxS]
                          + m_pdS[nIdxS] * pdDeltasTmp[nIdxS];

        instdata.ApplyBoundaryConditionToSensitivityRHS( pdTmp.Get() );

        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          m_pdSensitivities[nIdxD] += pdAPricesTmp[nIdxS] * pdTmpTmp[nIdxS];

        break;
      }
      case SensitivityType_JumpIntensity:
        {
        // intensity
        data.m_pSparseMatrix->ProductMatrixVector( pdPrices, pdTmp.Get() );

        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          pdTmpTmp[nIdxS] += - pdPricesTmp[nIdxS]
                          - data.m_dAmplitude * m_pdS[nIdxS] * pdDeltasTmp[nIdxS];

        instdata.ApplyBoundaryConditionToSensitivityRHS( pdTmp.Get() );

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          m_pdSensitivities[nIdxD] += pdAPrices[nIdx] * pdTmp[nIdx];

        break;
        }
      case SensitivityType_JumpAmplitude:
        {
        // amplitude
        data.m_pSparseMatrix->ProductMatrixVector( pdPrices, pdTmp.Get() );

        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          pdTmpTmp[nIdxS] -= data.m_dIntensity * m_pdS[nIdxS] * pdDeltasTmp[nIdxS];

        instdata.ApplyBoundaryConditionToSensitivityRHS( pdTmp.Get() );

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          m_pdSensitivities[nIdxD] += pdAPrices[nIdx] * pdTmp[nIdx];

        break;
        }
      } // end the switch
    }

    // Dividend caused interpolation
    if ( n > 0 && m_pAdjointDatas[n - 1].m_pInterpMatrix )
    {
      swap(pdOldOldAPrices, pdAPrices);
      swap(pdOldOldAPrices, pdOldAPrices);

      double dOldTimeWeight = m_pAdjointDatas[n].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if ( n + 1 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight;

      if ( dOldOldTimeWeight == 0 )
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdAPrices[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
      else
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdAPrices[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                          + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
    }
  }

  if ( !m_bSensitivityOnObjectif )
    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
      m_pdSensitivities[nIdxD] /= m_dPrice;
}

void BackwardNumOutput::SolveTransposeSystemFE
     (BackwardInstData& instdata, const double* pdPrices0, double* pdPrices)
{
  const TridiagonalMatrix* pMassMatrix = instdata.m_pMassMatrix;
  const MorseMatrix* pSparseMatrix = instdata.m_pSparseMatrix;
  
  Array<double> pdAPrices(m_nNbX);
  Array<double> pdOldAPrices(m_nNbX);
  Array<double> pdOldOldAPrices(m_nNbX);
  
  Array<double> pdRHS(m_nNbX);
  Array<double> pdTmp(m_nNbX);
  Array<double> pdS2(m_nNbS);
  for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    pdS2[nIdxS] = m_pdS[nIdxS] * m_pdS[nIdxS];

  TriSparseLinearSolverGmres solver(m_nNbX);
  solver.EnableTransposeMatrixSolving();

  // setup the initial guess
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdAPrices[nIdx] = pdOldAPrices[nIdx] = pdOldOldAPrices[nIdx] = pdPrices0[nIdx];

  size_t nNbSystems = m_pAdjointDatas.size();

  // UpdateMe updates for initial values which is not associated to a system
  for (size_t n = nNbSystems - 1; n < nNbSystems; n--)
  {    
    swap(pdOldOldAPrices, pdAPrices);
    swap(pdOldOldAPrices, pdOldAPrices);

    const TridiagonalMatrix* pMatrix = m_pAdjointDatas[n].m_pMatrix.get();

    if ( n == nNbSystems - 1 )
      pMassMatrix->ProductTransposeMatrixVector( pdAPrices.Get(), pdRHS.Get() );
    else
    {
      double dOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if (  n + 2 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 2].m_dOldOldTimeWeight;

      if ( n + 1 < nNbSystems && m_pAdjointDatas[n + 1].m_pInterpMatrix )
        dOldOldTimeWeight = 0.;

      const InterpolationMatrix* 
        pInterpMatrix = m_pAdjointDatas[n].m_pInterpMatrix.get();

      if ( pInterpMatrix )
      {
        pInterpMatrix->ProductTransposeMatrixVector
                       ( pdOldAPrices.Get(), pdRHS.Get() );

        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight; 

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdTmp[nIdx] = pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
       
        pMassMatrix->ProductTransposeMatrixVector
                     ( pdTmp.Get(), pdRHS.Get(), true );
      }
      else
      {
        if ( dOldOldTimeWeight == 0 )
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdTmp[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
        else
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdTmp[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                        + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;

        pMassMatrix->ProductTransposeMatrixVector( pdTmp.Get(), pdRHS.Get() ); 
      }
    }

    if ( m_pAdjointDatas[n].m_pRMatrix )
      m_pAdjointDatas[n].m_pRMatrix->ProductTransposeMatrixVector
          (m_pAdjointDatas[n].m_pdResiduals.Get(), pdRHS.Get(), true);

    solver.Solve( *pMatrix, *pSparseMatrix, pdRHS.Get(), pdAPrices.Get() );

    // Dividend caused interpolation
    if ( n > 0 && m_pAdjointDatas[n - 1].m_pInterpMatrix )
    {
      swap(pdOldOldAPrices, pdAPrices);
      swap(pdOldOldAPrices, pdOldAPrices);

      double dOldTimeWeight = m_pAdjointDatas[n].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if ( n + 1 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight;

      if ( dOldOldTimeWeight == 0 )
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdAPrices[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
      else
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdAPrices[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                          + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
    }
  }

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdPrices[nIdx] = m_pAdjointDatas[0].m_dOldTimeWeight * pdAPrices[nIdx];
}

void BackwardNumOutput::SensitivityByAdjointFE(BackwardInstData& instdata)
{  
  const size_t nNbSensitivities = instdata.m_pSensitivityData.size();
  
  if ( !nNbSensitivities )
    return;

  m_pdSensitivities.resize(nNbSensitivities);
  for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
    m_pdSensitivities[nIdxD] = 0.;

  // Get the original pointers
  const TridiagonalMatrix* pMassMatrix = instdata.m_pMassMatrix;
  const MorseMatrix* pSparseMatrix = instdata.m_pSparseMatrix;

  // memory allocation
  Array<double> pdAPrices(m_nNbX);
  Array<double> pdOldAPrices(m_nNbX);
  Array<double> pdOldOldAPrices(m_nNbX);
  
  Array<double> pdRHS(m_nNbX);
  Array<double> pdTmp(m_nNbX);
 
  TriSparseLinearSolverGmres solver(m_nNbX);
  solver.EnableTransposeMatrixSolving();

  // setup the initial guess
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    pdAPrices[nIdx] = pdOldAPrices[nIdx] = pdOldOldAPrices[nIdx] = 0.;

  const size_t nNbSystems = m_pAdjointDatas.size();
  
  for (size_t n = nNbSystems - 1; n < nNbSystems; n--)
  {  
    swap(pdOldOldAPrices, pdAPrices);
    swap(pdOldOldAPrices, pdOldAPrices);

    // Pointers that depends on each step
    const TridiagonalMatrix* pMatrix = m_pAdjointDatas[n].m_pMatrix.get();
    const double* pdPrices = m_pAdjointDatas[n].m_pdPrices.Get();

    if ( n == nNbSystems - 1 )
    {
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        pdRHS[nIdx] = 0.;
    }
    else
    {
      double dOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if ( n + 2 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 2].m_dOldOldTimeWeight;

      if ( n + 1 < nNbSystems && m_pAdjointDatas[n + 1].m_pInterpMatrix )
        dOldOldTimeWeight = 0.;

      const InterpolationMatrix* 
        pInterpMatrix = m_pAdjointDatas[n].m_pInterpMatrix.get();

      if ( pInterpMatrix )
      {
        pInterpMatrix->ProductTransposeMatrixVector
                       ( pdOldAPrices.Get(), pdRHS.Get() );

        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight; 

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdTmp[nIdx] = pdOldOldAPrices[nIdx] * dOldOldTimeWeight;
       
        pMassMatrix->ProductTransposeMatrixVector
                     ( pdTmp.Get(), pdRHS.Get(), true );
      }
      else
      {
        if ( dOldOldTimeWeight == 0 )
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdTmp[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
        else
          for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
            pdTmp[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                        + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;

        pMassMatrix->ProductTransposeMatrixVector( pdTmp.Get(), pdRHS.Get() ); 
      }
    }

    if ( m_pAdjointDatas[n].m_pRMatrix )
      m_pAdjointDatas[n].m_pRMatrix->ProductTransposeMatrixVector
          (m_pAdjointDatas[n].m_pdResiduals.Get(), pdRHS.Get(), true);

    solver.Solve( *pMatrix, *pSparseMatrix, pdRHS.Get(), pdAPrices.Get() );

    double dConst = m_pAdjointDatas[n].m_dRecoveryValue;

    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
    {
      SensitivityData& data = instdata.m_pSensitivityData[nIdxD];
      data.m_pSparseMatrix->ProductMatrixVector( pdPrices, pdTmp.Get() );

      if ( data.m_sensitivityType == SensitivityType_DefaultIntensity )
      {
        Array<double> pdBP(m_nNbX);

        for (size_t nIdxR = 0; nIdxR < instdata.m_nNbRegimes; nIdxR++)
        {
          double* pdBPTmp = pdBP.Get() + nIdxR * m_nNbS;
          if (nIdxR == data.m_nRegime1)
            for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
              pdBPTmp[nIdxS] = dConst;
          else
            for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
              pdBPTmp[nIdxS] = 0.;        
        }

        pMassMatrix->ProductMatrixVector( pdBP.Get(), pdTmp.Get(), true );
      }

      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        m_pdSensitivities[nIdxD] += pdAPrices[nIdx] * pdTmp[nIdx];
    }  
    
    // Dividend caused interpolation
    if ( n > 0 && m_pAdjointDatas[n - 1].m_pInterpMatrix )
    {
      swap(pdOldOldAPrices, pdAPrices);
      swap(pdOldOldAPrices, pdOldAPrices);

      double dOldTimeWeight = m_pAdjointDatas[n].m_dOldTimeWeight;

      double dOldOldTimeWeight = 0;    
      if ( n + 1 < nNbSystems )
        dOldOldTimeWeight = m_pAdjointDatas[n + 1].m_dOldOldTimeWeight;

      if ( dOldOldTimeWeight == 0 )
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdTmp[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight;
      else
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          pdTmp[nIdx] = pdOldAPrices[nIdx] * dOldTimeWeight
                          + pdOldOldAPrices[nIdx] * dOldOldTimeWeight;

      pMassMatrix->ProductTransposeMatrixVector
                   ( pdTmp.Get(), pdAPrices.Get() ); 
    }
  }

  if ( !m_bSensitivityOnObjectif )
    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
      m_pdSensitivities[nIdxD] /= m_dPrice;
}

} // namespace hg

} // namespace ito33
