/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/varianceswapstepper.cpp
// Purpose:     implementation of HG stepper class for variance swaps
// Created:     2006/03/05
// RCS-ID:      $Id: varianceswapstepper.cpp,v 1.2 2006/04/10 11:59:55 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"

#include "hg/model.h"
#include "hg/varianceswapstepper.h"
#include "hg/sensitivitymethod.h"

namespace ito33
{

  using namespace numeric;

namespace hg
{

void VarianceSwapStepper::Init()
{  
  VarianceSwapStepperBase::Init();

  m_pLinearSolver = AutoPtr<numeric::TriSparseLinearSolver>
     ( new TriSparseLinearSolverFixedPoint(m_nNbX) );  
     //( new TriSparseLinearSolverGmres(m_nNbX) );

}

void VarianceSwapStepper::Run()
{
  // Build the main system for price
  BuildMainSystem();

  // Solve
  m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                         m_sparseMatrix,
                         m_pdRHS.Get(),
                         m_instdata.m_pdPrices.Get());

  /*
  // Calculate sensitivities if required. Note that gamma is needed.
  // Pricing data comes from the solution computed above
  size_t nNbSensitivity = m_instdata.m_pSensitivityData.size();
  if ( nNbSensitivity 
    && m_instdata.m_sensitivityMethod == SensitivityMethod_PDE)
  {
    ComputeHelperArrays();

    for (size_t nIdxD = 0; nIdxD < nNbSensitivity; nIdxD++)
    {
      MakeSensitivityCoefficients(m_instdata.m_pSensitivityData[nIdxD]);

      BuildRHS(m_instdata.m_ppdOldSensitivities[nIdxD].Get(), 
               m_instdata.m_ppdOldOldSensitivities[nIdxD].Get());

      m_pLinearSolver->Solve(
                      m_tridiagonalMatrix,
                      m_sparseMatrix,
                      m_pdRHS.Get(),
                      m_instdata.m_ppdSensitivities[nIdxD].Get());

    } // loop over the model params
  } // if computing sensitivities

  if ( m_instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
  {
    SetupSensivityByAdjointData();

    // constant used to modify the matrix according to the constraint flags
    const double dLarge = 1.e8;

    if ( m_instdata.m_pConstraints )
    {
      double* pdB = m_instdata.m_aData.m_pMatrix->GetB();

      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        if ( m_instdata.m_piFrozenFlags[nIdx] )
          pdB[nIdx] += dLarge;
    }

    if ( m_instdata.m_pConstraints )
    {
      Array<int> piConstraintFlags(m_nNbX);
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        piConstraintFlags[nIdx] = m_instdata.m_piFrozenFlags[nIdx];
      
      m_instdata.m_aData.m_piConstraintFlags = piConstraintFlags;
    }
  }
  */
}


} // namespace hg

} // namespace ito33
