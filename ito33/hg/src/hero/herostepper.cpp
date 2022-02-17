/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/herostepper.cpp
// Purpose:     implementation of HG stepper class for HERO
// Created:     2005/09/26
// RCS-ID:      $Id: herostepper.cpp,v 1.4 2006/04/17 17:58:55 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

//#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"

#include "hg/model.h"
#include "hg/herostepper.h"

namespace ito33
{

namespace hg
{


void HeroStepper::Init()
{  
  HeroStepperBase::Init();

  m_pLinearSolver = AutoPtr<numeric::TriSparseLinearSolver>
    //( new numeric::TriSparseLinearSolverFixedPoint(m_nNbX) );  
    ( new numeric::TriSparseLinearSolverGmres(m_nNbX) );

  // The PDE coefficients are not quite right for the HERO. The drift
  // term should be beta_k, which includes the Sharpe ratio times total vol.
  std::vector<double> pdTotalVols = m_pModel->ComputeTotalVolatilities();
  double dSR = m_pModel->GetSharpeRatio();

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    m_pdCoe1st0[nIdxR1] += dSR * pdTotalVols[nIdxR1];
}

void HeroStepper::MakeCoefficients()
{

  // Let the base class do most of the work.  More importantly, this
  // makes it easier to support FD and FE methods
  HeroStepperBase::MakeCoefficients();

  // Must undo the addition of the interest rate in the zero coeffficient,
  // and add in the RHS of the PDE in section 4.2 of HERO notes to the 
  // constant term
  double dRate = m_instdata.m_dRate;

  // Get the matrix and vectors at the current time required for
  // the hero calculation
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    m_pdCoeZero[nIdxR] -= dRate;
 
    size_t nOffset = nIdxR * m_nNbS;
    double* pdCoeConst = m_pdCoeConst.Get() + nOffset;
    for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
      pdCoeConst[nIdx] += m_instdata.m_pdHeroTerms[nIdx + nOffset];
  }
}

void HeroStepper::Run()
{
  // Build the main system for price
  BuildMainSystem();

  m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                         m_sparseMatrix,
                         m_pdRHS.Get(),
                         m_instdata.m_pdPrices.Get());

}


} // namespace hg

} // namespace ito33
