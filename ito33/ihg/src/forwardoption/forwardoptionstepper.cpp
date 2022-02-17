/////////////////////////////////////////////////////////////////////////////
// Name:        forwardoption/forwardoptionstepper.cpp
// Purpose:     implementation of option stepper class using forward PDE 
// Author:      Wang
// RCS-ID:      $Id: forwardoptionstepper.cpp,v 1.11 2005/01/28 19:47:17 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/tridiagonalpreconditioner.h"
#include "ito33/numeric/gmressolver.h"
#include "ito33/numeric/schemetype.h"

#include "ihg/forwardoptionstepper.h"

using ito33::ihg::ForwardOptionStepper;

CVecteurDouble ForwardOptionStepper::operator*(CVecteurDouble &X) const
{
  // NOTE: This code is duplicated in the ComputeIntegralTerm function below

  CVecteurDouble F( m_nNbX );

  m_tridiagonalMatrix.ProductMatrixVector(&X[0], &F[0]);

  // Rewrite the forward equation, try to get better integration term handle
  // testing showed that it is quicker to allocate this space each time instead
  // of using a member variable
  Array<double> pdGammas(m_nNbX);

  // Compute gamma using finite differences, as in the discretization
  numeric::ComputeGammaFD(m_pdX, &X[0], m_nNbX, pdGammas.Get());

  // Compute the integral  int_K^\infty C_{SS} * lambda dS, where C is the current 
  // solution, and lambda is the hazard rate
  double dA;
  double dB = pdGammas[m_nNbX - 1] * m_instdata.m_pdHazardRates[m_nNbX - 1];
  double dSum = 0.0;  
  size_t nIdx;

  // handle Crank-Nicolson timestepping. For efficiency, keep the if
  // statement outside of the for loop
  if ( m_instdata.m_schemeType == ito33::numeric::SchemeType_CrankNicolson)
  {
    for (nIdx = m_nNbX - 2; nIdx < m_nNbX - 1; nIdx--)
    {
      dA = pdGammas[nIdx] * m_instdata.m_pdHazardRates[nIdx];
      dSum += 0.5 * (dA + dB) * (m_pdX[nIdx+1] - m_pdX[nIdx]);
      F[nIdx] -= 0.5 * dSum * m_pdX[nIdx];
      dB = dA;
    }
  }
  else
  {
    for (nIdx = m_nNbX - 2; nIdx < m_nNbX - 1; nIdx--)
    {
      dA = pdGammas[nIdx] * m_instdata.m_pdHazardRates[nIdx];
      dSum += 0.5 * (dA + dB) * (m_pdX[nIdx+1] - m_pdX[nIdx]);
      F[nIdx] -= dSum * m_pdX[nIdx];
      dB = dA;
    }
  }

  return F;
}


void ForwardOptionStepper::Run()
{

  // Make the PDE coefficients
  MakeCoefficients();

  // Make the alpha/beta coefficients
  BuildAlphaBeta();
 
  // Build the tridiagonal matrix
  BuildMatrix();

  // Build the discrete system
  BuildRHS(m_instdata.m_pdOldPrices.Get(),
           m_instdata.m_pdOldOldPrices.Get());
   
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_instdata.m_pdPrices[nIdx] = m_instdata.m_pdOldPrices[nIdx];

  // For an explicit solve of the integral term, use a linear solver

  // If the hazard rate is time only, then we can do a simple linear
  // solve
  if (m_instdata.m_bIsHRTimeOnly)
  {
    m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                           m_pdRHS.Get(),
                           m_instdata.m_pdPrices.Get());
  }
  else
  {
    // this is essentially an iterative solve
    m_pGmresSolver->Solve(*this, *m_pPreconditioner, 
                          m_pdRHS.Get(),
                          m_instdata.m_pdPrices.Get());

    // Save the integral term in case Crank-Nicolson timestepping is used
    ComputeIntegralTerm(m_instdata.m_pdPrices.Get());
  }

  // Now that solving is done, save values in case Crank-Nicolson is used
  swap(m_pdCoeConst, m_pdCoeConstOld);
  swap(m_pdAlpha, m_pdAlphaOld);
  swap(m_pdBeta, m_pdBetaOld);
  swap(m_pdCoeZero, m_pdCoeZeroOld);
}


void ForwardOptionStepper::Init()
{  
  // Get the grid size 
  m_nNbX = m_instdata.m_nNbS;

  // allocate storage
  Alloc(m_nNbX);
   
  // Storage for the integral terms, as computed by ComputeIntegralTerm
  //m_pdIntegralValues = Array<double>(m_nNbX);

  CalculateAreaArrays(m_instdata.m_pdS, m_nNbX);
  
  m_pdX = m_instdata.m_pdS;

  size_t nRestart = m_nNbX;
  if (nRestart > 200)
    nRestart = 200;

  m_pGmresSolver = AutoPtr<numeric::GmresSolver>
                   (
                     new numeric::GmresSolver(m_nNbX, nRestart, 400, 1.e-14)
                   );

  m_pPreconditioner = AutoPtr<numeric::TridiagonalPreconditioner>
                      (
                        new numeric::TridiagonalPreconditioner
                            (
                              m_tridiagonalMatrix
                            )
                      );

  m_pLinearSolver = AutoPtr<numeric::TridiagonalSolver>
                    (
                      new numeric::TridiagonalSolver(m_nNbX)
                    );

}

void ForwardOptionStepper::MakeCoefficients()
{
  size_t nIdx;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  // Recovery value is zero, no need to take into account

  // Setup the coefficients
  if (m_instdata.m_bIsHRTimeOnly == true)
  {
    // For time only hazard rate, use the form of equations where the
    // integral uses the second derivative of the hazard rate
    // In this case, the integral is zero, since the 2nd deriv of the hazard
    // rates is zero. Further, the 1st deriv of the hazard rate is one.
    double* pdHazardRates = m_instdata.m_pdHazardRates.Get();
    for (nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = 0.5 * m_instdata.m_pdVolsSquared[nIdx];
      m_pdCoe1st[nIdx] = - (dRate - dForeignRate + pdHazardRates[nIdx]);
      m_pdCoeZero[nIdx] = dForeignRate;
      m_pdCoeConst[nIdx] = 0.;
    }
  }
  else
  {
    // Use the form of equations where the intergal term uses the 2nd 
    // deriv of the solution (integrate the usual equations by parts
    // twice).  However, these equations will be solved by GMRES, so
    // set the integral term (coeConst) to zero. It is computed in the
    // matrix-vector multiply.
    for (nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = 0.5 * m_instdata.m_pdVolsSquared[nIdx];
      m_pdCoe1st[nIdx] = - (dRate - dForeignRate);
      m_pdCoeZero[nIdx] = dForeignRate;
      m_pdCoeConst[nIdx] = 0.;
    }
  }

  for (nIdx = 0; nIdx < m_nNbX; nIdx++)
  {
    m_pdCoe2nd[nIdx] *= m_instdata.m_pdS[nIdx] * m_instdata.m_pdS[nIdx];
    m_pdCoe1st[nIdx] *= m_instdata.m_pdS[nIdx];
  }
}


void ForwardOptionStepper::ComputeIntegralTerm(double* pdPrices)
{

  // NOTE: This code is duplicated in the matrix-vector function above

  size_t nIdx;
  
  // Initial testing showed that it is quicker to allocate memory here
  // instead of as a member variable (caching issue?)
  Array<double> pdGammas(m_nNbX);

  // Section used for two integration by parts (full cancellation) 
  // Rewrite the forward equation, try to get better integration term handle
  numeric::ComputeGammaFD(m_instdata.m_pdS, pdPrices, m_instdata.m_nNbS,
                          pdGammas.Get());

  // Integral must go to zero as K->infinity.  Otherwise unstable.
  //m_pdIntegralValues[m_nNbX - 1] = 0.0;
  m_pdCoeConst[m_nNbX - 1] = 0.0;

  // Compute the integral  int_K^\infty C_{SS} * lambda dS, where C is the current 
  // solution, and lambda is the hazard rate
  double dA;
  double dB = pdGammas[m_nNbX-1] * m_instdata.m_pdHazardRates[m_nNbX - 1];
  double dSum = 0.0;  
  for (nIdx = m_nNbX - 2; nIdx < m_nNbX - 1; nIdx--)
  {
    dA = pdGammas[nIdx] * m_instdata.m_pdHazardRates[nIdx];
    dSum += 0.5 * (dA + dB) * (m_pdX[nIdx+1] - m_pdX[nIdx]);
    //m_pdIntegralValues[nIdx] = dSum * m_pdX[nIdx];
    m_pdCoeConst[nIdx] = dSum * m_pdX[nIdx];
    dB = dA;
  }

}
