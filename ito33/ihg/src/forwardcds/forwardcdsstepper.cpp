/////////////////////////////////////////////////////////////////////////////
// Name:        forwardcds/forwardcdsstepper.cpp
// Purpose:     implementation of cds stepper class using forward PDE 
// Author:      David
// RCS-ID:      $Id: forwardcdsstepper.cpp,v 1.14 2006/06/13 15:44:45 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
#include "ihg/forwardcdsstepper.h"

#include "ito33/error.h"
#include "ito33/useexception.h"

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/schemetype.h"

extern const ito33::Error ITO33_UNEXPECTED;

using ito33::numeric::ComputeGamma;
using ito33::numeric::ComputeGammaFD;
using ito33::ihg::ForwardCDSStepper;

CVecteurDouble ForwardCDSStepper::operator*(CVecteurDouble &X) const
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


void ForwardCDSStepper::Run()
{
  // Everything except the last line is solving the forward option problem

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

  CalculateCDSTerms();
}

void ForwardCDSStepper::Init()
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

void ForwardCDSStepper::MakeCoefficients()
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


void ForwardCDSStepper::ComputeIntegralTerm(double* pdPrices)
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


void ForwardCDSStepper::CalculateCDSTerms()
{

  // The 2nd derivatives of the forward call option prices can be used to get the
  // risky yield curves, which are then used to get the forward CDS prices.
  Array<double> pdGammas(m_nNbX);
  ComputeGammaFD(m_pdX, m_instdata.m_pdPrices.Get(), m_nNbX, 
                 pdGammas.Get());
 
  // reset the gammas values for the first two points so that round off error
  // on boundary will be ignored. This is important when the hazard rate is 
  // very large at the boundary. In general it's not needed for right boundary.
  pdGammas[0] = pdGammas[1] = 0.;
 
  // This integral is used twice, but only needs to be computed once
  double dSum = 0.0;
  double dA = m_instdata.m_pdHazardRates[0] * pdGammas[0];
  double dB;

  for (size_t nIdx = 1; nIdx < m_nNbX; nIdx++)
  {
    dB = m_instdata.m_pdHazardRates[nIdx] * pdGammas[nIdx];
    dSum += 0.5 * (dA + dB) * (m_pdX[nIdx] - m_pdX[nIdx - 1]);
    dA = dB;
  }

  // Recovery and accrued terms (simple ODEs, or integrals, depeding on
  // viewpoint)
  // Note that the spread term is 1 the instant before the spread date, 
  // then jumps to zero just after the spread date.  This discontinuity is
  // handled by making the accrued fraction 1 at the spread date, and using
  // implicit timestepping after each spread date
  double dRecoveryTerm;
  double dAccruedTerm;
  double dRecoveryRHS = m_instdata.m_dRecovery * dSum;
  double dAccruedRHS = m_instdata.m_dAccruedFraction * m_instdata.m_dSpread * dSum;

  switch(m_instdata.m_schemeType)
  {
  case ito33::numeric::SchemeType_Implicit:
    {
      dRecoveryTerm = m_instdata.m_dRecoveryTermOld
                    + m_instdata.m_dTimeStep * dRecoveryRHS;

      dAccruedTerm = m_instdata.m_dAccruedTermOld
                   + m_instdata.m_dTimeStep * dAccruedRHS;
    break;
    }
  case ito33::numeric::SchemeType_ThreeLevel:
  case ito33::numeric::SchemeType_CrankNicolson:
    {    
      dRecoveryTerm = m_instdata.m_dRecoveryTermOld + m_instdata.m_dTimeStep 
                    * 0.5 * (dRecoveryRHS + m_instdata.m_dRecoveryRHSOld);

      dAccruedTerm = m_instdata.m_dAccruedTermOld + m_instdata.m_dTimeStep 
                   * 0.5 * (dAccruedRHS + m_instdata.m_dAccruedRHSOld);
    break;
    }
  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, 
      "Unknown scheme type in forward cds stepper");
  }

  // Spared term (just a number).  Update the data if at new spread date
  if (m_instdata.m_dAccruedFraction == 1.0)
  {
    // The forward option code assume a linear boundary condition. Hence,
    // we can just use one-sided forward differencing to get the first 
    // derivative at the left boundary point (assumed to be zero or close
    // enough to zero to not make a difference).
    double d1stDerivAtZero = (m_instdata.m_pdPrices[1] - m_instdata.m_pdPrices[0])
                           / (m_pdX[1] - m_pdX[0]);

    // Subtract from the sum since the formula needs -d1stDerivAtZero
    m_instdata.m_dSpreadTermOld -= m_instdata.m_dSpread * d1stDerivAtZero;
  }
  double dSpreadTerm = m_instdata.m_dSpreadTermOld;

  // Combine to get the CDS price at this timestep
  m_instdata.m_dCDSPrice = dRecoveryTerm - dSpreadTerm - dAccruedTerm;
  m_instdata.m_dAccruedTerm = dAccruedTerm;
  m_instdata.m_dRecoveryTerm = dRecoveryTerm;
  m_instdata.m_dSpreadTerm = dSpreadTerm;

  // Update the data in instdata in case of multistep methods (well, just
  // modified Euler)
  m_instdata.m_dRecoveryRHSOld = dRecoveryRHS;
  m_instdata.m_dAccruedRHSOld = dAccruedRHS;

  m_instdata.m_dRecoveryTermOld = dRecoveryTerm;
  m_instdata.m_dAccruedTermOld = dAccruedTerm;
}

