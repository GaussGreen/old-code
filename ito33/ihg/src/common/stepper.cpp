/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/stepper.cpp
// Purpose:     base stepper class for ihg (finite differences)
// Author:      Zhang
// Created:     2004/01/21
// RCS-ID:      $Id: stepper.cpp,v 1.14 2006/06/13 15:45:11 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/boundary1d.h"

#include "ihg/stepper.h"
  
extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::numeric;


namespace ito33
{

namespace ihg
{

Stepper::Stepper(InstData& instdata, const finance::ComputationalFlags& flags)
               : StepperTriDiag(flags), m_instdata(instdata)
{
}

void Stepper::Alloc(size_t nNb)
{  
  Stepper::BaseClass::Alloc(nNb);

  // Compute the grid spacings for later use
  m_pdDeltaX = Array<double>(nNb);
  m_pdAreaX = Array<double>(nNb);
  
  m_pdInverseDeltaX = Array<double>(nNb);
  m_pdInverseAreaX = Array<double>(nNb);

  // Coefficient arrays of the particular PDE being solved
  m_pdCoe2nd = Array<double>(nNb);
  m_pdCoe1st = Array<double>(nNb);
  m_pdCoeZero = Array<double>(nNb);
  m_pdCoeConst = Array<double>(nNb);

  // Discretization arrays, constant for the price and Greek PDEs
  m_pdAlpha = Array<double>(nNb);
  m_pdBeta = Array<double>(nNb);

  // For Crank-Nicolson, we save the old values
  m_pdAlphaOld = Array<double>(nNb);
  m_pdBetaOld = Array<double>(nNb);
  m_pdCoeZeroOld = Array<double>(nNb);
  m_pdCoeConstOld = Array<double>(nNb);

  // Also need to save some vega values for Crank-Nicolson
  m_pdVegaCoeConstOld = Array<double>(nNb);

}

void Stepper::CalculateAreaArrays(const double *pdX, size_t nNbX)
{
  // The mesh and its size must be initialized by the derived class, since
  // this class doesn't always know what type of grid is being used
  size_t nIdxS;

  m_pdDeltaX[0] = pdX[1] - pdX[0];
  m_pdAreaX[0] = m_pdDeltaX[0];
  for (nIdxS = 1; nIdxS < nNbX - 1; nIdxS++)
  {
    m_pdDeltaX[nIdxS] = pdX[nIdxS + 1] - pdX[nIdxS];
    m_pdAreaX[nIdxS] = (pdX[nIdxS + 1] - pdX[nIdxS-1])*0.5;
  }
  m_pdAreaX[nNbX - 1] = m_pdDeltaX[nNbX-2];

  for(nIdxS = 0; nIdxS < nNbX - 1; nIdxS++)
  {
    m_pdInverseAreaX[nIdxS] = 1. / m_pdAreaX[nIdxS];
    m_pdInverseDeltaX[nIdxS] = 1. / m_pdDeltaX[nIdxS];
  }
  m_pdInverseAreaX[nIdxS] = 1. / m_pdAreaX[nIdxS];

  m_dLBC_CoefLeft = m_pdInverseDeltaX[0];
  m_dLBC_CoefRight = m_pdInverseDeltaX[nNbX - 2];
}


void Stepper::ModifyLinearBoundaryConditionCoef(const double *pdX, size_t nNbX)
{
  m_dLBC_CoefLeft = 1. / (pdX[1] - pdX[0]);
  m_dLBC_CoefRight = 1. / (pdX[nNbX - 1] - pdX[nNbX - 2]);
}


// Build the the right hand side
void Stepper::BuildRHS(const double* pdPrice1, const double* pdPriceM)
{
  switch(m_instdata.m_schemeType)
  {
  case ito33::numeric::SchemeType_Implicit:
    {
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdRHS[nIdx] = pdPrice1[nIdx] * m_instdata.m_dOldTimeWeight
                      + m_pdCoeConst[nIdx];
  
    break;
    }
  case ito33::numeric::SchemeType_ThreeLevel:
    {
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdRHS[nIdx] = pdPrice1[nIdx] * m_instdata.m_dOldTimeWeight +
                      pdPriceM[nIdx] * m_instdata.m_dOldOldTimeWeight +
                      m_pdCoeConst[nIdx];
    break;
    }
  case ito33::numeric::SchemeType_CrankNicolson:
    {    
    double dOneMinusTheta = 0.5;
    double dTheta = 0.5;

    m_pdRHS[0] = pdPrice1[0] * m_instdata.m_dOldTimeWeight 
      + dTheta * ( pdPrice1[0] * (-m_pdBetaOld[0] - m_pdCoeZeroOld[0])
      + pdPrice1[1] * m_pdBetaOld[0] + m_pdCoeConstOld[0])
      + dOneMinusTheta * m_pdCoeConst[0];
    for (size_t nIdx = 1; nIdx < m_nNbX-1; nIdx++)
    {
      double dTmp1 = -pdPrice1[nIdx] * (m_pdAlphaOld[nIdx] 
                  + m_pdBetaOld[nIdx] + m_pdCoeZeroOld[nIdx]);
      double dTmp2 = pdPrice1[nIdx-1] * m_pdAlphaOld[nIdx];
      double dTmp3 = pdPrice1[nIdx+1] * m_pdBetaOld[nIdx];
      m_pdRHS[nIdx] = pdPrice1[nIdx] * m_instdata.m_dOldTimeWeight 
                    + dTheta * (dTmp1 + dTmp2 + dTmp3 + m_pdCoeConstOld[nIdx])
                    + dOneMinusTheta * m_pdCoeConst[nIdx];
    }
    m_pdRHS[m_nNbX-1] = pdPrice1[m_nNbX-1] * m_instdata.m_dOldTimeWeight 
      + dTheta * ( pdPrice1[m_nNbX-1] * 
                    (-m_pdAlphaOld[m_nNbX-1] - m_pdCoeZeroOld[m_nNbX-1])
                  + pdPrice1[m_nNbX-2] * m_pdAlphaOld[m_nNbX-1]
                  + m_pdCoeConstOld[m_nNbX-1])
      + dOneMinusTheta * m_pdCoeConst[m_nNbX-1]; 
    break;
    }
  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, 
      "Stepper doesn't handle the current scheme type yet!");
  }
  
  // handle boundary condition
  if (m_instdata.GetBoundaryCondition().GetLeftType() == BCType_Dirichlet)
    m_pdRHS[0] = m_instdata.GetBoundaryCondition().GetLeftValue();
  if (m_instdata.GetBoundaryCondition().GetRightType() == BCType_Dirichlet)
    m_pdRHS[m_nNbX-1] = m_instdata.GetBoundaryCondition().GetRightValue();
}


void Stepper::BuildMatrix()
{
  const numeric::Boundary1D 
    &boundaryCondition = m_instdata.GetBoundaryCondition();
  
  double
    dOneMinusTheta = 1.0;

  if(m_instdata.m_schemeType == ito33::numeric::SchemeType_CrankNicolson)
    dOneMinusTheta = 0.5;

  // Get the diagonals of the tridiagonal matrix for easy access
  double* pdA = m_tridiagonalMatrix.GetA();
  double* pdB = m_tridiagonalMatrix.GetB();
  double* pdC = m_tridiagonalMatrix.GetC();

  // Construct the tridiagonal matrix
  // Handle left boundary
  if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
  {
    pdA[0] = 0.0;
    pdB[0] = 1.0;
    pdC[0] = 0.0;
  }
  else if (boundaryCondition.GetLeftType() == BCType_Gamma)
  {
    pdA[0] = 0.0;
    pdB[0] = m_instdata.m_dTimeWeight + dOneMinusTheta*(m_pdBeta[0] + m_pdCoeZero[0]);
    pdC[0] = -dOneMinusTheta * m_pdBeta[0];
  }
  else
  {
    ASSERT_MSG(false, "Unknown boundary type");
  }

  // handle all interior points. Alpha and beta should already be constructed.
  // m_pdCoeZero entries change. The matrix should be an M-matrix
  size_t nIdx;
  for (nIdx=1; nIdx < m_nNbX-1; nIdx++)
  {
    pdA[nIdx] = -dOneMinusTheta * m_pdAlpha[nIdx];
    pdB[nIdx] = m_instdata.m_dTimeWeight + 
                dOneMinusTheta
                  *(m_pdAlpha[nIdx] + m_pdBeta[nIdx] + m_pdCoeZero[nIdx]);
    pdC[nIdx] = -dOneMinusTheta * m_pdBeta[nIdx];
  }

  // handle right boundary
  if (boundaryCondition.GetRightType() == BCType_Dirichlet)
  {
    pdA[m_nNbX-1] = 0.0;
    pdB[m_nNbX-1] = 1.0;
    pdC[m_nNbX-1] = 0.0;
  }
  else if (boundaryCondition.GetRightType() == BCType_Gamma)
  {
    pdA[m_nNbX-1] = -dOneMinusTheta * m_pdAlpha[m_nNbX-1];
    pdB[m_nNbX-1] = m_instdata.m_dTimeWeight + 
                    dOneMinusTheta * 
                        (m_pdAlpha[m_nNbX-1] + m_pdCoeZero[m_nNbX-1]);
    pdC[m_nNbX-1] = 0.0;
  }
  else
  {
    ASSERT_MSG(false, "Unknown boundary type");
  }
    
}


// Build the alpha and beta entries.  These stay the same for the price
// and Greek PDEs, so only need to be made once
void Stepper::BuildAlphaBeta()
{
  // handle all interior points
  size_t nIdx;
  for (nIdx=1; nIdx < m_nNbX-1; nIdx++)
  {
    // Construct coefficients using central weighting. Alpha is for U_{i-1},
    // while beta is for U_{i+1}
    
    // Try central weighting first, since it is 2nd order accurate
    m_pdAlpha[nIdx] = (m_pdCoe2nd[nIdx]*m_pdInverseDeltaX[nIdx-1] 
                    - m_pdCoe1st[nIdx]*0.5)*m_pdInverseAreaX[nIdx];
    m_pdBeta[nIdx]  = (m_pdCoe2nd[nIdx]*m_pdInverseDeltaX[nIdx]   
                    + m_pdCoe1st[nIdx]*0.5)*m_pdInverseAreaX[nIdx];

     // Switch to forward differencing, if necessary
    if (m_pdAlpha[nIdx] < 0.0)
    {
      // Upstream weighting can introduce oscillations in delta and gamma. 
      // Use forward differencing
      m_pdAlpha[nIdx] = (m_pdCoe2nd[nIdx] * m_pdInverseDeltaX[nIdx-1])
                          * m_pdInverseAreaX[nIdx];
      m_pdBeta[nIdx]  = (m_pdCoe2nd[nIdx] * m_pdInverseAreaX[nIdx] 
                      + m_pdCoe1st[nIdx]) * m_pdInverseDeltaX[nIdx];
    }

    // Switch to backward differencing, if necessary
    if (m_pdBeta[nIdx] < 0.0)
    {
      // Use backward differencing, instead of downstream weighting
      m_pdAlpha[nIdx] = (m_pdCoe2nd[nIdx] * m_pdInverseAreaX[nIdx] 
                      - m_pdCoe1st[nIdx]) * m_pdInverseDeltaX[nIdx-1];
      m_pdBeta[nIdx]  = (m_pdCoe2nd[nIdx] * m_pdInverseDeltaX[nIdx])
                        * m_pdInverseAreaX[nIdx];
    }

  } // loop over internal nodes

  // Assume linear boundary condition. The Crank-Nicolson method assumes
  // that these values are saved
  m_pdAlpha[0] = 0.0;
  m_pdBeta[0] = m_pdCoe1st[0] * m_dLBC_CoefLeft;
  
  m_pdAlpha[m_nNbX-1] = -m_pdCoe1st[m_nNbX-1] * m_dLBC_CoefRight;
  m_pdBeta[m_nNbX-1] = 0.0;
}


} // namespace ihg

} // namespace ito33
