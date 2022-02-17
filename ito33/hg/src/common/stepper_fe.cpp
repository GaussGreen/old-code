/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/stepper_fe.cpp
// Purpose:     base stepper class for HG (finite element)
// Created:     2005/01/18
// RCS-ID:      $Id: stepper_fe.cpp,v 1.26 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/boundary1d.h"

#include "hg/model.h"
#include "hg/stepper_fe.h"

namespace ito33
{
  using namespace numeric;

namespace hg
{


inline double Sqr(double d) { return d * d; }

inline double Max(double dA, double dB)
{
  return dA > dB ? dA : dB;
}

inline double Min(double dA, double dB)
{
  return dA > dB ? dB : dA;
}

// Compute the integral \int_{dMin}^{dMax} (x - dX1) (x - dX2) dx
// Two points gauss formula gives here the exact result
double Integrate(double dMin, double dMax, double dX1, double dX2)
{
  double dTmp = ((dMax - dMin) * (-0.577350269189626) + dMin + dMax) * 0.5;
  double dTmp1 = ((dMax - dMin) * 0.577350269189626 + dMin + dMax) * 0.5;

  return   ( (dTmp - dX1) * (dTmp - dX2) + (dTmp1 - dX1) * (dTmp1 - dX2) )
         * 0.5 * (dMax - dMin);
}

// Compute the integral \int_{dMin}^{dMax} (dCoe exp(x) - dX1) (x - dX2) dx
// This is just an estimation
double IntegrateExp(double dMin, double dMax, 
                    double dCoe, double dX1, double dX2)
{
  double dTmp = ((dMax - dMin) * (-0.577350269189626) + dMin + dMax) * 0.5;
  double dTmp1 = ((dMax - dMin) * 0.577350269189626 + dMin + dMax) * 0.5;
  
  return   (  (dCoe * exp(dTmp) - dX1) * (dTmp - dX2)
            + (dCoe * exp(dTmp1) - dX1) * (dTmp1 - dX2) )
         * 0.5 * (dMax - dMin);
}

/*
// Compute the integral \int_{dMin}^{dMax} (dCoe exp(x) - dX1) (x - dX2) dx
// The exact formula, souldn't be really needed
double IntegrateExp(double dMin, double dMax, 
                    double dCoe, double dX1, double dX2)
{
  return - 0.5 * dX1 * (Sqr(dMax) - Sqr(dMin)) + dX1 * dX2 * (dMax - dMin)
         + dCoe * (  exp(dMax) * (dMax - dX2 - 1.) 
                   - exp(dMin) * (dMin - dX2 - 1.) );
}
*/

void StepperFE::Init()
{
  const Model& model = m_instdata.GetModel();

  // Store the pointer to model for quick access
  // We can't use const reference since it would then need to be done in ctor
  m_pModel = &model;

  // store some vectors for quick access
  m_nNbRegimes = model.GetNbRegimes();

  m_pdVols = model.GetVolatilities();

  m_pdDefaultIntensities = model.GetJumpsToDefault();

  // Time independent coefficient arrays of the particular PDE being solved
  m_pdCoe2nd0 = Array<double>(m_nNbRegimes);
  m_pdCoe1st0 = Array<double>(m_nNbRegimes);
  m_pdCoeZero0 = Array<double>(m_nNbRegimes);

  // Backward equation, forward equation has different coefficients
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    m_pdCoe2nd0[nIdxR1] = 0.5 * m_pdVols[nIdxR1] * m_pdVols[nIdxR1];

    m_pdCoe1st0[nIdxR1] = 0.;
    m_pdCoeZero0[nIdxR1] = 0.0;
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_pModel->GetJumps(nIdxR1, nIdxR2);

      Jumps::const_iterator pJump;
      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
      {
        m_pdCoe1st0[nIdxR1] += pJump->GetIntensity() * pJump->GetAmplitude();
        m_pdCoeZero0[nIdxR1] += pJump->GetIntensity();  
      }                                     
    }

    double dDefaultIntensity = m_pdDefaultIntensities[nIdxR1];
    m_pdCoe1st0[nIdxR1] += dDefaultIntensity * (- 1.);
    m_pdCoeZero0[nIdxR1] += dDefaultIntensity;
  }

  // Coefficent array of the PDE at each time step 
  m_pdCoe2nd = Array<double>(m_nNbRegimes);
  m_pdCoe1st = Array<double>(m_nNbRegimes);
  m_pdCoeZero = Array<double>(m_nNbRegimes);

  m_bSensitivityRHS = false;
}

void StepperFE::Alloc(size_t nNb)
{  
  pricing::StepperTriDiag::Alloc(nNb);

  // Allocate memory for mass matrix
  m_massMatrix.SetDimension(nNb);

  // Helper arrays related to the grid for discretization
  // We are assuming that regimes use the same grid
  m_pdDeltaX = Array<double>(m_nNbS - 1);
  m_pdInverseDeltaX = Array<double>(m_nNbS - 1);

  // Helper array for building right hand side
  m_pdCoeConst = Array<double>(nNb);

  // Helper array for building right hand side for sensitivity system
  // Will be added directly to the right hand side of sensitivity system,
  // unlike the previous array
  m_pdCoeConstSensitivity = Array<double>(nNb);
}

void StepperFE::BuildMassMatrix()
{  
  double *pdA, *pdB, *pdC;

  pdA = m_massMatrix.GetA();
  pdB = m_massMatrix.GetB();
  pdC = m_massMatrix.GetC();
  
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  { 
    // pdATmp[nIdx] = \int_E_{i - 1}(x) E_i(x) dx
    // pdBTmp[nIdx] = \int_E_i(x) E_i(x) dx
    // pdCTmp[nIdx] = \int_E_{i + 1}(x) E_i(x) dx
    
    size_t nOffset = nIdxR * m_nNbS;
    double* pdATmp = pdA + nOffset;
    double* pdBTmp = pdB + nOffset;
    double* pdCTmp = pdC + nOffset;

    size_t nIdx = 0;

    pdATmp[nIdx] = 0.;
    pdBTmp[nIdx] = 0.;
    pdCTmp[nIdx] = 0.; 
      
    for (nIdx = 1; nIdx < m_nNbS - 1; nIdx++)
    {
      pdATmp[nIdx] = m_pdDeltaX[nIdx - 1] / 6.;
      pdBTmp[nIdx] = (m_pdDeltaX[nIdx - 1] + m_pdDeltaX[nIdx]) / 3.;
      pdCTmp[nIdx] = m_pdDeltaX[nIdx] / 6.;
    }
    
    pdATmp[nIdx] = 0.;
    pdBTmp[nIdx] = 0.;
    pdCTmp[nIdx] = 0.;
  }

  // Values for gamma boundary condition
  // Will be used to linearly (on S) interpolate the boundary value
  const double* pdS = m_instdata.m_pdS;
  m_dLBC_CoefLeft = (pdS[0] - pdS[2]) / (pdS[1] - pdS[2]);
  m_dLBC_CoefRight = (pdS[m_nNbS - 1] - pdS[m_nNbS - 3])
                   / (pdS[m_nNbS - 2] - pdS[m_nNbS - 3]);
}   

void StepperFE::BuildMatrix()
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    double dCoe2nd = m_pdCoe2nd[nIdxR];
    double dCoe1st = m_pdCoe1st[nIdxR];
    double dCoeZero = m_pdCoeZero[nIdxR];

    double* pdA = m_tridiagonalMatrix.GetA() + nIdxR * m_nNbS;
    double* pdB = m_tridiagonalMatrix.GetB() + nIdxR * m_nNbS;
    double* pdC = m_tridiagonalMatrix.GetC() + nIdxR * m_nNbS;

    size_t nIdx = 0;
 
    // Modification due to boundary condition
    // At the boundary, pde is forced to be linear or constant in S
    pdA[nIdx] = 0.;
    pdB[nIdx] = 1.;

    if ( boundaryCondition.GetLeftType() == BCType_Gamma )
      pdC[nIdx] = - m_dLBC_CoefLeft;
    else if ( boundaryCondition.GetLeftType() == BCType_Dirichlet )
      pdC[nIdx] = 0.;    

    for (nIdx = 1; nIdx < m_nNbS - 1; nIdx++)
    {  
      double dDeltaX0 = m_pdDeltaX[nIdx - 1];  
      double dInvDeltaX0 = m_pdInverseDeltaX[nIdx - 1];
      double dDeltaX = m_pdDeltaX[nIdx];  
      double dInvDeltaX = m_pdInverseDeltaX[nIdx];   
      
      pdA[nIdx] = - dCoe2nd * dInvDeltaX0 - 0.5 * dCoe1st
                + dCoeZero * dDeltaX0 / 6.;      
      pdB[nIdx] = dCoe2nd * (dInvDeltaX0 + dInvDeltaX)
                + dCoeZero * (dDeltaX0 + dDeltaX) / 3.; 
      pdC[nIdx] = - dCoe2nd * dInvDeltaX + 0.5 * dCoe1st 
                + dCoeZero * dDeltaX / 6.;  
    }

    // Modification due to boundary condition
    // At the boundary, pde is forced to be linear or constant in S
    if ( boundaryCondition.GetRightType() == BCType_Gamma )
      pdA[nIdx] = - m_dLBC_CoefRight;
    else if ( boundaryCondition.GetRightType() == BCType_Dirichlet )
      pdA[nIdx] = 0.;

    pdB[nIdx] = 1.;
    pdC[nIdx] = 0.;
  }
}
    
// Build the the right hand side
void StepperFE::BuildRHS(const double* pdOldPrice, const double* pdOldOldPrice)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  // only implicit or three level will be implemented
  ASSERT(   m_instdata.m_schemeType == SchemeType_Implicit
         || m_instdata.m_schemeType == SchemeType_ThreeLevel );

  Array<double> pdTmp(m_nNbX);

  switch(m_instdata.m_schemeType)
  {
  case SchemeType_Implicit:
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      pdTmp[nIdx] = pdOldPrice[nIdx] * m_instdata.m_dOldTimeWeight
                  + m_pdCoeConst[nIdx];
    break;

  case SchemeType_ThreeLevel:
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      pdTmp[nIdx] = pdOldPrice[nIdx] * m_instdata.m_dOldTimeWeight
                  + pdOldOldPrice[nIdx] * m_instdata.m_dOldOldTimeWeight
                  + m_pdCoeConst[nIdx];
    break;
  }

  m_massMatrix.ProductMatrixVector( pdTmp.Get(), m_pdRHS.Get() );  
  
  // Add the contribution of sensitivity specific const coefficient
  if (m_bSensitivityRHS)
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdRHS[nIdx] += m_pdCoeConstSensitivity[nIdx];
    
  // Modification due to boundary condition
  // Can't be done before the product matrix vector
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    size_t nOffset = nIdxR * m_nNbS;
    double* pdRHS = m_pdRHS.Get() + nOffset;

    if (boundaryCondition.GetLeftType() == BCType_Dirichlet )
      pdRHS[0] = boundaryCondition.GetLeftValue();

    if ( boundaryCondition.GetRightType() == BCType_Dirichlet )
      pdRHS[m_nNbS - 1] = boundaryCondition.GetRightValue();
  }
}

// A row or stupidest implementation of the jump system. It's quite difficult
// to use it to compute the sensitivity jump system but it can serve as a 
// reference that other implementation can be checked against
void StepperFE::BuildJumpSystem()
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  // Used to find the columns of the sparse matrix. Lists are used since
  // we don't know how many nonzero entries appear in each row
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  // Pass 1. Get the structure of the sparse matrix
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& pJumps = m_pModel->GetJumps(nIdxR1, nIdxR2);
     
      // loop over all jumps from regime I to regime J
      Jumps::const_iterator pJump;
      for (pJump = pJumps.begin(); pJump != pJumps.end(); ++pJump)
      {
        const double dAmplitude = pJump->GetAmplitude();
        
        BuildSensitivityJumpSystemStructure
        ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

      } // loop over jumps from regime I to regime J
    } // inner regime loop
  } // outer regime loop

  // A weird implementation of the gamma boundary condition
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    if ( boundaryCondition.GetLeftType() == BCType_Gamma )
    {
      const size_t nIdxRow = nIdxR1 * m_nNbS; 
      ppnColumnLists[nIdxRow].push_back(nIdxRow + 2);
    }

    if ( boundaryCondition.GetRightType() == BCType_Gamma )
    {
      const size_t nIdxRow = nIdxR1 * m_nNbS + m_nNbS - 1; 
      ppnColumnLists[nIdxRow].push_back(nIdxRow - 2);
    }
  }

  shared_ptr<MorseStruct> 
    pMS( new MorseStruct( ppnColumnLists.Get(), m_nNbX ) );

  m_sparseMatrix.Init(pMS, 0.);   

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      // loop over all jumps from regime I to regime J
      const Jumps& pJumps = m_pModel->GetJumps(nIdxR1, nIdxR2);   
  
      Jumps::const_iterator pJump;
      for (pJump = pJumps.begin(); pJump != pJumps.end(); ++pJump)
      {
        const double dIntensity = pJump->GetIntensity();
        const double dAmplitude = pJump->GetAmplitude();

        BuildPartialSensitivityJumpSystem
        ( nIdxR1, nIdxR2, dAmplitude, - dIntensity, m_sparseMatrix);

      } // loop over jumps from regime I to regime J
    } // inner regime loop
  } // outer regime loop

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    if (boundaryCondition.GetLeftType() == BCType_Gamma)
    {
      const size_t nIdxRow = nIdxR1 * m_nNbS; 
      m_sparseMatrix(nIdxRow, nIdxRow + 2) = - (1. - m_dLBC_CoefLeft);
    }

    if (boundaryCondition.GetRightType() == BCType_Gamma )
    {
      const size_t nIdxRow = nIdxR1 * m_nNbS + m_nNbS - 1; 
      m_sparseMatrix(nIdxRow, nIdxRow - 2) = - (1. - m_dLBC_CoefRight);
    }
  }
}

void
StepperFE::BuildSensitivityJumpSystemStructure
           (size_t nIdxR1, size_t nIdxR2, double dAmplitude, 
            std::list<size_t>* ppnColumnLists)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  Array<double> pdNewX(m_nNbS);
  Array<int> piIV(m_nNbS);

  const int iNbIntervals = int(m_nNbS) - 1;

  // Pass 1. Get the structure of the sparse matrix
  const size_t nOffsetRow = nIdxR1 * m_nNbS;
  const size_t nOffsetColumn = nIdxR2 * m_nNbS;
      
  const double dSize = log(1. + dAmplitude);

  for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    pdNewX[nIdxS] = m_pdX[nIdxS] + dSize;

  SearchIntervals(m_pdX, m_nNbS, pdNewX.Get(), piIV.Get(), m_nNbS);

  // do only internal points, leave alone the boundary points
  for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
  {
    const size_t nIdxRow = nOffsetRow + nIdxS;
    
    // The support of the shifted base function \phi_{nIdx}
    const int iIV0 = piIV[nIdxS - 1];
    const int iIV2 = piIV[nIdxS + 1];

    // The indexs of base functions that have a non empty intersection 
    // with the above function
    size_t nIdxJMin = (iIV0 >= 0) ? iIV0 : 0;
    size_t nIdxJMax = (iIV2 <= iNbIntervals - 1) ? iIV2 + 1 : m_nNbS - 1;

    // If the left boundary is gamma condition, then at least 0,1 should
    // be counted for the columns
    if (   iIV0 == -1
        && boundaryCondition.GetLeftType() == BCType_Gamma )
    {
      nIdxJMax = (nIdxJMax > 1) ? nIdxJMax : 1;
    }

    // If the right boundary is gamma condition, then at least m_nNbS - 2
    // and m_nNbS - 1 should be counted for the columns
    if (   iIV2 == iNbIntervals
        && boundaryCondition.GetRightType() == BCType_Gamma)
    {
      nIdxJMin = (nIdxJMin > m_nNbS - 2) ? m_nNbS - 2 : nIdxJMin;
    }

    // just add those indexs to the column of the current row.
    for (size_t nIdxJ = nIdxJMin; nIdxJ <= nIdxJMax; nIdxJ++)
      ppnColumnLists[nIdxRow].push_back(nOffsetColumn + nIdxJ);
  
  } // loop over space nodes
}

AutoPtr<MorseMatrix>
StepperFE::BuildIntensitySensitivityJumpSystem
(size_t nIdxR1, size_t nIdxR2, double dAmplitude, double /* dIntensity */)
{
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  BuildSensitivityJumpSystemStructure
  ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

  shared_ptr<MorseStruct> pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );

  AutoPtr<MorseMatrix> pSparseMatrix(new MorseMatrix( pMS ) ); 

  pSparseMatrix->Init(pMS, 0);
  
  BuildPartialSensitivityJumpSystem
  ( nIdxR1, nIdxR2, dAmplitude, 1., *pSparseMatrix);

  return pSparseMatrix;
}

// A row or stupidest implementation of the jump system for sensitivity.
// It may be merged into the building of the main jump system
void
StepperFE::BuildPartialSensitivityJumpSystem
           (size_t nIdxR1, size_t nIdxR2,
            double dAmplitude, double dMultiplier,
            MorseMatrix& sparseMatrix)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  const size_t nOffsetRow = nIdxR1 * m_nNbS;
  const size_t nOffsetColumn = nIdxR2 * m_nNbS;

  const double dSize = log(1. + dAmplitude);

  const double dS0 = m_instdata.m_pdS[0], dS1 = m_instdata.m_pdS[1],
               dSN2 = m_instdata.m_pdS[m_nNbS - 2], 
               dSN1 = m_instdata.m_pdS[m_nNbS - 1],
               dInvDeltaS0 = 1. / (dS1 - dS0), 
               dInvDeltaSN2 = 1. / (dSN1 - dSN2); 

  double dMin, dMax, dTmp;
  size_t nIdxI, nIdxJ;
  double dNewXJm, dNewXJ, dNewXJp;
  double dXIm, dXI, dXIp;
  double dInvDeltaColumn;
  
  size_t nIdxColumn;

  const double dCoe = 1. + dAmplitude;

  double dInvDeltaRow;
       
  // The interval points for the row, we don't care about the boundary 
  // points
  for (nIdxI = 1; nIdxI < m_nNbS - 1; nIdxI++)
  {
    // The current row
    const size_t nIdxRow = nOffsetRow + nIdxI;

    // left support of the row base function
    // the base function is (dX - dXIm) / (dXI - dXIm)
    dXIm = m_pdX[nIdxI - 1];
    dXI = m_pdX[nIdxI];

    dInvDeltaRow = dMultiplier * m_pdInverseDeltaX[nIdxI - 1];
    
    // the left boundary
    nIdxJ = 0;
    nIdxColumn = nOffsetColumn + nIdxJ;

    // dNewXJm = - infinity;
    dNewXJ = m_pdX[nIdxJ] - dSize;

    dMin = dXIm; // Max(dNewXJm, dXIm);
    dMax = Min(dNewXJ, dXI);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
      {
        // base function is f_0 
        sparseMatrix(nIdxRow, nIdxColumn)
            += dInvDeltaRow
             * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
      }
      else if (boundaryCondition.GetLeftType() == BCType_Gamma)
      {
        // base function is f_0 * (dCoe * dS - dS1) / (dS0 - dS1)
        //                + f_1 * (dCoe * dS - dS0) / (dS1 - dS0)
        dTmp = dInvDeltaRow * dInvDeltaS0;

        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * IntegrateExp(dMin, dMax, dCoe, dS1, dXIm);

        sparseMatrix(nIdxRow, nIdxColumn + 1)
            += dTmp * IntegrateExp(dMin, dMax, dCoe, dS0, dXIm);  
      }
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    dNewXJp = m_pdX[nIdxJ + 1] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

    dMin = Max(dNewXJ, dXIm);
    dMax = Min(dNewXJp, dXI);

    if (dMin < dMax)
    {
      // base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
      dTmp = dInvDeltaRow * dInvDeltaColumn;
      sparseMatrix(nIdxRow, nIdxColumn)
          -= dTmp * Integrate(dMin, dMax, dNewXJp, dXIm);
    }

    // the right boundary
    nIdxJ = m_nNbS - 1;
    nIdxColumn = nOffsetColumn + nIdxJ;

    dNewXJm = m_pdX[nIdxJ - 1] - dSize;
    dNewXJ = m_pdX[nIdxJ] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];

    dMin = Max(dNewXJm, dXIm);
    dMax = Min(dNewXJ, dXI);

    if (dMin < dMax)
    {
      // base function is f_{N1} * (dX - dNewXJm) / (dNewXJ - dNewXJm)
      dTmp = dInvDeltaRow * dInvDeltaColumn;
      sparseMatrix(nIdxRow, nIdxColumn)
          += dTmp * Integrate(dMin, dMax, dNewXJm, dXIm);
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    // dNewXJp = + infinity;

    dMin = Max(dNewXJ, dXIm);
    dMax = dXI; // Min(dNewXJp, dXI);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetRightType() == BCType_Dirichlet)
      {
        // base function is f_{N1}
        sparseMatrix(nIdxRow, nIdxColumn)
            += dInvDeltaRow 
             * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
      }
      else if (boundaryCondition.GetRightType() == BCType_Gamma)
      {
        // base function is f_{N1} * (dCoe * dS - dSN2) / (dSN1 - dSN2)
        //                + f_{N2} * (dCoe * dS - dSN1) / (dSN2 - dSN1)
        dTmp = dInvDeltaRow * dInvDeltaSN2;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * IntegrateExp(dMin, dMax, dCoe, dSN2, dXIm);

        sparseMatrix(nIdxRow, nIdxColumn - 1)
            -= dTmp * IntegrateExp(dMin, dMax, dCoe, dSN1, dXIm);  
      }
    }

    // the internal points of the column
    for (nIdxJ = 1; nIdxJ < m_nNbS - 1; nIdxJ++)
    {
      nIdxColumn = nOffsetColumn + nIdxJ;
      
      dNewXJm = m_pdX[nIdxJ - 1] - dSize;
      dNewXJ = m_pdX[nIdxJ] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];
     
      dMin = Max(dNewXJm, dXIm);
      dMax = Min(dNewXJ, dXI);

      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJm) / (dNewXJ - dNewXJm)
        dTmp = dInvDeltaRow * dInvDeltaColumn;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * Integrate(dMin, dMax, dNewXJm, dXIm);
      }

      dNewXJ = m_pdX[nIdxJ] - dSize;
      dNewXJp = m_pdX[nIdxJ + 1] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

      dMin = Max(dNewXJ, dXIm);
      dMax = Min(dNewXJp, dXI);
     
      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
        dTmp = dInvDeltaRow * dInvDeltaColumn;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * Integrate(dMin, dMax, dNewXJp, dXIm);
      }
    }


    // right support of the row base function
    // the base function is (dX - dXIp) / (dXI - dXIp)
    dXI = m_pdX[nIdxI];
    dXIp = m_pdX[nIdxI + 1];

    dInvDeltaRow = dMultiplier * m_pdInverseDeltaX[nIdxI];
    
    // the left boundary
    nIdxJ = 0;
    nIdxColumn = nOffsetColumn + nIdxJ;

    // dNewXJm = - infinity;
    dNewXJ = m_pdX[nIdxJ] - dSize;

    dMin = dXI; // Max(dNewXJm, dXI); 
    dMax = Min(dNewXJ, dXIp);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
      {
        // base function is f_0 
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dInvDeltaRow 
             * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
      }
      else if (boundaryCondition.GetLeftType() == BCType_Gamma)
      {
        // base function is f_0 * (dCoe * dS - dS1) / (dS0 - dS1)
        //                + f_1 * (dCoe * dS - dS0) / (dS1 - dS0)
        dTmp = dInvDeltaRow * dInvDeltaS0;

        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * IntegrateExp(dMin, dMax, dCoe, dS1, dXIp);

        sparseMatrix(nIdxRow, nIdxColumn + 1)
            -= dTmp * IntegrateExp(dMin, dMax, dCoe, dS0, dXIp);  
      }
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    dNewXJp = m_pdX[nIdxJ + 1] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

    dMin = Max(dNewXJ, dXI);
    dMax = Min(dNewXJp, dXIp);

    if (dMin < dMax)
    {
      // base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
      dTmp = dInvDeltaRow * dInvDeltaColumn;
      sparseMatrix(nIdxRow, nIdxColumn)
          += dTmp * Integrate(dMin, dMax, dNewXJp, dXIp);
    }

    // the right boundary
    nIdxJ = m_nNbS - 1;
    nIdxColumn = nOffsetColumn + nIdxJ;

    dNewXJm = m_pdX[nIdxJ - 1] - dSize;
    dNewXJ = m_pdX[nIdxJ] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];

    dMin = Max(dNewXJm, dXI);
    dMax = Min(dNewXJ, dXIp);

    if (dMin < dMax)
    {
      // base function is f_{N1} * (dX - dNewXJm) / (dNewXJ - dNewXJm)
      dTmp = dInvDeltaRow * dInvDeltaColumn;
      sparseMatrix(nIdxRow, nIdxColumn)
          -= dTmp * Integrate(dMin, dMax, dNewXJm, dXIp);
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    // dNewXJp = + inifinity;

    dMin = Max(dNewXJ, dXI);
    dMax = dXIp; // Min(dNewXJp, dXIp);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetRightType() == BCType_Dirichlet)
      {
        // base function is f_{N1}
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dInvDeltaRow 
             * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
      }
      else if (boundaryCondition.GetRightType() == BCType_Gamma)
      {
        // base function is f_{N1} * (dCoe * dS - dSN2) / (dSN1 - dSN2)
        //                + f_{N2} * (dCoe * dS - dSN1) / (dSN2 - dSN1)
        dTmp = dInvDeltaRow * dInvDeltaSN2;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * IntegrateExp(dMin, dMax, dCoe, dSN2, dXIp);

        sparseMatrix(nIdxRow, nIdxColumn - 1)
            += dTmp * IntegrateExp(dMin, dMax, dCoe, dSN1, dXIp);  
      }
    }

    // the internal points of the column
    for (nIdxJ = 1; nIdxJ < m_nNbS - 1; nIdxJ++)
    { 
      nIdxColumn = nOffsetColumn + nIdxJ;
      
      dNewXJm = m_pdX[nIdxJ - 1] - dSize;
      dNewXJ = m_pdX[nIdxJ] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];
   
      dMin = Max(dNewXJm, dXI);
      dMax = Min(dNewXJ, dXIp);

      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJm) / (dNewXJ - dNewXJm)
        dTmp = dInvDeltaRow * dInvDeltaColumn;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * Integrate(dMin, dMax, dNewXJm, dXIp);
      }

      dNewXJ = m_pdX[nIdxJ] - dSize;
      dNewXJp = m_pdX[nIdxJ + 1] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

      dMin = Max(dNewXJ, dXI);
      dMax = Min(dNewXJp, dXIp);
     
      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
        dTmp = dInvDeltaRow * dInvDeltaColumn;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * Integrate(dMin, dMax, dNewXJp, dXIp);
      }
    }
  } // loop over space nodes
}

// A row or stupidest implementation of the jump system for sensitivity.
// It may be merged into the building of the main jump system
void
StepperFE::BuildAmplitudeSensitivityJumpSystem
           (size_t nIdxR1, size_t nIdxR2,
            double dAmplitude, double dIntensity,
            MorseMatrix& sparseMatrix)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();  

  const size_t nOffsetRow = nIdxR1 * m_nNbS;
  const size_t nOffsetColumn = nIdxR2 * m_nNbS;

  const double dSize = log(1. + dAmplitude);

  const double dS0 = m_instdata.m_pdS[0], dS1 = m_instdata.m_pdS[1],
               dSN2 = m_instdata.m_pdS[m_nNbS - 2], 
               dSN1 = m_instdata.m_pdS[m_nNbS - 1],
               dInvDeltaS0 = 1. / (dS1 - dS0), 
               dInvDeltaSN2 = 1. / (dSN1 - dSN2); 

  double dMin, dMax, dTmp;
  size_t nIdxI, nIdxJ;
  double dNewXJm, dNewXJ, dNewXJp;
  double dXIm, dXI, dXIp;
  double dInvDeltaColumn;
  
  size_t nIdxColumn;

  double dInvDeltaRow;
  const double dInvCoe = 1. / (1 + dAmplitude);
       
  // The interval points for the row, we don't care about the boundary 
  // points
  for (nIdxI = 1; nIdxI < m_nNbS - 1; nIdxI++)
  {
    // The current row
    const size_t nIdxRow = nOffsetRow + nIdxI;

    // left support of the row base function
    // the base function is (dX - dXIm) / (dXI - dXIm)
    dXIm = m_pdX[nIdxI - 1];
    dXI = m_pdX[nIdxI];

    dInvDeltaRow = dIntensity * m_pdInverseDeltaX[nIdxI - 1];
    
    // the left boundary
    nIdxJ = 0;
    nIdxColumn = nOffsetColumn + nIdxJ;

    // dNewXJm = - infinity;
    dNewXJ = m_pdX[nIdxJ] - dSize;

    dMin = dXIm; // Max(dNewXJm, dXIm);
    dMax = Min(dNewXJ, dXI);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
      {
        // nothing to add
      }
      else if (boundaryCondition.GetLeftType() == BCType_Gamma)
      {
        // base function is f_0 * (dCoe * dS - dS1) / (dS0 - dS1)
        //                + f_1 * (dCoe * dS - dS0) / (dS1 - dS0)
        dTmp = dInvDeltaRow * dInvDeltaS0;

        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIm);

        sparseMatrix(nIdxRow, nIdxColumn + 1)
            += dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIm);  
      }
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    dNewXJp = m_pdX[nIdxJ + 1] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

    dMin = Max(dNewXJ, dXIm);
    dMax = Min(dNewXJp, dXI);

    if (dMin < dMax)
    {
      // base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
      dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
      sparseMatrix(nIdxRow, nIdxColumn)
          -= dTmp * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
    }

    // the right boundary
    nIdxJ = m_nNbS - 1;
    nIdxColumn = nOffsetColumn + nIdxJ;

    dNewXJm = m_pdX[nIdxJ - 1] - dSize;
    dNewXJ = m_pdX[nIdxJ] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];

    dMin = Max(dNewXJm, dXIm);
    dMax = Min(dNewXJ, dXI);

    if (dMin < dMax)
    {
      // base function is f_{N1} * (dX - dNewXJm) / (dNewXJ - dNewXJm)
      dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
      sparseMatrix(nIdxRow, nIdxColumn)
          += dTmp * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    // dNewXJp = + infinity;

    dMin = Max(dNewXJ, dXIm);
    dMax = dXI; // Min(dNewXJp, dXI);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetRightType() == BCType_Dirichlet)
      {
        // nothing to add
      }
      else if (boundaryCondition.GetRightType() == BCType_Gamma)
      {
        // base function is f_{N1} * (dCoe * dS - dSN2) / (dSN1 - dSN2)
        //                + f_{N2} * (dCoe * dS - dSN1) / (dSN2 - dSN1)
        dTmp = dInvDeltaRow * dInvDeltaSN2;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIm);

        sparseMatrix(nIdxRow, nIdxColumn - 1)
            -= dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIm);  
      }
    }

    // the internal points of the column
    for (nIdxJ = 1; nIdxJ < m_nNbS - 1; nIdxJ++)
    {
      nIdxColumn = nOffsetColumn + nIdxJ;
      
      dNewXJm = m_pdX[nIdxJ - 1] - dSize;
      dNewXJ = m_pdX[nIdxJ] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];
  
      dMin = Max(dNewXJm, dXIm);
      dMax = Min(dNewXJ, dXI);

      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJm) / (dNewXJ - dNewXJm)
        dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
      }

      dNewXJ = m_pdX[nIdxJ] - dSize;
      dNewXJp = m_pdX[nIdxJ + 1] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

      dMin = Max(dNewXJ, dXIm);
      dMax = Min(dNewXJp, dXI);
      
      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
        dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
      }
    }


    // right support of the row base function
    // the base function is (dX - dXIp) / (dXI - dXIp)
    dXI = m_pdX[nIdxI];
    dXIp = m_pdX[nIdxI + 1];

    dInvDeltaRow = dIntensity * m_pdInverseDeltaX[nIdxI];
    
    // the left boundary
    nIdxJ = 0;
    nIdxColumn = nOffsetColumn + nIdxJ;

    // dNewXJm = - infinity;
    dNewXJ = m_pdX[nIdxJ] - dSize;

    dMin = dXI; // Max(dNewXJm, dXI); 
    dMax = Min(dNewXJ, dXIp);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
      {
        // nothing to add
      }
      else if (boundaryCondition.GetLeftType() == BCType_Gamma)
      {
        // base function is f_0 * (dCoe * dS - dS1) / (dS0 - dS1)
        //                + f_1 * (dCoe * dS - dS0) / (dS1 - dS0)
        dTmp = dInvDeltaRow * dInvDeltaS0;

        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIp);

        sparseMatrix(nIdxRow, nIdxColumn + 1)
            -= dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIp);  
      }
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    dNewXJp = m_pdX[nIdxJ + 1] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

    dMin = Max(dNewXJ, dXI);
    dMax = Min(dNewXJp, dXIp);

    if (dMin < dMax)
    {
      // base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
      dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
      sparseMatrix(nIdxRow, nIdxColumn)
          += dTmp * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
    }

    // the right boundary
    nIdxJ = m_nNbS - 1;
    nIdxColumn = nOffsetColumn + nIdxJ;

    dNewXJm = m_pdX[nIdxJ - 1] - dSize;
    dNewXJ = m_pdX[nIdxJ] - dSize;
    dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];

    dMin = Max(dNewXJm, dXI);
    dMax = Min(dNewXJ, dXIp);

    if (dMin < dMax)
    {
      // base function is f_{N1} * (dX - dNewXJm) / (dNewXJ - dNewXJm)
      dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
      sparseMatrix(nIdxRow, nIdxColumn)
          -= dTmp * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
    }

    dNewXJ = m_pdX[nIdxJ] - dSize;
    // dNewXJp = + inifinity;

    dMin = Max(dNewXJ, dXI);
    dMax = dXIp; // Min(dNewXJp, dXIp);

    if (dMin < dMax)
    {
      if (boundaryCondition.GetRightType() == BCType_Dirichlet)
      {
        // nothing to add
      }
      else if (boundaryCondition.GetRightType() == BCType_Gamma)
      {
        // base function is f_{N1} * (dCoe * dS - dSN2) / (dSN1 - dSN2)
        //                + f_{N2} * (dCoe * dS - dSN1) / (dSN2 - dSN1)
        dTmp = dInvDeltaRow * dInvDeltaSN2;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIp);

        sparseMatrix(nIdxRow, nIdxColumn - 1)
            += dTmp * IntegrateExp(dMin, dMax, 1, 0, dXIp);  
      }
    }

    // the internal points of the column
    for (nIdxJ = 1; nIdxJ < m_nNbS - 1; nIdxJ++)
    {
      nIdxColumn = nOffsetColumn + nIdxJ;
      
      dNewXJm = m_pdX[nIdxJ - 1] - dSize;
      dNewXJ = m_pdX[nIdxJ] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ - 1];    
      
      dMin = Max(dNewXJm, dXI);
      dMax = Min(dNewXJ, dXIp);
      
      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJm) / (dNewXJ - dNewXJm)
        dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
      }

      dNewXJ = m_pdX[nIdxJ] - dSize;
      dNewXJp = m_pdX[nIdxJ + 1] - dSize;
      dInvDeltaColumn = m_pdInverseDeltaX[nIdxJ];

      dMin = Max(dNewXJ, dXI);
      dMax = Min(dNewXJp, dXIp);
      
      if (dMin < dMax)
      {
        // Base function is f_J * (dX - dNewXJp) / (dNewXJ - dNewXJp)
        dTmp = dInvDeltaRow * dInvDeltaColumn * dInvCoe;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
      }
    }
  } // loop over space nodes

}

AutoPtr<MorseMatrix>
StepperFE::BuildAmplitudeSensitivityJumpSystem
(size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity)
{
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  BuildSensitivityJumpSystemStructure
  ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

  shared_ptr<MorseStruct> 
    pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );

  AutoPtr<MorseMatrix> pSparseMatrix(new MorseMatrix);
  MorseMatrix& sparseMatrix = *pSparseMatrix;

  sparseMatrix.Init(pMS, 0.);   

  BuildAmplitudeSensitivityJumpSystem
  ( nIdxR1, nIdxR2, dAmplitude, dIntensity, sparseMatrix);

  return pSparseMatrix;
}

void StepperFE::ComputeHelperArrays()
{
  
}


} // namespace hg

} // namespace ito33
