/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/stepper_fd.cpp
// Purpose:     base stepper class for HG (finite differences)
// Created:     2005/01/13
// RCS-ID:      $Id: stepper_fd.cpp,v 1.24 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/tridiagonalmatrix.h"
#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/boundary1d.h"

#include "hg/model.h"
#include "hg/stepper.h"

namespace ito33
{

namespace hg
{
  using namespace numeric;


void Stepper::Init()
{
  const Model& model = m_instdata.GetModel();

  // Store the pointer to model for quick access
  // We can't use const reference since it would then need to be done in ctor
  m_pModel = &(m_instdata.GetModel());

  // store some vectors for quick access
  m_nNbRegimes = model.GetNbRegimes();

  m_pdVols = model.GetVolatilities();

  m_pdDefaultIntensities = model.GetJumpsToDefault();

  // Time independent cefficient arrays of the particular PDE being solved
  m_pdCoe2nd0 = Array<double>(m_nNbRegimes);
  m_pdCoe1st0 = Array<double>(m_nNbRegimes);
  m_pdCoeZero0 = Array<double>(m_nNbRegimes);
 
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    m_pdCoe2nd0[nIdxR1] = 0.5 * m_pdVols[nIdxR1] * m_pdVols[nIdxR1];

    m_pdCoe1st0[nIdxR1] = 0.0;
    m_pdCoeZero0[nIdxR1] = 0.0;
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_pModel->GetJumps(nIdxR1, nIdxR2);

      Jumps::const_iterator pJump;
      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
      {
        m_pdCoe1st0[nIdxR1] -= pJump->GetIntensity() * pJump->GetAmplitude();
        m_pdCoeZero0[nIdxR1] += pJump->GetIntensity();  
      }                                     
    }

    double dIntensity = m_pdDefaultIntensities[nIdxR1];
    m_pdCoe1st0[nIdxR1] += dIntensity;
    m_pdCoeZero0[nIdxR1] += dIntensity;
  }

  // Coefficent array of the PDE at each time step 
  m_pdCoe2nd = Array<double>(m_nNbRegimes);
  m_pdCoe1st = Array<double>(m_nNbRegimes);
  m_pdCoeZero = Array<double>(m_nNbRegimes);
}

void Stepper::Alloc(size_t nNb)
{  
  pricing::StepperTriDiag::Alloc(nNb);

  // Helper arrays related to the grid for finite difference discretization
  // We are assuming that regimes use the same grid
  m_pdDeltaX = Array<double>(m_nNbS - 1);
  m_pdAreaX = Array<double>(m_nNbS);
  
  m_pdInverseDeltaX = Array<double>(m_nNbS - 1);
  m_pdInverseAreaX = Array<double>(m_nNbS);
  m_pdS2 = Array<double>(m_nNbS);

  // Helper array for building right hand side
  m_pdCoeConst = Array<double>(nNb);

  m_piFD = Array<int>(nNb);

  // Discretization arrays, constant for the price and Greek PDEs
  m_pdAlpha = Array<double>(nNb);
  m_pdBeta = Array<double>(nNb);
}

void Stepper::CalculateAreaArrays(const double* pdX, size_t nNbX)
{
  // The mesh and its size must be initialized by the derived class, since
  // this class doesn't always know what type of grid is being used
  size_t nIdxS;

  m_pdDeltaX[0] = pdX[1] - pdX[0];
  m_pdAreaX[0] = m_pdDeltaX[0];
  for (nIdxS = 1; nIdxS < nNbX - 1; nIdxS++)
  {
    m_pdDeltaX[nIdxS] = pdX[nIdxS + 1] - pdX[nIdxS];
    m_pdAreaX[nIdxS] = (pdX[nIdxS + 1] - pdX[nIdxS-1]) * 0.5;
  }
  m_pdAreaX[nNbX - 1] = m_pdDeltaX[nNbX-2];

  for (nIdxS = 0; nIdxS < nNbX - 1; nIdxS++)
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
  // only implicit or three level will be implemented
  ASSERT(   m_instdata.m_schemeType == SchemeType_Implicit
         || m_instdata.m_schemeType == SchemeType_ThreeLevel );

  switch(m_instdata.m_schemeType)
  {
  case numeric::SchemeType_Implicit:
    {
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdRHS[nIdx] = pdPrice1[nIdx] * m_instdata.m_dOldTimeWeight
                    + m_pdCoeConst[nIdx];
  
    break;
    }
  case numeric::SchemeType_ThreeLevel:
    {
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_pdRHS[nIdx] = pdPrice1[nIdx] * m_instdata.m_dOldTimeWeight
                    + pdPriceM[nIdx] * m_instdata.m_dOldOldTimeWeight
                    + m_pdCoeConst[nIdx];
    break;
    }
  }

  const numeric::Boundary1D &
    boundaryCondition = m_instdata.GetBoundaryCondition();

  // handle boundary condition
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  { 
    // This can't be done before the above loop, otherwise, it will be reset.
    if (m_instdata.GetBoundaryCondition().GetLeftType() == BCType_Dirichlet)
      m_pdRHS[nIdxR * m_nNbS] = boundaryCondition.GetLeftValue();

    if (m_instdata.GetBoundaryCondition().GetRightType() == BCType_Dirichlet)
      m_pdRHS[(nIdxR + 1) * m_nNbS - 1] = boundaryCondition.GetRightValue();
  }
}

void Stepper::BuildMatrix()
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    // Get the diagonals of the tridiagonal matrix for easy access
    size_t nOffSet = nIdxR * m_nNbS;

    double* pdA = m_tridiagonalMatrix.GetA() + nOffSet;
    double* pdB = m_tridiagonalMatrix.GetB() + nOffSet;
    double* pdC = m_tridiagonalMatrix.GetC() + nOffSet;
    double* pdAlpha = m_pdAlpha.Get() + nOffSet;
    double* pdBeta = m_pdBeta.Get() + nOffSet; 
    double dCoeZero =  m_pdCoeZero[nIdxR]; 

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
      pdB[0] = m_instdata.m_dTimeWeight + pdBeta[0] + dCoeZero;
      pdC[0] = - pdBeta[0];
    }
    else
    {
      ASSERT_MSG(false, "Unknown boundary type");
    }

    // handle all interior points. Alpha and beta should already be constructed.
    // m_pdCoeZero entries change. The matrix should be an M-matrix
    size_t nIdx;
    for (nIdx = 1; nIdx < m_nNbS - 1; nIdx++)
    {
      pdA[nIdx] = - pdAlpha[nIdx];
      pdB[nIdx] = m_instdata.m_dTimeWeight
                + pdAlpha[nIdx] + pdBeta[nIdx] + dCoeZero;
      pdC[nIdx] = - pdBeta[nIdx];
    }

    // handle right boundary
    if (boundaryCondition.GetRightType() == BCType_Dirichlet)
    {
      pdA[nIdx] = 0.0;
      pdB[nIdx] = 1.0;
      pdC[nIdx] = 0.0;
    }
    else if (boundaryCondition.GetRightType() == BCType_Gamma)
    {
      pdA[nIdx] = - pdAlpha[nIdx];
      pdB[nIdx] = m_instdata.m_dTimeWeight + pdAlpha[nIdx] + dCoeZero;
      pdC[nIdx] = 0.0;
    }
    else
    {
      ASSERT_MSG(false, "Unknown boundary type");
    }
  }
}

// Build the alpha and beta entries.  These stay the same for the price
// and Greek PDEs, so only need to be made once
void Stepper::BuildAlphaBeta()
{
  // handle all interior points
  size_t nIdx;

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    double dCoe2nd = m_pdCoe2nd[nIdxR];
    double dCoe1st = m_pdCoe1st[nIdxR];
    
    double *pdAlpha = m_pdAlpha.Get() + nIdxR * m_nNbS;
    double *pdBeta = m_pdBeta.Get() + nIdxR * m_nNbS;
    int *piFD = m_piFD.Get() + nIdxR * m_nNbS;
    for (nIdx = 0; nIdx < m_nNbS; nIdx++)
      piFD[nIdx] = 1;
    
    piFD[0] = 2;
    piFD[m_nNbS - 1] = 0;

    // Assume linear boundary condition at left
    nIdx = 0;
    pdAlpha[nIdx] = 0.0;
    pdBeta[nIdx] = dCoe1st * m_pdX[nIdx] * m_dLBC_CoefLeft;

    for (nIdx = 1; nIdx < m_nNbS - 1; nIdx++)
    {
      // Construct coefficients using central weighting. Alpha is for U_{i-1},
      // while beta is for U_{i+1}
      
      // Try central weighting first, since it is 2nd order accurate
      double dTmp1 = dCoe2nd * m_pdS2[nIdx];
      double dTmp2 = dCoe1st * m_pdX[nIdx] * 0.5;

      pdAlpha[nIdx] = (dTmp1 * m_pdInverseDeltaX[nIdx - 1] - dTmp2) 
                    * m_pdInverseAreaX[nIdx];

      pdBeta[nIdx] = (dTmp1 * m_pdInverseDeltaX[nIdx] + dTmp2)
                   * m_pdInverseAreaX[nIdx];

      // Switch to forward differencing, if necessary
      if ( pdAlpha[nIdx] < 0.0 )
      {
        // Upstream weighting can introduce oscillations in delta and gamma. 
        // Use forward differencing
        pdAlpha[nIdx] = dTmp1 * m_pdInverseDeltaX[nIdx - 1]
                      * m_pdInverseAreaX[nIdx];

        pdBeta[nIdx] = (dTmp1 * m_pdInverseAreaX[nIdx] + dCoe1st * m_pdX[nIdx])
                     * m_pdInverseDeltaX[nIdx];

        piFD[nIdx] = 2;
      }


      // Switch to backward differencing, if necessary
      if ( pdBeta[nIdx] < 0.0 )
      {
        // Use backward differencing, instead of downstream weighting
        pdAlpha[nIdx] = (dTmp1 * m_pdInverseAreaX[nIdx] - dCoe1st*m_pdX[nIdx])
                      * m_pdInverseDeltaX[nIdx-1];

        pdBeta[nIdx] = dTmp1 * m_pdInverseDeltaX[nIdx]
                     * m_pdInverseAreaX[nIdx];

        piFD[nIdx] = 0;
      }

    } // loop over internal nodes
 
    // Assume linear boundary condition at right
    pdAlpha[nIdx] = - dCoe1st * m_pdX[nIdx] * m_dLBC_CoefRight;
    pdBeta[nIdx] = 0.0;
  }
}

void Stepper::BuildJumpSystem()
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  // Used to find the columns of the sparse matrix. Lists are used since
  // we don't know how many nonzero entries appear in each row
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  int iNbInterval = (int) m_nNbS - 1;

  // Pass 1. Get the structure of the sparse matrix
  size_t nIdxR1, nIdxR2;
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& pJumps = m_pModel->GetJumps(nIdxR1, nIdxR2);
      size_t nOffset = nIdxR2 * m_nNbS;

      // loop over all jumps from regime I to regime J
      Jumps::const_iterator pJump;
      for (pJump = pJumps.begin(); pJump != pJumps.end(); ++pJump)
      {
        double dAmplitude = pJump->GetAmplitude();

        // loop over all nodes in the space mesh, checking where they jump to
        // assume all regimes use the same grid
        size_t nIdxS = 0;
        if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
          nIdxS = 1;

        size_t nEnd = m_nNbS;
        if (boundaryCondition.GetRightType() == BCType_Dirichlet)
          nEnd = m_nNbS - 1;

        // The interval that a point jumps to
        int iInterval = -1;

        // The current row in the matrix
        size_t nRow = nIdxR1 * m_nNbS + nIdxS;

        for (; nIdxS < nEnd; nIdxS++)
        {
          // find where the point jumps to
          double dNewX = m_pdX[nIdxS] * (1.0 + dAmplitude);
 
          // get the interval. Make sure to not go past the end
          if (iInterval < iNbInterval)
          {
            while (dNewX > m_pdX[iInterval + 1])
            {
              iInterval++;
              if (iInterval == iNbInterval)
                break;
            }
          }

          // Add to the list for the current row
          if (iInterval == -1) // to the left of the grid
          {
            ppnColumnLists[nRow].push_back(nOffset);
            if (boundaryCondition.GetLeftType() == BCType_Gamma)
              ppnColumnLists[nRow].push_back(nOffset + 1);
          }
          else if (iInterval == iNbInterval) // to the right of the grid
          {
            ppnColumnLists[nRow].push_back(nOffset + m_nNbS - 1);
            if (boundaryCondition.GetRightType() == BCType_Gamma)
              ppnColumnLists[nRow].push_back(nOffset + m_nNbS - 2);
          }
          else // internal
          {
            ppnColumnLists[nRow].push_back(nOffset + iInterval);
            ppnColumnLists[nRow].push_back(nOffset + iInterval + 1);
            // extra for quadratic interpolation
            if (iInterval < iNbInterval - 1)
              ppnColumnLists[nRow].push_back(nOffset + iInterval + 2);
            else
              ppnColumnLists[nRow].push_back(nOffset + iInterval - 1);
          }

          // Move to next row
          nRow++;

        } // loop over space nodes
      } // loop over jumps from regime I to regime J
    } // inner regime loop
  } // outer regime loop

  shared_ptr<MorseStruct> 
    pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );

  m_sparseMatrix.Init(pMS, 0.);   

  // Pass 2. Fill the sparse matrix
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& pJumps = m_pModel->GetJumps(nIdxR1, nIdxR2);
      size_t nOffset = nIdxR2 * m_nNbS;

      // loop over all jumps from regime I to regime J
      Jumps::const_iterator pJump;

      for (pJump = pJumps.begin(); pJump != pJumps.end(); ++pJump)
      {
        double dAmplitude = pJump->GetAmplitude();
        double dIntensity = pJump->GetIntensity();

        // loop over all nodes in the space mesh, checking where they jump to
        // assume all regimes use the same grid
        size_t nIdxS = 0;
        if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
          nIdxS = 1;

        size_t nEnd = m_nNbS;
        if (boundaryCondition.GetRightType() == BCType_Dirichlet)
          nEnd = m_nNbS - 1;
        
        // The interval that a point jumps to
        int iInterval = -1;

        // The current row in the matrix
        size_t nRow = nIdxR1 * m_nNbS + nIdxS;

        for (; nIdxS < nEnd; nIdxS++)
        {
          // find where the point jumps to
          double dNewX = m_pdX[nIdxS] * (1.0 + dAmplitude);
 
          // get the interval
          if (iInterval < iNbInterval)
          {
            while (dNewX > m_pdX[iInterval + 1])
            {
              iInterval++;
              if (iInterval == iNbInterval)
                break;
            }
          }

          // use quadratic interpolation and add to the matrix
          double dX1, dX2;
          if (iInterval == -1) // to the left of the grid
          {
            if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
              m_sparseMatrix(nRow, nOffset) -= dIntensity;
            else if (boundaryCondition.GetLeftType() == BCType_Gamma)
            {
              dX1 = m_pdX[0];
              dX2 = m_pdX[1];
              double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
              double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
              m_sparseMatrix(nRow, nOffset) -= dIntensity * dWeight1;
              m_sparseMatrix(nRow, nOffset + 1) -= dIntensity * dWeight2;
            }
          }
          else if (iInterval == iNbInterval) // to the right of the grid
          {                  
            if (boundaryCondition.GetRightType() == BCType_Dirichlet)
              m_sparseMatrix(nRow, nOffset + iInterval) -= dIntensity;
            else if (boundaryCondition.GetRightType() == BCType_Gamma)
            {
              dX1 = m_pdX[m_nNbS - 1];
              dX2 = m_pdX[m_nNbS - 2];
              double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
              double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
              m_sparseMatrix(nRow, nOffset + m_nNbS - 1) -= dIntensity * dWeight1;
              m_sparseMatrix(nRow, nOffset + m_nNbS - 2) -= dIntensity * dWeight2;
            }              
          }
          else // internal
          {
#if 1
            int iIntervalStart = iInterval;
            if (iInterval >= iNbInterval - 1)
              iIntervalStart = iInterval - 1;

            dX1 = m_pdX[iIntervalStart];
            dX2 = m_pdX[iIntervalStart + 1];
            double dX3 = m_pdX[iIntervalStart + 2];

            double dWeight1 
                 = (dNewX - dX2) * (dNewX - dX3) / ((dX1 - dX2) * (dX1 - dX3));

            double dWeight2 
                 = (dNewX - dX1) * (dNewX - dX3) / ((dX2 - dX1) * (dX2 - dX3));

            double dWeight3 
                 = (dNewX - dX1) * (dNewX - dX2) / ((dX3 - dX1) * (dX3 - dX2));

            m_sparseMatrix(nRow, nOffset + iIntervalStart) 
              -= dIntensity * dWeight1;

            m_sparseMatrix(nRow, nOffset + iIntervalStart + 1) 
              -= dIntensity * dWeight2;

            m_sparseMatrix(nRow, nOffset + iIntervalStart + 2) 
              -= dIntensity * dWeight3;
#else
            // linear
            dX1 = m_pdX[iInterval];
            dX2 = m_pdX[iInterval + 1];

            double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
            double dWeight2 = (dNewX - dX1) / (dX2 - dX1);

            m_sparseMatrix(nRow, nOffset + iInterval)
              -= dIntensity * dWeight1;

            m_sparseMatrix(nRow, nOffset + iInterval + 1) 
              -= dIntensity * dWeight2;
#endif
          }
           
          // Move to next row
          nRow++;

        } // loop over space nodes
      } // loop over jumps from regime I to regime J
    } // inner regime loop
  } // outer regime loop

  // m_sparseMatrix.Dump();
}


AutoPtr<numeric::MorseMatrix>
Stepper::BuildJumpSensitivitySystem
(size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity,
 bool bIsAmplitude)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  // Used to find the columns of the sparse matrix. Lists are used since
  // we don't know how many nonzero entries appear in each row
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  int iNbInterval = (int) m_nNbS - 1;

  // Pass 1. Get the structure of the sparse matrix
  size_t nOffset = nIdxR2 * m_nNbS;

  // loop over all nodes in the space mesh, checking where they jump to
  // assume all regimes use the same grid
  size_t nIdxS = 0;
  if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
    nIdxS = 1;

  size_t nEnd = m_nNbS;
  if (boundaryCondition.GetRightType() == BCType_Dirichlet)
    nEnd = m_nNbS - 1;

  // The interval that a point jumps to
  int iInterval = -1;

  // The current row in the matrix
  size_t nRow = nIdxR1 * m_nNbS + nIdxS;
 
  for (; nIdxS < nEnd; nIdxS++)
  {
    // find where the point jumps to
    double dNewX = m_pdX[nIdxS] * (1.0 + dAmplitude);
 
    // get the interval. Make sure to not go past the end
    if (iInterval < iNbInterval)
    {
      while (dNewX > m_pdX[iInterval + 1])
      {
        iInterval++;
        if (iInterval == iNbInterval)
          break;
      }
    }

    // Add to the list for the current row
    if (iInterval == -1) // to the left of the grid
    {
      ppnColumnLists[nRow].push_back(nOffset);
      if (boundaryCondition.GetLeftType() == BCType_Gamma)
        ppnColumnLists[nRow].push_back(nOffset + 1);
    }
    else if (iInterval == iNbInterval) // to the right of the grid
    {
      ppnColumnLists[nRow].push_back(nOffset + m_nNbS - 1);
      if (boundaryCondition.GetRightType() == BCType_Gamma)
        ppnColumnLists[nRow].push_back(nOffset + m_nNbS - 2);
    }
    else // internal
    {
      ppnColumnLists[nRow].push_back(nOffset + iInterval);
      ppnColumnLists[nRow].push_back(nOffset + iInterval + 1);
      // extra for quadratic interpolation
      if (iInterval < iNbInterval - 1)
        ppnColumnLists[nRow].push_back(nOffset + iInterval + 2);
      else
        ppnColumnLists[nRow].push_back(nOffset + iInterval - 1);
    }

    // Move to next row
    nRow++;

  } // loop over space nodes

  shared_ptr<MorseStruct> 
    pMS( new MorseStruct( ppnColumnLists.Get(), m_nNbX) );

  AutoPtr<MorseMatrix> pSparseMatrix( new MorseMatrix(pMS) );   

  // Pass 2. Fill the sparse matrix

  // loop over all nodes in the space mesh, checking where they jump to
  // assume all regimes use the same grid
  nIdxS = 0;
  if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
    nIdxS = 1;

  nEnd = m_nNbS;
  if (boundaryCondition.GetRightType() == BCType_Dirichlet)
    nEnd = m_nNbS - 1;
        
  // The interval that a point jumps to
  iInterval = -1;

  // The current row in the matrix
  nRow = nIdxR1 * m_nNbS + nIdxS;

  for (; nIdxS < nEnd; nIdxS++)
  {
    // find where the point jumps to
    double dX = m_pdX[nIdxS];
    double dNewX = dX * (1.0 + dAmplitude);
 
    // get the interval
    if (iInterval < iNbInterval)
    {
      while (dNewX > m_pdX[iInterval + 1])
      {
        iInterval++;
        if (iInterval == iNbInterval)
          break;
      }
    }

    // use quadratic interpolation and add to the matrix
    double dX1, dX2, dX3;
    if (iInterval == -1) // to the left of the grid
    {
      if (boundaryCondition.GetLeftType() == BCType_Dirichlet)
      {
        if (bIsAmplitude)
          (*pSparseMatrix)(nRow, nOffset) = 0.0;
        else
          (*pSparseMatrix)(nRow, nOffset) = 1.0;
      }
      else if (boundaryCondition.GetLeftType() == BCType_Gamma)
      {
        if (bIsAmplitude)
        {
          dX1 = m_pdX[0];
          dX2 = m_pdX[1];
          double dWeight1 = dX / (dX1 - dX2);
          double dWeight2 = dX / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset) = dIntensity * dWeight1;
          (*pSparseMatrix)(nRow, nOffset + 1) = dIntensity * dWeight2;
        }
        else
        {
          dX1 = m_pdX[0];
          dX2 = m_pdX[1];
          double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
          double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset) = dWeight1;
          (*pSparseMatrix)(nRow, nOffset + 1) = dWeight2;
        }
      }
    }
    else if (iInterval == iNbInterval) // to the right of the grid
    {                  
      if (boundaryCondition.GetRightType() == BCType_Dirichlet)
        if (bIsAmplitude)
          (*pSparseMatrix)(nRow, nOffset + iInterval) = 0.0;
        else
          (*pSparseMatrix)(nRow, nOffset + iInterval) = 1.0;
      else if (boundaryCondition.GetRightType() == BCType_Gamma)
      {
        if (bIsAmplitude)
        {
          dX1 = m_pdX[m_nNbS - 1];
          dX2 = m_pdX[m_nNbS - 2];
          double dWeight1 = dX / (dX1 - dX2);
          double dWeight2 = dX / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 1) = dIntensity * dWeight1;
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 2) = dIntensity * dWeight2;
        }
        else
        {
          dX1 = m_pdX[m_nNbS - 1];
          dX2 = m_pdX[m_nNbS - 2];
          double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
          double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 1) = dWeight1;
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 2) = dWeight2;
        }
      }              
    }
    else // internal
    {
      int iIntervalStart = iInterval;
      if (iInterval >= iNbInterval - 1)
        iIntervalStart = iInterval - 1;

      dX1 = m_pdX[iIntervalStart];
      dX2 = m_pdX[iIntervalStart + 1];
      dX3 = m_pdX[iIntervalStart + 2];

      if (bIsAmplitude)
      {
        double dWeight1 
           = ( dX * (2.0*dNewX - dX3 - dX2) ) / ((dX1 - dX2) * (dX1 - dX3));

        double dWeight2 
           = ( dX * (2.0*dNewX - dX1 - dX3) ) / ((dX2 - dX1) * (dX2 - dX3));

        double dWeight3 
           = ( dX * (2.0*dNewX - dX1 - dX2) ) / ((dX3 - dX1) * (dX3 - dX2));

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart)    
          = dIntensity * dWeight1;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 1)
          = dIntensity * dWeight2;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 2)
          = dIntensity * dWeight3;
      }
      else
      {
        // Must be intensity
        double dWeight1 
             = (dNewX - dX2) * (dNewX - dX3) / ((dX1 - dX2) * (dX1 - dX3));

        double dWeight2 
             = (dNewX - dX1) * (dNewX - dX3) / ((dX2 - dX1) * (dX2 - dX3));

        double dWeight3 
             = (dNewX - dX1) * (dNewX - dX2) / ((dX3 - dX1) * (dX3 - dX2));

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart) = dWeight1;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 1) = dWeight2;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 2) = dWeight3;
      }

    }
           
    // Move to next row
    nRow++;

  } // loop over space nodes

  return pSparseMatrix;
}

void Stepper::ComputeHelperArrays()
{
  const numeric::Boundary1D &
    boundaryCondition = m_instdata.GetBoundaryCondition();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    double* pdPrices = m_instdata.m_pdPrices.Get() + nIdxR * m_nNbS;
    double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;
    double* pdGammas = m_instdata.m_pdGammas.Get() + nIdxR * m_nNbS;
    int* piFD = m_piFD.Get() + nIdxR * m_nNbS;
    
    ComputeGammaFD(m_pdX, pdPrices, m_nNbS, pdGammas);

    if ( boundaryCondition.GetLeftType() == BCType_Gamma )
      pdGammas[0] = 0.;
    if ( boundaryCondition.GetRightType() == BCType_Gamma )
      pdGammas[m_nNbS - 1] = 0.;

    ComputeDelta(m_pdX, piFD, pdPrices, m_nNbS, pdDeltas);
  }
}


} // namespace hg

} // namespace ito33
