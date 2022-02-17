/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/forward/forwardoptionstepper_fd.cpp
// Purpose:     implementation of HG stepper class for forward options
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionstepper_fd.cpp,v 1.5 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"

#include "ito33/numeric/boundary1d.h"

#include "hg/model.h"
#include "hg/forwardoptionstepper_fd.h"

namespace ito33
{

namespace hg
{
  using namespace numeric;


void ForwardOptionStepperFD::MakeCoefficients()
{
  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
    m_pdCoe1st[nIdxR] = -(m_pdCoe1st0[nIdxR] + dRate - dForeignRate);
    m_pdCoeZero[nIdxR] = m_pdCoe1st[nIdxR] + m_pdCoeZero0[nIdxR] + dRate;
 
    
    // recovery is zero for calls
    double dTmp = m_pdDefaultIntensities[nIdxR] * m_instdata.m_dRecoveryValue;

    double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;
    for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
      pdCoeConst[nIdx] = dTmp;
  }
}

void ForwardOptionStepperFD::MakeSensitivityCoefficients
                             ( const SensitivityData& data )
{
  // m_pdCoe2nd, m_pdCoe1st and m_pdCoeZero are the same as for pricing.
  // Looking at the actual HG PDE should make the code below easy to
  // understand. Just derive the equation with respect to the appropriate
  // variable

  // Since the amplitude and intensity parameters can apply to more than one
  // regime, clear the constant coefficient array.  The different cases
  // below can then just update the appropriate regime.  Perhaps not the
  // most efficient approach, but is hopefully clear.
  for (size_t nIdx = 0; nIdx < m_nNbS * m_nNbRegimes; nIdx++)
    m_pdCoeConst[nIdx] = 0.0;

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++) 
  { 
    // get the constant coefficients for this regime
    double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;

    // Determine what type of sensitivity is being computed.
    switch (data.m_sensitivityType)
    {
    case SensitivityType_Volatility:
      {
      if (data.m_nRegime1 != nIdxR)
        break;

      // volatility coefficient
      double* pdGammas = m_instdata.m_pdGammas.Get() + nIdxR * m_nNbS;

      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] = m_pdVols[nIdxR] * pdGammas[nIdx] * m_pdS2[nIdx];

      break;
      }
    case SensitivityType_DefaultIntensity:
    {
      if (data.m_nRegime1 != nIdxR)
        break;

      double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;

      // default intensity coefficient

      for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] = - m_pdX[nIdx] * pdDeltas[nIdx];

      break;
    }
    case SensitivityType_JumpIntensity:
      {
      // intensity
      if (data.m_nRegime1 == nIdxR)
      {
        double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;
        double* pdPrices = m_instdata.m_pdPrices.Get() + nIdxR * m_nNbS;

        for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
          pdCoeConst[nIdx] += - pdPrices[nIdx]
                            + data.m_dAmplitude * m_pdX[nIdx] * pdDeltas[nIdx]
                            - data.m_dAmplitude * pdPrices[nIdx]; 
      }
      
      if (data.m_nRegime2 == nIdxR)
        data.m_pSparseMatrix->ProductMatrixVector
             ( m_instdata.m_pdPrices.Get(), m_pdCoeConst.Get(), true );
  
      break;
      }
    case SensitivityType_JumpAmplitude:
      {
      // amplitude
      if (data.m_nRegime1 == nIdxR)
      {
        double* pdDeltas = m_instdata.m_pdDeltas.Get() + nIdxR * m_nNbS;
        double* pdPrices = m_instdata.m_pdPrices.Get() + nIdxR * m_nNbS;

        for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
        pdCoeConst[nIdx] += data.m_dIntensity * m_pdX[nIdx] * pdDeltas[nIdx]
                          - data.m_dIntensity * pdPrices[nIdx];
      }
     
      if (data.m_nRegime2 == nIdxR)
        data.m_pSparseMatrix->ProductMatrixVector
             ( m_instdata.m_pdPrices.Get(), m_pdCoeConst.Get(), true );
      
      break;
      }
    } // end the switch

  } // loop over regimes
}


void ForwardOptionStepperFD::BuildJumpSystem()
{
  const numeric::Boundary1D &
    boundaryCondition = m_instdata.GetBoundaryCondition();

  // Used to find the columns of the sparse matrix. Lists are used since
  // we don't know how many nonzero entries appear in each row
  Array< std::list<size_t> > 
    ppnColumnLists = Array< std::list<size_t> >(m_nNbX);

  int iNbInterval = (int) m_nNbS - 1;

  // Pass 1. Get the structure of the sparse matrix
  size_t nIdxR1, nIdxR2;
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& pJumps = m_pModel->GetJumps(nIdxR1, nIdxR2);
      size_t nOffset = nIdxR1 * m_nNbS;

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

        // The current row in the matrix. Note the use of R2.  For the
        // forward equation, we want the row of the matrix that is being
        // jumped to (instead of from).
        size_t nRow = nIdxR2 * m_nNbS + nIdxS;

        for (; nIdxS < nEnd; nIdxS++)
        {
          // find where the point jumps to
          double dNewX = m_pdX[nIdxS] / (1.0 + dAmplitude);
 
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
    pMS( new MorseStruct( ppnColumnLists.Get(), m_nNbX) );

  m_sparseMatrix.Init(pMS, 0.);   

  // Pass 2. Fill the sparse matrix
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& pJumps = m_pModel->GetJumps(nIdxR1, nIdxR2);
      size_t nOffset = nIdxR1 * m_nNbS;

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

        // The current row in the matrix. Note the use of R2.  For the
        // forward equation, we want the row of the matrix that is being
        // jumped to (instead of from).
        size_t nRow = nIdxR2 * m_nNbS + nIdxS;

        for (; nIdxS < nEnd; nIdxS++)
        {
          // find where the point jumps to
          double dNewX = m_pdX[nIdxS] / (1.0 + dAmplitude);
 
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
              m_sparseMatrix(nRow, nOffset) -= dIntensity * (1.0 + dAmplitude);
            else if (boundaryCondition.GetLeftType() == BCType_Gamma)
            {
              dX1 = m_pdX[0];
              dX2 = m_pdX[1];
              double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
              double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
              m_sparseMatrix(nRow, nOffset) -= dIntensity * (1.0 + dAmplitude) * dWeight1;
              m_sparseMatrix(nRow, nOffset + 1) -= dIntensity * (1.0 + dAmplitude) * dWeight2;
            }
          }
          else if (iInterval == iNbInterval) // to the right of the grid
          {                  
            if (boundaryCondition.GetRightType() == BCType_Dirichlet)
              m_sparseMatrix(nRow, nOffset + iInterval) -= dIntensity * (1.0 + dAmplitude);
            else if (boundaryCondition.GetRightType() == BCType_Gamma)
            {
              dX1 = m_pdX[m_nNbS - 1];
              dX2 = m_pdX[m_nNbS - 2];
              double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
              double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
              m_sparseMatrix(nRow, nOffset + m_nNbS - 1) -= dIntensity * (1.0 + dAmplitude) * dWeight1;
              m_sparseMatrix(nRow, nOffset + m_nNbS - 2) -= dIntensity * (1.0 + dAmplitude) * dWeight2;
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

            double dWeight1 
                 = (dNewX - dX2) * (dNewX - dX3) / ((dX1 - dX2) * (dX1 - dX3));

            double dWeight2 
                 = (dNewX - dX1) * (dNewX - dX3) / ((dX2 - dX1) * (dX2 - dX3));

            double dWeight3 
                 = (dNewX - dX1) * (dNewX - dX2) / ((dX3 - dX1) * (dX3 - dX2));

            m_sparseMatrix(nRow, nOffset + iIntervalStart) 
              -= dIntensity * (1.0 + dAmplitude) * dWeight1;

            m_sparseMatrix(nRow, nOffset + iIntervalStart + 1) 
              -= dIntensity * (1.0 + dAmplitude) * dWeight2;

            m_sparseMatrix(nRow, nOffset + iIntervalStart + 2) 
              -= dIntensity * (1.0 + dAmplitude) * dWeight3;

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
ForwardOptionStepperFD::BuildJumpSensitivitySystem(size_t nIdxR1, 
                                                 size_t nIdxR2,
                                                 double dAmplitude, 
                                                 double dIntensity,
                                                 bool bIsAmplitude)
{
  const numeric::Boundary1D &
    boundaryCondition = m_instdata.GetBoundaryCondition();

  // Used to find the columns of the sparse matrix. Lists are used since
  // we don't know how many nonzero entries appear in each row
  Array< std::list<size_t> > 
    ppnColumnLists = Array< std::list<size_t> >(m_nNbX);

  int iNbInterval = (int) m_nNbS - 1;

  // Pass 1. Get the structure of the sparse matrix
  size_t nOffset = nIdxR1 * m_nNbS;

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
  size_t nRow = nIdxR2 * m_nNbS + nIdxS;
 
  for (; nIdxS < nEnd; nIdxS++)
  {
    // find where the point jumps to
    double dNewX = m_pdX[nIdxS] / (1.0 + dAmplitude);
 
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

  AutoPtr<numeric::MorseMatrix> pSparseMatrix = AutoPtr<numeric::MorseMatrix>
    (new numeric::MorseMatrix( pMS ) );   

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
  nRow = nIdxR2 * m_nNbS + nIdxS;

  for (; nIdxS < nEnd; nIdxS++)
  {
    // find where the point jumps to
    double dX = m_pdX[nIdxS];
    double dNewX = dX / (1.0 + dAmplitude);
 
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
          (*pSparseMatrix)(nRow, nOffset) = dIntensity;
        else
          (*pSparseMatrix)(nRow, nOffset) = (1.0 + dAmplitude);
      }
      else if (boundaryCondition.GetLeftType() == BCType_Gamma)
      {
        if (bIsAmplitude)
        {
          dX1 = m_pdX[0];
          dX2 = m_pdX[1];
          double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
          double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
          double dDerivWeight1 = -dNewX / (1.0 + dAmplitude) / (dX1 - dX2);
          double dDerivWeight2 = -dNewX / (1.0 + dAmplitude) / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset) = dIntensity * dWeight1
                                          + dIntensity * (1.0 + dAmplitude) * dDerivWeight1;;
          (*pSparseMatrix)(nRow, nOffset + 1) = dIntensity * dWeight2
                                              + dIntensity * (1.0 + dAmplitude) * dDerivWeight2;
        }
        else
        {
          dX1 = m_pdX[0];
          dX2 = m_pdX[1];
          double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
          double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset) = (1.0 + dAmplitude) * dWeight1;
          (*pSparseMatrix)(nRow, nOffset + 1) = (1.0 + dAmplitude) * dWeight2;
        }
      }
    }
    else if (iInterval == iNbInterval) // to the right of the grid
    {                  
      if (boundaryCondition.GetRightType() == BCType_Dirichlet)
        if (bIsAmplitude)
          (*pSparseMatrix)(nRow, nOffset + iInterval) = dIntensity;
        else
          (*pSparseMatrix)(nRow, nOffset + iInterval) = (1.0 + dAmplitude);
      else if (boundaryCondition.GetRightType() == BCType_Gamma)
      {
        if (bIsAmplitude)
        {
          dX1 = m_pdX[m_nNbS - 1];
          dX2 = m_pdX[m_nNbS - 2];
          double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
          double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
          double dDerivWeight1 = -dNewX / (1.0 + dAmplitude) / (dX1 - dX2);
          double dDerivWeight2 = -dNewX / (1.0 + dAmplitude) / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 1) = dIntensity * dWeight1
                                                       + dIntensity * (1.0 + dAmplitude) * dDerivWeight1;
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 2) = dIntensity * dWeight2
                                                       + dIntensity * (1.0 + dAmplitude) * dDerivWeight2;
        }
        else
        {
          dX1 = m_pdX[m_nNbS - 1];
          dX2 = m_pdX[m_nNbS - 2];
          double dWeight1 = (dNewX - dX2) / (dX1 - dX2);
          double dWeight2 = (dNewX - dX1) / (dX2 - dX1);
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 1) = (1.0 + dAmplitude) * dWeight1;
          (*pSparseMatrix)(nRow, nOffset + m_nNbS - 2) = (1.0 + dAmplitude) * dWeight2;
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
           = (dNewX - dX2) * (dNewX - dX3) / ((dX1 - dX2) * (dX1 - dX3));

        double dWeight2 
           = (dNewX - dX1) * (dNewX - dX3) / ((dX2 - dX1) * (dX2 - dX3));

        double dWeight3 
           = (dNewX - dX1) * (dNewX - dX2) / ((dX3 - dX1) * (dX3 - dX2));

        double dDerivWeight1 = ( -dNewX * (2.0*dNewX - dX3 - dX2) ) / 
                               ( (dX1 - dX2) * (dX1 - dX3) * (1 + dAmplitude) );

        double dDerivWeight2 = ( -dNewX * (2.0*dNewX - dX1 - dX3) ) / 
                               ( (dX2 - dX1) * (dX2 - dX3) * (1 + dAmplitude) );

        double dDerivWeight3 = ( -dNewX * (2.0*dNewX - dX1 - dX2) ) / 
                               ( (dX3 - dX1) * (dX3 - dX2) * (1 + dAmplitude) );


        (*pSparseMatrix)(nRow, nOffset + iIntervalStart) 
          = dIntensity * dWeight1
          + dIntensity * (1.0 + dAmplitude) * dDerivWeight1;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 1) 
          = dIntensity * dWeight2
          + dIntensity * (1.0 + dAmplitude) * dDerivWeight2;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 2) 
          = dIntensity * dWeight3
          + dIntensity * (1.0 + dAmplitude) * dDerivWeight3;

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

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart) = (1.0 + dAmplitude) * dWeight1;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 1) = (1.0 + dAmplitude) * dWeight2;

        (*pSparseMatrix)(nRow, nOffset + iIntervalStart + 2) = (1.0 + dAmplitude) * dWeight3;
      }

    }
           
    // Move to next row
    nRow++;

  } // loop over space nodes

  return pSparseMatrix;
}


} // namespace hg

} // namespace ito33
