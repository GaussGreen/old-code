/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/forward/forwardoptionstepper_fe.cpp
// Purpose:     FE implementation of HG stepper class for forward options
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionstepper_fe.cpp,v 1.12 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"
#include "ito33/numeric/trisparseconstraintsolver_penalty.h"
#include "ito33/numeric/trisparseconstraintsolver_frozen.h"

#include "ito33/numeric/boundary1d.h"

#include "hg/model.h"
#include "hg/forwardoptionstepper_fe.h"

using namespace ito33;
using namespace ito33::numeric;

inline double Sqr(double d) { return d * d; }

inline double Max(double dA, double dB)
{
  return dA > dB ? dA : dB;
}

inline double Min(double dA, double dB)
{
  return dA > dB ? dB : dA;
}


namespace ito33
{

namespace hg
{

extern double Integrate(double dMin, double dMax, double dX1, double dX2);

extern double IntegrateExp(double dMin, double dMax, 
                           double dCoe, double dX1, double dX2);

void ForwardOptionStepperFE::MakeCoefficients()
{
  m_bSensitivityRHS = false;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    m_pdCoe2nd[nIdxR] = m_pdCoe2nd0[nIdxR];
    m_pdCoe1st[nIdxR] = m_pdCoe2nd0[nIdxR] - m_pdCoe1st0[nIdxR]
                      + dRate - dForeignRate;
    m_pdCoeZero[nIdxR] = m_pdCoeZero0[nIdxR] + m_pdCoe1st0[nIdxR]
                       + dForeignRate + m_instdata.m_dTimeWeight;
 
    double dTmp = m_pdDefaultIntensities[nIdxR] * m_instdata.m_dRecoveryValue;

    double* pdCoeConst = m_pdCoeConst.Get() + nIdxR * m_nNbS;
    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      pdCoeConst[nIdxS] = dTmp;
  }
}

void ForwardOptionStepperFE::BuildDerivedJumpSystem()
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
        
        sparseMatrix(nIdxRow, nIdxRow - 1) = 0.5;      
        sparseMatrix(nIdxRow, nIdxRow) = 0.; 
        sparseMatrix(nIdxRow, nIdxRow + 1) = - 0.5;  
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
      
      double dCoeZero = 1. + dAmplitude;
     
      for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
      {  
        const size_t nIdxRow = nOffSetRow + nIdxS;  
        
        double dDeltaX0 = m_pdDeltaX[nIdxS - 1];  
        double dDeltaX = m_pdDeltaX[nIdxS];  
        
        sparseMatrix(nIdxRow, nIdxRow - 1) = - 0.5 * dAmplitude 
                                           - dCoeZero * dDeltaX0 / 6.;      
        sparseMatrix(nIdxRow, nIdxRow) = - dCoeZero * (dDeltaX0 + dDeltaX) / 3.; 
        sparseMatrix(nIdxRow, nIdxRow + 1) = 0.5 * dAmplitude
                                           - dCoeZero * dDeltaX / 6.;  
      }
      
      BuildPartialJumpSystem
      ( nIdxR1, nIdxR2, dAmplitude, dCoeZero, sparseMatrix );

      break;
      }
    case SensitivityType_JumpAmplitude:
      {
      BuildSensitivityJumpSystemStructure
      ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

      shared_ptr<MorseStruct> 
        pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );      
      sparseMatrix.Init(pMS, 0);
      
      const double dCoeZero = dIntensity;
      
      for (size_t nIdxS = 1; nIdxS < m_nNbS - 1; nIdxS++)
      {  
        const size_t nIdxRow = nOffSetRow + nIdxS;  

        double dDeltaX0 = m_pdDeltaX[nIdxS - 1];  
        double dDeltaX = m_pdDeltaX[nIdxS];  
        
        sparseMatrix(nIdxRow, nIdxRow - 1) = - 0.5 * dIntensity
                                           - dCoeZero * dDeltaX0 / 6.;      
        sparseMatrix(nIdxRow, nIdxRow) = - dCoeZero * (dDeltaX0 + dDeltaX) / 3.;
        sparseMatrix(nIdxRow, nIdxRow + 1) = 0.5 * dIntensity
                                           - dCoeZero * dDeltaX / 6.;  
      }
      
      BuildPartialJumpSystem
      ( nIdxR1, nIdxR2, dAmplitude, dCoeZero, sparseMatrix );

      BuildAmplitudeSensitivityJumpSystem
      ( nIdxR1, nIdxR2, dAmplitude, dIntensity, sparseMatrix);

      break;
      }
   } // end the switch
  }
}

void ForwardOptionStepperFE::MakeSensitivityCoefficients
     (const SensitivityData& data)
{
  m_bSensitivityRHS = true;

  for (size_t nIdxX = 0; nIdxX < m_nNbX; nIdxX++)
    m_pdCoeConst[nIdxX] = 0.;

  data.m_pSparseMatrix->ProductMatrixVector
                        (m_instdata.m_pdPrices.Get(), 
                         m_pdCoeConstSensitivity.Get());
}


// A row or stupidest implementation of the jump system. It's quite difficult
// to use it to compute the sensitivity jump system but it can serve as a 
// reference that other implementation can be checked against
void ForwardOptionStepperFE::BuildJumpSystem()
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

        const double dMultiplier = - dIntensity * (1.0 + dAmplitude);
        BuildPartialJumpSystem
        ( nIdxR1, nIdxR2, dAmplitude, dMultiplier, m_sparseMatrix );

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

  // m_sparseMatrix.Dump();
}

void
ForwardOptionStepperFE::BuildSensitivityJumpSystemStructure
                        (size_t nIdxR1, size_t nIdxR2, double dAmplitude,
                         std::list<size_t>* ppnColumnLists)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  Array<double> pdNewX(m_nNbS);
  Array<int> piIV(m_nNbS);

  const int iNbIntervals = int(m_nNbS) - 1;

  const size_t nOffsetColumn = nIdxR1 * m_nNbS;
  const size_t nOffsetRow = nIdxR2 * m_nNbS;
      
  const double dSize = - log(1. + dAmplitude);

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

// A row or stupidest implementation of the jump system
// It may be merged into the building of the main jump system
void
ForwardOptionStepperFE::BuildPartialJumpSystem
           (size_t nIdxR1, size_t nIdxR2,
            double dAmplitude, double dMultiplier,
            MorseMatrix& sparseMatrix)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  const size_t nOffsetColumn = nIdxR1 * m_nNbS;
  const size_t nOffsetRow = nIdxR2 * m_nNbS;
      
  const double dSize = - log(1. + dAmplitude);

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

  const double dCoe = 1. / (1. + dAmplitude);

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

// Just use the helper function
AutoPtr<MorseMatrix>
ForwardOptionStepperFE::BuildIntensitySensitivityJumpSystem
           (size_t nIdxR1, size_t nIdxR2,
            double dAmplitude, double /* dIntensity */)
{
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  BuildSensitivityJumpSystemStructure
  ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

  shared_ptr<MorseStruct> 
    pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );

  AutoPtr<MorseMatrix> pSparseMatrix(new MorseMatrix);
  MorseMatrix& sparseMatrix = *pSparseMatrix;

  sparseMatrix.Init(pMS, 0.);   

  BuildPartialJumpSystem
  ( nIdxR1, nIdxR2, dAmplitude, 1. + dAmplitude, sparseMatrix );

  return pSparseMatrix;
}

// A row or stupidest implementation of the jump system for sensitivity.
// It may be merged into the building of the main jump system
void
ForwardOptionStepperFE::BuildAmplitudeSensitivityJumpSystem
(size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity,
 numeric::MorseMatrix& sparseMatrix)
{
  const Boundary1D& boundaryCondition = m_instdata.GetBoundaryCondition();

  const size_t nOffsetColumn = nIdxR1 * m_nNbS;
  const size_t nOffsetRow = nIdxR2 * m_nNbS;
      
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

  const double dSize = - log(1. + dAmplitude);
  const double dCoe = 1. / (1 + dAmplitude);
       
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

    dInvDeltaRow = - dIntensity * m_pdInverseDeltaX[nIdxI - 1];
    
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
        dTmp = dInvDeltaRow * dInvDeltaS0 * dCoe;

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
      dTmp = dInvDeltaRow * dInvDeltaColumn;
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
      dTmp = dInvDeltaRow * dInvDeltaColumn;
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
        dTmp = dInvDeltaRow * dInvDeltaSN2 * dCoe;
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
        dTmp = dInvDeltaRow * dInvDeltaColumn;
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
        dTmp = dInvDeltaRow * dInvDeltaColumn;
        sparseMatrix(nIdxRow, nIdxColumn)
            -= dTmp * 0.5 * ( Sqr(dMax - dXIm) - Sqr(dMin - dXIm) );
      }
    }


    // right support of the row base function
    // the base function is (dX - dXIp) / (dXI - dXIp)
    dXI = m_pdX[nIdxI];
    dXIp = m_pdX[nIdxI + 1];

    dInvDeltaRow = - dIntensity * m_pdInverseDeltaX[nIdxI];
    
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
        dTmp = dInvDeltaRow * dInvDeltaS0 * dCoe;

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
      dTmp = dInvDeltaRow * dInvDeltaColumn;
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
      dTmp = dInvDeltaRow * dInvDeltaColumn;
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
        dTmp = dInvDeltaRow * dInvDeltaSN2 * dCoe;
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
        dTmp = dInvDeltaRow * dInvDeltaColumn;
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
        dTmp = dInvDeltaRow * dInvDeltaColumn;
        sparseMatrix(nIdxRow, nIdxColumn)
            += dTmp * 0.5 * ( Sqr(dMax - dXIp) - Sqr(dMin - dXIp) );
      }
    }
  } // loop over space nodes
}

// A row or stupidest implementation of the jump system for sensitivity.
// It may be merged into the building of the main jump system
AutoPtr<MorseMatrix>
ForwardOptionStepperFE::BuildAmplitudeSensitivityJumpSystem
(size_t nIdxR1, size_t nIdxR2, double dAmplitude, double dIntensity)
{
  Array< std::list<size_t> > ppnColumnLists(m_nNbX);

  BuildSensitivityJumpSystemStructure
  ( nIdxR1, nIdxR2, dAmplitude, ppnColumnLists.Get() );

  shared_ptr<MorseStruct> 
    pMS( new MorseStruct(ppnColumnLists.Get(), m_nNbX) );

  AutoPtr<MorseMatrix> pSparseMatrix(new MorseMatrix);
  MorseMatrix& sparseMatrix = *pSparseMatrix; 

  sparseMatrix.Init(pMS, 0);

  BuildPartialJumpSystem
  ( nIdxR1, nIdxR2, dAmplitude, dIntensity, sparseMatrix );

  BuildAmplitudeSensitivityJumpSystem
  ( nIdxR1, nIdxR2, dAmplitude, dIntensity, sparseMatrix );

  return pSparseMatrix;
}


} // namespace hg

} // namespace ito33
