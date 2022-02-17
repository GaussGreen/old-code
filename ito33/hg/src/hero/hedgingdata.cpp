/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/hedgingdata.cpp
// Purpose:     implementation of the contract class for HERO
// Created:     2005/09/26
// RCS-ID:      $Id: hedgingdata.cpp,v 1.8 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/array.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivative.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain.h"

#include "ito33/numeric/qp_nag.h"
#include "ito33/numeric/densematrix.h"

#include "ito33/hg/underlyingprocess.h"
#include "ito33/hg/modeloutput.h"

#include "hg/backwardnumoutput.h"
#include "hg/hedgingdata.h"


namespace ito33
{

namespace hg
{

// Default to SVD for hero calculations
bool HedgingData::bForceNAG = false;


HedgingData::HedgingData(const finance::Derivative& derivative, 
                         shared_ptr<UnderlyingProcess> pRealProcess,
                         size_t nNbHedges)
: m_target(derivative), m_nNbHedges(nNbHedges), m_pRealUP(pRealProcess)
{
  // Base class data
  m_dMaturityTime = GetDoubleFrom( derivative.GetMaturityDate() );
  m_pDerivativeCurve = derivative.GetSessionData()->GetYieldCurve();

  // Iinitialize all the data rquired for hero and hedge calculations
  m_nNbRegimes = pRealProcess->GetNbRegimes();

  // deltas for the target and each hedge contract: [hedge][grid]
  m_ppdDeltas.resize(nNbHedges+1);

  // original prices and spots (ie. not interpolated onto the PDE grid)
  // for the target and hedge contracts: [hedge][original grid]
  m_ppdOrigSpots.resize(nNbHedges+1);
  m_ppdOrigPrices.resize(nNbHedges+1);

  // final price and delta at the spot for the target and hedges
  m_pdPrices.resize(nNbHedges+1);
  m_pdDeltas.resize(nNbHedges+1);

  // default values for the target and hedges
  m_pdValuesAfterDefault.resize(nNbHedges+1);  

  // inverses of the squared total vols
  m_pdInvSqrTotalVols.resize(m_nNbRegimes);    
  std::vector<double> pdTotalVols = m_pRealUP->ComputeTotalVolatilities();
  for (size_t nIdx = 0; nIdx < m_nNbRegimes; nIdx++)
    m_pdInvSqrTotalVols[nIdx] = 1.0 / (pdTotalVols[nIdx] * pdTotalVols[nIdx]);

  // Create vectors of amplitudes and intensities from each regime. Does
  // not matter where they jump to, so long as they are consistent with the
  // jump differences computed in functions below. Indexed by [regime][jump]
  m_ppdJumpIntensities.resize(m_nNbRegimes);
  m_ppdJumpAmplitudes.resize(m_nNbRegimes);

  const std::vector<double>& 
    pdDefaultIntensities = m_pRealUP->GetJumpsToDefault();

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {    
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_pRealUP->GetJumps(nIdxR1, nIdxR2);
     
      Jumps::const_iterator pJump;
      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
      {
        const double dIntensity = pJump->GetIntensity();
        const double dAmplitude = pJump->GetAmplitude();

        m_ppdJumpIntensities[nIdxR1].push_back(dIntensity);
        m_ppdJumpAmplitudes[nIdxR1].push_back(dAmplitude);

      } // loop over jumps   
    } // loop over regime 2

    // jump to default
    m_ppdJumpIntensities[nIdxR1].push_back(pdDefaultIntensities[nIdxR1]);
    m_ppdJumpAmplitudes[nIdxR1].push_back(-1.0);

  } // loop over regime 1

  // allocate the known sizes of the jump differences
  m_ppppdJumpDiffs.resize(nNbHedges+1);
  for (size_t nIdxH = 0; nIdxH < nNbHedges + 1; nIdxH++)
  {
    m_ppppdJumpDiffs[nIdxH].resize(m_nNbRegimes);
  
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      size_t nNbJumps = m_ppdJumpAmplitudes[nIdxR].size();

      m_ppppdJumpDiffs[nIdxH][nIdxR].resize(nNbJumps);
    }
  }

  // Both the hedge and hero functions store data in a matrix.
  // Allocate the space here
  m_matrix.Init( m_nNbHedges + 1, m_nNbHedges + 1 );

  // Both the hedge and hero functions need to solve a matrix.
  // Allocate the space here
  if ( m_nNbHedges )
    m_matrixToSolve.Init(m_nNbHedges, m_nNbHedges);

  // If solving by SVD, need W vector and V matrix.  U matrix will be
  // stored in ppdMatrixToSolve
  pdW.resize(m_nNbHedges);

  if ( m_nNbHedges )
    m_matrixV.Init(m_nNbHedges, m_nNbHedges);
}


void HedgingData::ComputeFinalHedgeRatios(double& dTargetRatio, 
                              std::vector<double>& pdHedgeRatios) const
{
  // Should be the correct size already, but resizing makes sure
  pdHedgeRatios.clear();
  pdHedgeRatios.resize(m_nNbHedges);
  
  // Get the price data at the analysis time, which is assumed to be
  // set to the valuation date
  for (size_t nIdxH = 0; nIdxH < m_nNbHedges + 1; nIdxH++)
  {
    shared_ptr< finance::ModelOutput > pMO;
    
    if (nIdxH == m_nNbHedges)
      pMO = m_pTargetOutput;
    else
      pMO = m_ppHedgeOutputs[nIdxH];
      
    shared_ptr<BackwardNumOutput>
      pNumOutput( static_pointer_cast<BackwardNumOutput>(pMO->GetNumOutput()) );

    m_ppdOrigSpots[nIdxH] = pNumOutput->GetSpotsAtAnalysisDate();
    m_ppdOrigPrices[nIdxH] = pNumOutput->GetPricesAtAnalysisDate();

    m_pdPrices[nIdxH] = pMO->GetPrice();
    m_pdDeltas[nIdxH] = pMO->GetDelta();

    m_pdValuesAfterDefault[nIdxH] = pMO->GetValueAfterDefault();
  }

  // Allocate the 'delta' or 'A' vector: delta (dF^x)
  Array<double> pdA(m_nNbHedges + 1);

  // Actually construct the data. Matrix allocated in constructor.
  ComputeHedgeData(m_matrix.Get(), pdA.Get());

  // hedges = V^(-1) C.  Solve using nag, so redundant contracts are zeroed.
  if (m_nNbHedges > 0)
    Solve(m_matrix.Get(), m_matrix[m_nNbHedges], 
          &pdHedgeRatios[0], true);

  // Compute the hedge on underlying 
  dTargetRatio = pdA[m_nNbHedges];
  for (size_t nIdxH1 = 0; nIdxH1 < m_nNbHedges; nIdxH1++)
    dTargetRatio -= pdHedgeRatios[nIdxH1] * pdA[nIdxH1];

}



void HedgingData::ComputeHedgeData(double** ppdVMatrix, 
                                   double* pdA) const
{

  // Final hedge data is specific to regime zero, and to the spot.
  // This is the primary difference between this function and the 
  // hero PDE function below.  It might be possible to avoid
  // some of the code duplication.  For now, don't worry about it.

  // Some helper arrays
  Array<double> pdVHedge(m_nNbHedges + 1);
  Array< std::vector<double> > ppdJumpPrices(m_nNbHedges + 1);

  // The final hedges are only needed for the first regime
  const size_t nIdxR1 = 0;

  const double dSpot = m_target.GetSessionData()->GetSpotSharePrice();

  const double dDefaultIntensity = m_pRealUP->GetJumpsToDefault()[nIdxR1];
  const double dVol = m_pRealUP->GetVolatilities()[nIdxR1];

  // include the target in the helper arrays
  for (size_t nIdxH = 0; nIdxH < m_nNbHedges + 1; nIdxH++)
  {
    // Get the appropriate data for this hedging instrument. Spot, price and 
    // delta vectors should have been set before calling this function.
    // Use the original grid to help with interpolation
    const size_t nNbS = m_ppdOrigSpots[nIdxH].size();

    const double* pdS = &m_ppdOrigSpots[nIdxH][0];

    pdVHedge[nIdxH] = dSpot * m_pdDeltas[nIdxH] * dVol;
    pdA[nIdxH] = pdVHedge[nIdxH] * dVol;

    double dJumpPrice;
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      // Get the prices for the regime we are jumping to.
      const double* pdPrices = &m_ppdOrigPrices[nIdxH][nIdxR2 * nNbS];
       
      // Jumps to non default regimes
      const Jumps& jumps = m_pRealUP->GetJumps(nIdxR1, nIdxR2);

      Jumps::const_iterator pJump;
      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
      {
        const double dIntensity = pJump->GetIntensity();
        const double dAmplitude = pJump->GetAmplitude();
        const double dS = dSpot * (1. + dAmplitude);
    
        // Contracts with Dirichlet nodes should have had extra points added
        // to the grid so that linear extrapolation will always work
        numeric::Interpolate(pdS, pdPrices, nNbS, &dS, &dJumpPrice, 1);
       
        // JumpPrice becomes  F(T, Spot*(1-amplitude), 0) - F(T, Spot, 0) 
        dJumpPrice -= m_pdPrices[nIdxH];
        
        ppdJumpPrices[nIdxH].push_back(dJumpPrice);

        pdA[nIdxH] += dIntensity * dAmplitude * dJumpPrice;
      }
    }

    // Jumps to default
    dJumpPrice = m_pdValuesAfterDefault[nIdxH] - m_pdPrices[nIdxH];

    ppdJumpPrices[nIdxH].push_back(dJumpPrice);

    pdA[nIdxH] += dDefaultIntensity * (- 1.) * dJumpPrice;
  }

  
  // Build the matrix. Include the target contract
  for (size_t nIdxH1 = 0; nIdxH1 < m_nNbHedges + 1; nIdxH1++)
  {
    for (size_t nIdxH2 = 0; nIdxH2 < m_nNbHedges + 1; nIdxH2++)
    {
      ppdVMatrix[nIdxH1][nIdxH2] = pdVHedge[nIdxH1] * pdVHedge[nIdxH2]
                                 - pdA[nIdxH1] * pdA[nIdxH2] 
                                   * m_pdInvSqrTotalVols[nIdxR1];
      
      size_t nIdxJump = 0;
      for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      {
        const Jumps& jumps = m_pRealUP->GetJumps(nIdxR1, nIdxR2);

        Jumps::const_iterator pJump;
        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump, nIdxJump++)
        {    
          const double dIntensity = pJump->GetIntensity();
          ppdVMatrix[nIdxH1][nIdxH2] += dIntensity 
                                      * ppdJumpPrices[nIdxH1][nIdxJump]
                                      * ppdJumpPrices[nIdxH2][nIdxJump];
        }
      }

      ppdVMatrix[nIdxH1][nIdxH2] += dDefaultIntensity 
                                  * ppdJumpPrices[nIdxH1][nIdxJump]
                                  * ppdJumpPrices[nIdxH2][nIdxJump];

    } // column loop (includes target contract)

  } // row loop (includes target contract)


  // vector 'A' was not fully constructed while being used to make 'V'
  for (size_t nIdxH1 = 0; nIdxH1 < m_nNbHedges + 1; nIdxH1++)
    pdA[nIdxH1] *= m_pdInvSqrTotalVols[nIdxR1] / dSpot;
 
}


std::vector<double> 
HedgingData::ComputePDETermsAt(double dTime, 
                               const double* pdSpots, 
                               size_t nNbS) const
{

  // surface interpolation needs a vector for the spots.  Copy.
  std::vector<double> pdVectorSpots(nNbS);
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdVectorSpots[nIdx] = pdSpots[nIdx];

  // For each regime, get where each jump takes the spots.  If the spots
  // are constant, this data can be precomputed (does not need to
  // be recomputed each timestep)
  std::vector< std::vector< std::vector<double> > > pppdNewSpots(m_nNbRegimes);
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    size_t nNbJumps = m_ppdJumpAmplitudes[nIdxR].size();
    pppdNewSpots[nIdxR].resize(nNbJumps);
    
    for (size_t nIdxJump = 0; nIdxJump < nNbJumps; nIdxJump++)
    {
      pppdNewSpots[nIdxR][nIdxJump].resize( nNbS );
      
      double dAmplitude = m_ppdJumpAmplitudes[nIdxR][nIdxJump];
      double dShift = (1.0 + dAmplitude);
      
      for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)        
        pppdNewSpots[nIdxR][nIdxJump][nIdxS] = pdSpots[nIdxS] * dShift;

    } // loop over jumps
  } // loop over regimes


  // Get the price diffs and deltas for each hedge contract at the given time.
  // The price surface accepts vectors only, not double pointers.
  std::vector<double> pdValuesTmp(nNbS);
  std::vector<double> pdJumpValuesTmp(nNbS);

  for (size_t nIdxH = 0; nIdxH < m_nNbHedges + 1; nIdxH++)
  {

    // Get the numerical output from the model output. Then get the surface.
    shared_ptr< finance::ModelOutput > pMO;
    
    if (nIdxH == m_nNbHedges)
      pMO = m_pTargetOutput;
    else
      pMO = m_ppHedgeOutputs[nIdxH];
      
    // resize the delta vector. Could also have extra indirection for regime
    m_ppdDeltas[nIdxH].resize(nNbS * m_nNbRegimes);

    shared_ptr<BackwardNumOutput>
      pNumOutput( static_pointer_cast<BackwardNumOutput>(pMO->GetNumOutput()) );

    // Get the recovery value at this time    
    const std::vector<double>& 
      pdTimes = pNumOutput->GetTimes(),
      pdRecoveryValues = pNumOutput->GetRecoveryValues();

    double dRecovery = 0.0;
    if ( numeric::IsEqualOrBefore(dTime, pdTimes[0]) )
    {

      size_t nIdxTime = 1;
      while ( pdTimes[nIdxTime] > dTime && nIdxTime < pdTimes.size() - 1 )
        nIdxTime++;

      // time should be between times[nIdxTime-1] and times[nIdxTime]
      double dTime1 = pdTimes[nIdxTime-1];
      double dTime2 = pdTimes[nIdxTime];
      double dWeight = (dTime - dTime1) / ( dTime2 - dTime1);

      double dV1 = pdRecoveryValues[nIdxTime-1];
      double dV2 = pdRecoveryValues[nIdxTime];

      dRecovery = dWeight * dV2 + (1. - dWeight) * dV1;
    }

    // Get the data for each regime
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {      
      shared_ptr<BackwardNumOutput>
        pNumOutput(static_pointer_cast<BackwardNumOutput>(pMO->GetNumOutput()));

      // The price surface for this contract and this regime
      shared_ptr<numeric::SurfaceGeneral> 
        pPriceSurface( pNumOutput->GetPriceSurface(nIdxR) );

      // Get the prices on the PDE grid    
      pPriceSurface->GetFirstValuesAt(dTime, pdVectorSpots, pdValuesTmp);

      // Get the deltas on the PDE grid    
      size_t nOffset = nIdxR * nNbS;
      numeric::ComputeDelta(&pdVectorSpots[0], &pdValuesTmp[0], nNbS, 
                            &m_ppdDeltas[nIdxH][nOffset]);

      // Loop over the jumps, computing the price differences
      // F(t, S(1+y), l) - F(t, S, k)
      size_t nNbJumps = m_ppdJumpAmplitudes[nIdxR].size();
    
      for (size_t nIdxJump = 0; nIdxJump < nNbJumps; nIdxJump++)
      {
        m_ppppdJumpDiffs[nIdxH][nIdxR][nIdxJump].resize( nNbS );

        // Handle default jump differently
        if ( nIdxJump == nNbJumps - 1 )
        {
          for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
            m_ppppdJumpDiffs[nIdxH][nIdxR][nIdxJump][nIdxS] = 
              dRecovery - pdValuesTmp[nIdxS];
        }
        else
        {
          pPriceSurface->GetFirstValuesAt
            (dTime, pppdNewSpots[nIdxR][nIdxJump], pdJumpValuesTmp);

          for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
            m_ppppdJumpDiffs[nIdxH][nIdxR][nIdxJump][nIdxS] = 
              pdJumpValuesTmp[nIdxS] - pdValuesTmp[nIdxS];
        }
      } // loop over jumps

    } // loop over regimes    

  } // loop over hedge contracts

  // Construct the data using the data computed above
  return ComputePDETerms(pdSpots, nNbS);

}



std::vector<double> 
HedgingData::ComputePDETerms(const double* pdSpots, size_t nNbS) const
{

  // The return vector
  std::vector<double> pdPDETerms(nNbS * m_nNbRegimes);

  // vols are constant.  Total vols computed in constructor.
  const std::vector<double>& pdVols = m_pRealUP->GetVolatilities();

  // loop over the regimes
  size_t nIdxR1;
  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {

    // vol is constant in a regime
    double dVol = pdVols[nIdxR1];      

    // Offset to this regime into arrays containing all regimes
    size_t nOffsetR1 = nIdxR1 * nNbS;

    // Build the matrices, solve, and update the return coefficient.    
    for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    {
           
      for (size_t nIdxH1 = 0; nIdxH1 < m_nNbHedges + 1; nIdxH1++)
      {

        double* pdDeltasH1 = &m_ppdDeltas[nIdxH1][nOffsetR1];

        double dTerm1 = pdSpots[nIdxS] * pdDeltasH1[nIdxS] * dVol;

        for (size_t nIdxH2 = 0; nIdxH2 < m_nNbHedges + 1; nIdxH2++)
        {
          double* pdDeltasH2 = &m_ppdDeltas[nIdxH2][nOffsetR1];

          double dTerm2 = pdSpots[nIdxS] * pdDeltasH2[nIdxS] * dVol;
        
          double dSum1 = 0.0;
          double dSum2 = 0.0;
          double dSumCombined = 0.0;

          size_t nNbJumps = m_ppdJumpAmplitudes[nIdxR1].size();
          for (size_t nIdxJump = 0; nIdxJump < nNbJumps; nIdxJump++)
          {
            dSum1 += m_ppdJumpIntensities[nIdxR1][nIdxJump] 
                     * m_ppppdJumpDiffs[nIdxH1][nIdxR1][nIdxJump][nIdxS]
                     * m_ppdJumpAmplitudes[nIdxR1][nIdxJump];

            dSum2 += m_ppdJumpIntensities[nIdxR1][nIdxJump] 
                     * m_ppppdJumpDiffs[nIdxH2][nIdxR1][nIdxJump][nIdxS]                    
                     * m_ppdJumpAmplitudes[nIdxR1][nIdxJump];

            dSumCombined += m_ppdJumpIntensities[nIdxR1][nIdxJump] 
                          * m_ppppdJumpDiffs[nIdxH1][nIdxR1][nIdxJump][nIdxS]
                          * m_ppppdJumpDiffs[nIdxH2][nIdxR1][nIdxJump][nIdxS];

          } // loop over jumps

          double dTmp1 = dTerm1*dTerm2 + dSumCombined;
          double dTmp2 = (dTerm1 * dVol + dSum1) * (dTerm2 * dVol + dSum2) 
                       * m_pdInvSqrTotalVols[nIdxR1];
          
          m_matrix[nIdxH1][nIdxH2] = dTmp1 - dTmp2;
    
        } // column loop (includes target contract)

      } // row loop (includes target contract)

      double dHeroTerm = m_matrix[m_nNbHedges][m_nNbHedges];

      if (m_nNbHedges > 0)
      {
        // solve using svd, since it should be quicker than nag
        std::vector<double> pdX(m_nNbHedges);
        Solve(m_matrix.Get(), m_matrix[m_nNbHedges], &pdX[0], false);

        // Wanted C V^-1 C.  Have V^-1 C from the solve
        for (size_t nIdx = 0; nIdx < m_nNbHedges; nIdx++)
          dHeroTerm -= pdX[nIdx] * m_matrix[m_nNbHedges][nIdx];
        
      } // if there is a hedge instrument      

      // Hero part should always be positive. Avoid round-off error.
      if (dHeroTerm < 0.0)
        dHeroTerm = 0.0;

      pdPDETerms[nIdxS + nOffsetR1] = dHeroTerm;

    } // loop over S      
 
  } // loop over regimes

  return pdPDETerms;
}


void 
HedgingData::Solve(const double* const* ppdA, 
                   const double* pdB, 
                   double* pdX,
                   bool bUseNag) const
{

  // Copy the data into a temp matrix for 2 reasons
  // 1) Some solvers need contiguous memory (at least they appear to)
  // 2) Some solvers over-write the matrix

  for (size_t n1 = 0; n1 < m_nNbHedges; n1++)
    for (size_t n2 = 0; n2 < m_nNbHedges; n2++)
      m_matrixToSolve[n1][n2] = ppdA[n1][n2];

  // bForceNAG can be set by test code
  if ( bUseNag || bForceNAG )
  {

    // Need to change the sign of the RHS vector V(dF^x, dF) since QP
    // solves 1/2 x^T V x + C^T x  instead of  1/2 x^T V x - C^T x
    Array<double> pdC(m_nNbHedges);
    for (size_t nIdxH = 0; nIdxH < m_nNbHedges; nIdxH++)
      pdC[nIdxH] = - pdB[nIdxH];
   
    numeric::QPMinimizerNAG minimizer(m_nNbHedges);
    minimizer(m_matrixToSolve.Get(), pdC.Get(), pdX);
  }
  else
  {

    // Do the singular value decomposition
    numeric::SVD(m_matrixToSolve.Get(), m_nNbHedges, m_nNbHedges, &pdW[0], 
                 m_matrixV.Get() );

    // By default, the SVD solver returns the "shortest vector" as the solution
    // when the matrix is singular.  This is not always what we want.
    // But since the singular are ordered from large to small, we don't
    // know which equations are redundant
    numeric::SVD_Solve(m_matrixToSolve.Get(), &pdW[0], m_matrixV.Get(), 
                       m_nNbHedges, m_nNbHedges, pdB, pdX );
  }


}


} // namespace hg

} // namespace ito33
