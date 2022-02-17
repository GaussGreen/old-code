/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/underlyingprocess.cpp
// Purpose:     Implementation of homogeneous underlying process class
// Created:     2005/04/15
// RCS-ID:      $Id: underlyingprocess.cpp,v 1.12 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/numeric/bisecnewton.h"
#include "ito33/numeric/densematrix.h"
#include "ito33/numeric/computeexpmatrix.h"

#include "ito33/hg/error.h"
#include "ito33/hg/underlyingprocess.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"

#include "hg/xml/underlyingprocess.h"

extern const ito33::hg::Error ITO33_HG_REGIME,
                              ITO33_HG_REGIMEINDEX, ITO33_HG_REGIMESIZE,
                              ITO33_HG_INTENSITY;

namespace ito33
{

namespace hg
{
  

class RealUnderlyingProcessCalculator
{
public:

  RealUnderlyingProcessCalculator(double dSR) : m_dSR(dSR) { }

  void SetParams(size_t nNb, double* pdAmplitudes, double *pdIntensities, 
                 double dSquareVol)
  {
    m_dSquareVol = dSquareVol;
    m_nNbJumps = nNb;
    m_pdAmplitudes = pdAmplitudes;
    m_pdIntensities = pdIntensities;
  }

  void operator()(double dTotalVol, double& dF, double& dDeriv)
  {
    dF = m_dSquareVol - dTotalVol * dTotalVol;
    dDeriv = - 2. * dTotalVol;

    for (size_t nIdx = 0; nIdx < m_nNbJumps; nIdx++)
    {
      double dIntensity = m_pdIntensities[nIdx];

      if (dIntensity != 0.0)
      {
        double dAmpltitude = m_pdAmplitudes[nIdx];
        double dTmp = 1. / (dTotalVol - m_dSR * m_pdAmplitudes[nIdx]);
        double dTmp1 = dIntensity * dAmpltitude * dAmpltitude;

        dF += dTotalVol * dTmp1 * dTmp;
        dDeriv += dTmp1 * dTmp - dTotalVol * dTmp1 * dTmp * dTmp;
      }
    }
  }

private:

  size_t m_nNbJumps;
  double m_dSquareVol;
  double m_dSR;
  double* m_pdAmplitudes;
  double* m_pdIntensities;
};

UnderlyingProcess::UnderlyingProcess(size_t nNbRegimes, 
                                     const std::vector<double>& pdVols,
                                     const std::vector<double>& pdIntensities)
                                   : m_nNbRegimes(nNbRegimes)
{
  CHECK_COND(m_nNbRegimes > 0, ITO33_HG_REGIME);

  SetVolatilities(pdVols);

  SetJumpsToDefault(pdIntensities);
}

void UnderlyingProcess::CheckRegimeIndex(size_t nIdxR) const
{
  CHECK_COND(nIdxR < m_nNbRegimes, ITO33_HG_REGIMEINDEX);
}

void UnderlyingProcess::CheckRegimeSize(size_t nNbRegimes) const
{
  CHECK_COND(nNbRegimes == m_nNbRegimes, ITO33_HG_REGIMESIZE);
}

void 
UnderlyingProcess::SetJumps(size_t nIdxR1, size_t nIdxR2, 
                           const std::vector<double>& pdIntensities,
                           const std::vector<double>& pdAmplitudes)
{
  Jumps& jumps = m_ppJumps[nIdxR1][nIdxR2];
  
  jumps.clear();
  for (size_t nIdx = 0; nIdx < pdIntensities.size(); nIdx++)
    jumps.push_back( Jump(pdIntensities[nIdx], pdAmplitudes[nIdx]) );
}

void 
UnderlyingProcess::SetJumps(size_t nIdxR1, size_t nIdxR2, const Jumps& jumps)
{
  CheckRegimeIndex(nIdxR1);
  CheckRegimeIndex(nIdxR2);

  m_ppJumps[nIdxR1][nIdxR2] = jumps;
}

void 
UnderlyingProcess::SetJumpsToDefault(const std::vector<double>& pdIntensities)
{
  /*
  for (size_t nIdx = 0; nIdx < pdIntensities.size(); nIdx++)
    CHECK_COND(pdIntensities[nIdx] >= 0, ITO33_HG_INTENSITY);
  */

  CheckRegimeSize( pdIntensities.size() ); 
  
  m_pdDefaultIntensities = pdIntensities;
}

void 
UnderlyingProcess::SetVolatilities(const std::vector<double>& pdVols)
{
  CheckRegimeSize( pdVols.size() );
  
  m_pdVols = pdVols;
}

std::vector<double> 
UnderlyingProcess::ComputeTotalVolatilities() const
{
  std::vector<double> pdTotalVols(m_nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    double dTmp;

    dTmp = m_pdVols[nIdxR1] * m_pdVols[nIdxR1]
         + m_pdDefaultIntensities[nIdxR1];

    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_ppJumps[nIdxR1][nIdxR2];
      Jumps::const_iterator pJump;

      for (pJump = jumps.begin(); pJump < jumps.end(); ++pJump)
      {
        dTmp += pJump->GetIntensity() 
              * pJump->GetAmplitude() * pJump->GetAmplitude();
      }
    }

    pdTotalVols[nIdxR1] = sqrt(dTmp);
  }

  return pdTotalVols;
}

size_t UnderlyingProcess::GetNbParameters() const
{
  size_t nNbParams = 0;

  // Account for volatilities
  nNbParams += m_nNbRegimes;
  
  // Account for default intensities
  nNbParams += m_nNbRegimes;

  // Add the jumps between regimes (amplitudes and intensities)
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      nNbParams += 2 * m_ppJumps[nIdxR1][nIdxR2].size();

  // The post default volatility
  nNbParams++;

  return nNbParams;
}

shared_ptr<UnderlyingProcess> 
UnderlyingProcess::ComputeUnderlyingProcess(double dSR) const
{
  shared_ptr<UnderlyingProcess> pUPTmp( new UnderlyingProcess(*this) );

  RealUnderlyingProcessCalculator calculator(dSR);

  std::vector<double> pdIntensities;
  std::vector<double> pdAmplitudes;

  std::vector<double> pdDedefaultIntensitiesTmp(m_nNbRegimes);

  size_t nNbJumps;

  size_t nIdxR1, nIdxR2;

  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    // Count the number of jumps from nIdxR1
    nNbJumps = 0;
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      nNbJumps += m_ppJumps[nIdxR1][nIdxR2].size();
    
    // Add Jump to default
    nNbJumps++;

    // The ordered jumps
    pdIntensities.resize(nNbJumps);
    pdAmplitudes.resize(nNbJumps);
    
    double dSquareVol = m_pdVols[nIdxR1] * m_pdVols[nIdxR1];
   
    size_t nIdx = 0;
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_ppJumps[nIdxR1][nIdxR2];
      Jumps::const_iterator pJump;

      for (pJump = jumps.begin(); pJump < jumps.end(); ++pJump, nIdx++)
      {
        pdIntensities[nIdx] = pJump->GetIntensity();
        pdAmplitudes[nIdx] = pJump->GetAmplitude();
      }
    }
    
    // Jump to default
    pdIntensities[nIdx] = m_pdDefaultIntensities[nIdxR1];
    pdAmplitudes[nIdx] = - 1;

    // Temporary viariables help to determine the bound of the total vol 
    double dAmplitudeMax = 0.;
    double dIntensity1 = 0.;
    for (nIdx = 0; nIdx < nNbJumps; nIdx++)
      if (pdIntensities[nIdx] != 0.0)
        if (pdAmplitudes[nIdx] > dAmplitudeMax)
        {
          dAmplitudeMax = pdAmplitudes[nIdx];
          dIntensity1 = pdIntensities[nIdx];
        }
    
    
    // The lower bound
    double dLower;

    if (dAmplitudeMax != 0.)
      dLower = ( 0.5 * dSR + sqrt(0.25 * dSR * dSR + dIntensity1) )
             * dAmplitudeMax;
    else
      dLower = 0.;

    if (dLower < fabs(m_pdVols[nIdxR1]))
      dLower = fabs(m_pdVols[nIdxR1]);
    
    // The upper bound
    double dUpper = m_pdVols[nIdxR1] * m_pdVols[nIdxR1];
    for (nIdx = 0; nIdx < nNbJumps; nIdx++)
      dUpper += pdIntensities[nIdx] * pdAmplitudes[nIdx] * pdAmplitudes[nIdx];

    dUpper = sqrt(dUpper) + dSR * dAmplitudeMax;

    // solve the equation on totoal vol 
    calculator.SetParams( nNbJumps, &pdAmplitudes[0], &pdIntensities[0], 
                          dSquareVol );

    numeric::BisecNewton solver(1.e-10);

    double dTotalVol = solver(calculator, dLower, dUpper); 

    // Compute the intensities of the real underlying process
    nIdx = 0;
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_ppJumps[nIdxR1][nIdxR2];
      
      Jumps::const_iterator pJump;
      
      Jumps jumpsTmp;
      for (pJump = jumps.begin(); pJump < jumps.end(); ++pJump, nIdx++)
      {
        double dIntensity;

        dIntensity = dTotalVol * pdIntensities[nIdx]
                   / (dTotalVol - dSR * pdAmplitudes[nIdx]);

        jumpsTmp.push_back( Jump(dIntensity, pdAmplitudes[nIdx]) );
      }

      pUPTmp->SetJumps(nIdxR1, nIdxR2, jumpsTmp);
    }
    
    pdDedefaultIntensitiesTmp[nIdxR1] = dTotalVol * pdIntensities[nIdx] 
                                      / (dTotalVol - dSR * pdAmplitudes[nIdx]);

    nIdx++;
  }

  // The end result
  shared_ptr<UnderlyingProcess> 
    pUP( new UnderlyingProcess
             ( m_nNbRegimes, m_pdVols, pdDedefaultIntensitiesTmp ) );

  for (nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    for (nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      pUP->SetJumps( nIdxR1, nIdxR2, pUPTmp->GetJumps(nIdxR1, nIdxR2) );

  return pUP;
}

numeric::DenseMatrix
UnderlyingProcess::ComputeRegimeTransitionProba(double dT) const
{
  // Form the matrix
  numeric::DenseMatrix matrix(m_nNbRegimes, m_nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    matrix[nIdxR1][nIdxR1] = dT * (- m_pdDefaultIntensities[nIdxR1]);
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      if (nIdxR2 != nIdxR1)
      {
        const Jumps& jumps( GetJumps(nIdxR1, nIdxR2) );
        Jumps::const_iterator pJump;

        matrix[nIdxR1][nIdxR2] = 0;
        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          matrix[nIdxR1][nIdxR1] -= dT * pJump->GetIntensity();
          matrix[nIdxR1][nIdxR2] += dT * pJump->GetIntensity();
        }
      }
    }
  }

  numeric::DenseMatrix expMatrix(m_nNbRegimes, m_nNbRegimes + 1);
  numeric::ComputeExpMatrix(m_nNbRegimes, matrix.Get(), expMatrix.Get());
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    double dTmp = 1;
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      dTmp -= expMatrix[nIdxR1][nIdxR2];
 
    expMatrix[nIdxR1][m_nNbRegimes] = dTmp;
  }

  return expMatrix;
}

numeric::DenseMatrix
UnderlyingProcess::ComputeSpotRegimeTransitionProba(double dT) const
{
  // Form the matrix
  numeric::DenseMatrix matrix(m_nNbRegimes, m_nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    matrix[nIdxR1][nIdxR1] = 0;
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      if (nIdxR2 != nIdxR1)
      {
        const Jumps& jumps( GetJumps(nIdxR1, nIdxR2) );
        Jumps::const_iterator pJump;

        matrix[nIdxR1][nIdxR2] = 0;
        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          double dTmp = pJump->GetIntensity() * (1 + pJump->GetAmplitude() );
          matrix[nIdxR1][nIdxR1] -= dT * dTmp;
          matrix[nIdxR1][nIdxR2] += dT * dTmp;
        }
      }
    }
  }

  numeric::DenseMatrix expMatrix(m_nNbRegimes, m_nNbRegimes + 1);
  numeric::ComputeExpMatrix(m_nNbRegimes, matrix.Get(), expMatrix.Get());
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    double dTmp = 1;
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
      dTmp -= expMatrix[nIdxR1][nIdxR2];
 
    expMatrix[nIdxR1][m_nNbRegimes] = dTmp;
  }

  return expMatrix;
}

void UnderlyingProcess::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_UNDERLYING_PROCESS, tagParent);
    
  // dump the number of regimes
  tag.Element(XML_TAG_REGIME_NUMBER)(m_nNbRegimes);

  // dump the volatilities
  DumpVector( tag, m_pdVols,
              XML_TAG_VOLATILITIES, XML_TAG_VOLATILITY);

  // dump the default intensities
  DumpVector( tag, m_pdDefaultIntensities,
              XML_TAG_DEFAULTINTENSITIES, XML_TAG_DEFAULTINTENSITY);

  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {   
      const Jumps& jumps = m_ppJumps[nIdxR1][nIdxR2];

      if ( !jumps.empty() )
      {
        // dump the internal jumps
        ito33::XML::Tag tagJumps(XML_TAG_JUMPS, tag);
        tagJumps.Element(XML_TAG_JUMP_FROM)(nIdxR1);
        tagJumps.Element(XML_TAG_JUMP_TO)(nIdxR2);

        Jumps::const_iterator pJump;

        for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        {
          ito33::XML::Tag tagJump(XML_TAG_JUMP, tagJumps);

          tagJump.Element(XML_TAG_INTENSITY)( pJump->GetIntensity() );
          tagJump.Element(XML_TAG_AMPLITUDE)( pJump->GetAmplitude() );
        }
      }
    }

  // Dump the data of finance::UnderlyingProcess
  DumpMe(tag);
}     


} // namespace hg

} // namespace ito33
