/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/translator.cpp
// Purpose:     implementation of HG helper class Translator for calibration
// RCS-ID:      $Id: translator.cpp,v 1.10 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/autoptr.h"

#include "hg/translator.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(hg::Translator);

} // namespace ito33

namespace ito33
{

namespace hg
{


Translator::Translator(const UnderlyingProcess& underlyingProcess)
: m_underlyingProcess(underlyingProcess)
{
  size_t nNbRegimes = m_underlyingProcess.GetNbRegimes();

  m_nNbFlags = 0;

  // Account for volatilities
  m_nNbFlags += nNbRegimes;
  
  // Account for default intensities
  m_nNbFlags += nNbRegimes;

  // Add the jumps between regimes (amplitudes and intensities)
  size_t nIdxR1, nIdxR2;
  for (nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    for (nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
      m_nNbFlags += 2 * m_underlyingProcess.GetJumps(nIdxR1, nIdxR2).size();

  // Post default volatility
  m_nNbFlags++;

  m_flags.resize(m_nNbFlags);
  
  for (size_t n = 0; n < m_nNbFlags; n++)
    m_flags[n] = true;
}


shared_ptr<UnderlyingProcess> 
Translator::operator()(const double* pdXOriginal)
{
  size_t nIdxR;

  size_t nCounterP = 0;
  size_t nCounterF = 0;

  size_t nNbRegimes = m_underlyingProcess.GetNbRegimes();

  // Apply the bounds (note that the NAG optimizer can sometimes make a guess
  // that violates the bounds, especially the lower bounds)
  const std::vector<double> pdXNew = ApplyBounds(pdXOriginal);
  const double* pdX = &pdXNew[0];

  // Get the volatilities
  std::vector<double> pdProcessVols = m_underlyingProcess.GetVolatilities();
  std::vector<double> pdVols(nNbRegimes); 
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdVols[nIdxR] = pdX[nCounterP++];
    else
      pdVols[nIdxR] = pdProcessVols[nIdxR];
  }

  // The default intensities
  std::vector<double> pdProcessInt = m_underlyingProcess.GetJumpsToDefault();
  std::vector<double> pdDefaultIntensities(nNbRegimes);
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdDefaultIntensities[nIdxR] = pdX[nCounterP++];
    else
      pdDefaultIntensities[nIdxR] = pdProcessInt[nIdxR];
  }

  // Setup the underlying process
  shared_ptr<UnderlyingProcess> 
    pUnderlyingProcess( new UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );
  
  size_t nIdxR1, nIdxR2;
  for (nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    for (nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = m_underlyingProcess.GetJumps(nIdxR1, nIdxR2);
      Jumps::const_iterator pJump;
      Jumps jumpsTmp;

      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
      {
        double dIntensity;
        if ( m_flags[nCounterF++] )
          dIntensity = pdX[nCounterP++];
        else
          dIntensity = pJump->GetIntensity();
        
        double dAmplitude;
        
        if (m_flags[nCounterF++] )
          dAmplitude = pdX[nCounterP++];
        else
          dAmplitude = pJump->GetAmplitude();

        Jump jump(dIntensity, dAmplitude);

        jumpsTmp.push_back(jump);
    
      } // loop over jumps from regime i to j

      pUnderlyingProcess->SetJumps(nIdxR1, nIdxR2, jumpsTmp);

    } // inner regime loop

  } // outer regime loop

  // Post default volatility
  double dPostDefaultVol;
  if ( m_flags[nCounterF++] )
    dPostDefaultVol = pdX[nCounterP++];
  else
    dPostDefaultVol = m_underlyingProcess.GetPostDefaultVolatility();
  
  pUnderlyingProcess->SetPostDefaultVolatility(dPostDefaultVol);

  return pUnderlyingProcess;
}


shared_ptr<UnderlyingProcess> 
Translator::operator()(const std::vector<double>& pdX)
{
  return operator()(&pdX[0]);
}


void Translator::GetParameters(double* pdX) const
{
  size_t nIdxR;
 
  size_t nCounterF = 0;
  size_t nCounterP = 0;
 
  size_t nNbRegimes = m_underlyingProcess.GetNbRegimes();

  // volatilities
  std::vector<double> pdVols = m_underlyingProcess.GetVolatilities();
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdX[nCounterP++] = pdVols[nIdxR];
  }

  std::vector<double> pdDefaultIntensities = 
    m_underlyingProcess.GetJumpsToDefault();
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdX[nCounterP++] = pdDefaultIntensities[nIdxR];
  }

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      Jumps jumps = m_underlyingProcess.GetJumps(nIdxR1, nIdxR2);
      Jumps::const_iterator iter;
      for (iter = jumps.begin(); iter != jumps.end(); ++iter)
      {
        if ( m_flags[nCounterF++] )
          pdX[nCounterP++] = iter->GetIntensity();

        if ( m_flags[nCounterF++] )
          pdX[nCounterP++] = iter->GetAmplitude();
      } // iterate through the jumps from regime i to j

    } // inner regime loop

  } // outer regime loop

  if ( m_flags[nCounterF++] )
    pdX[nCounterP++] = m_underlyingProcess.GetPostDefaultVolatility();
}


std::vector<double> Translator::GetParameters() const
{
  std::vector<double> parameters( GetNbParameters() );
  
  GetParameters(&parameters[0]);

  return parameters;
}

std::vector<double> Translator::GetLowerBounds() const
{

  size_t nNbRegimes = m_underlyingProcess.GetNbRegimes();

  // Initialize the return vector
  size_t nNbParams = GetNbParameters();

  std::vector<double> pdLowerBounds(nNbParams);

  // loop through the params, checking for active status and type
  size_t nIdxR;
  
  size_t nCounterF = 0;
  size_t nCounterP = 0;
  
  // volatilities
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdLowerBounds[nCounterP++] = 0.0;
  }

  // default intensities
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdLowerBounds[nCounterP++] = 0.0;
  }

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      Jumps::const_iterator iter;
      Jumps jumps = m_underlyingProcess.GetJumps(nIdxR1, nIdxR2);
      for (iter = jumps.begin(); iter != jumps.end(); ++iter)
      {
        // intensity
        if ( m_flags[nCounterF++] )
          pdLowerBounds[nCounterP++] = 0.0;

        // amplitude
        if ( m_flags[nCounterF++] )
          pdLowerBounds[nCounterP++] = -0.99;
      } // iterate through the jumps from regime i to j

    } // inner regime loop

  } // outer regime loop

  // post default volatility
  if ( m_flags[nCounterF++] )
    pdLowerBounds[nCounterP++] = 0.0;

  return pdLowerBounds;

}


std::vector<double> Translator::GetUpperBounds() const
{

  size_t nNbRegimes = m_underlyingProcess.GetNbRegimes();

  // Initialize the return vector
  size_t nNbParams = GetNbParameters();

  std::vector<double> pdUpperBounds(nNbParams);

  // loop through the params, checking for active status and type
  size_t nIdxR;
  
  size_t nCounterF = 0;
  size_t nCounterP = 0;
  
  // volatilities
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdUpperBounds[nCounterP++] = 5.0;
  }

  // default intensities
  for (nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    if ( m_flags[nCounterF++] )
      pdUpperBounds[nCounterP++] = 10.0;
  }

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      Jumps jumps = m_underlyingProcess.GetJumps(nIdxR1, nIdxR2);
      Jumps::const_iterator iter;
      for (iter = jumps.begin(); iter != jumps.end(); ++iter)
      {
        // intensity
        if ( m_flags[nCounterF++] )
          pdUpperBounds[nCounterP++] = 10.0;

        // amplitude
        if ( m_flags[nCounterF++] )
          pdUpperBounds[nCounterP++] = 1.0;

      } // iterate through the jumps from regime i to j

    } // inner regime loop

  } // outer regime loop

  // post default volatility
  if ( m_flags[nCounterF++] )
    pdUpperBounds[nCounterP++] = 5.0;

  return pdUpperBounds;
}


} // namespace hg

} // namespace ito33
