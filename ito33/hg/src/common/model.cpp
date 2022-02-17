/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/model.cpp
// Purpose:     implement HG pricing Model class
// Created:     2005/01/13
// RCS-ID:      $Id: model.cpp,v 1.12 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/densematrix.h"

#include "hg/model.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::Model);
} // namespace ito33

namespace ito33
{

namespace hg
{


Model::Model(const UnderlyingProcess& underlyingProcess,
             const shared_ptr<UnderlyingProcess>& pUnderlyingProcessForMesh) 
           : m_pUnderlyingProcessForMesh(pUnderlyingProcessForMesh),
             m_underlyingProcess(underlyingProcess)
{
  m_dPostDefaultVolatility = m_underlyingProcess.GetPostDefaultVolatility();
}

double Model::GetSquaredTotalVol(double /*dMaturity*/, 
                                 double /*dSpot*/, 
                                 bool bForMesh) const
{

  // Find the max total volatility in the regimes.  If we wanted a sum,
  // or average, it would be easier to call ComputeTotalVolatilities().
  // For now, duplicate the code

  // The min value to be returned
  double dTmp = 0.01;

  // use m_pUnderlyingProcessForMesh for data if vol for mesh is requested
  const UnderlyingProcess* pProcess;
  if (bForMesh)
  {
    ASSERT_MSG( m_pUnderlyingProcessForMesh, 
      "Invalid process for mesh in Model::GetSquaredTotalVol" );

    pProcess = m_pUnderlyingProcessForMesh.get();
  }
  else
    pProcess = &m_underlyingProcess;

  size_t nNbRegimes = pProcess->GetNbRegimes();
  const std::vector<double>& pdVols = pProcess->GetVolatilities();
  const std::vector<double>& pdDefault = pProcess->GetJumpsToDefault();

  for (size_t nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    double dSQRTotalVol = pdVols[nIdxR] * pdVols[nIdxR];

    dSQRTotalVol += pdDefault[nIdxR];
      
    for (size_t nIdxR1 = 0; nIdxR1 < GetNbRegimes(); nIdxR1++)
    {
      const Jumps& jumps = pProcess->GetJumps(nIdxR, nIdxR1);
      Jumps::const_iterator pJump;
      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        dSQRTotalVol += pJump->GetIntensity() 
                      * pJump->GetAmplitude() * pJump->GetAmplitude();
    }

    // Find the max regime value
    if (dTmp < dSQRTotalVol)
      dTmp = dSQRTotalVol;
  }
  
  return dTmp;
}

double 
Model::GetConvection(double dMaturity, double dSpot) const
{
  double dHR = 0.;
  
  return pricing::Model::GetSquaredTotalVol(dMaturity, dSpot) - dHR * 2;
}

numeric::DenseMatrix Model::ComputeRegimeTransitionProba(double dT) const
{
  return m_underlyingProcess.ComputeRegimeTransitionProba(dT);
}

numeric::DenseMatrix Model::ComputeSpotRegimeTransitionProba(double dT) const
{
  return m_underlyingProcess.ComputeSpotRegimeTransitionProba(dT);
}

} // namespace hg

} // namespace ito33
