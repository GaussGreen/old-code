#include "ito33/beforestd.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/arrayutils.h"

#include "ito33/numeric/schemetype.h"

#include "ihg/optionstepperlogsgrid.h"

using namespace ito33;
using namespace ito33::numeric;
using ito33::ihg::OptionStepperLogSGrid;

void OptionStepperLogSGrid::Init()
{
  // Get the grid info
  m_nNbX = m_meshes.GetNbS();
  m_pdX = m_meshes.GetSMesh()->getLogS();

  // Need the real S grid to get hazard rates
  m_pdS = m_meshes.GetSMesh()->getS();

  OptionStepper::Init();
}


void OptionStepperLogSGrid::MakeCoefficients()
{
  size_t nIdx;

  // Get the various interest rates. These usually change at each timestep
  double dRate = m_instdata.m_dRate;
  double dForeignRate = m_instdata.m_dForeignRate;

  // Get the default recovery value
  double dRecovery = m_params.GetDefaultValue(m_instdata.m_dTime);

  // Get the volatility structure, and determine if it is flat
  shared_ptr<Volatility> pVol = m_params.GetVolatility();
  double dVolValue;
  bool bIsVolFlat = pVol->GetFlatVolatility(dVolValue);

  // Get the hazard rates
  for (nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdHazardRates[nIdx] = 0.0;
  m_params.GetHazardRate()->GetArray(m_instdata.m_dTime, m_pdS, 
    m_pdHazardRates.Get(), m_nNbX);

  if (bIsVolFlat)
  {
    double dTmp = 0.5 * dVolValue * dVolValue;
    for (nIdx = 0; nIdx < m_nNbX; nIdx++)
    {
      m_pdCoe2nd[nIdx] = dTmp;
      m_pdCoe1st[nIdx] = dRate - dForeignRate + m_pdHazardRates[nIdx] - dTmp;
      m_pdCoeZero[nIdx] = dRate + m_pdHazardRates[nIdx];
      m_pdCoeConst[nIdx] = m_pdHazardRates[nIdx] * dRecovery;
    }
  }
  else
  {
    ASSERT_MSG(false, "Vol surfaces are not yet supported in the stepper");
  }
}
