/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/backwardinstdata.h
// Purpose:     implementation for HG instdata class for backward solving
// Created:     2005/01/13
// RCS-ID:      $Id: backwardinstdata.cpp,v 1.22 2006/07/31 16:10:10 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/dividendevent.h"
#include "ito33/pricing/params.h"

#include "ito33/numeric/boundary1d.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/finance/computationalflags.h"

#include "hg/model.h"
#include "hg/sensitivitydata.h"
#include "hg/sensitivitymethod.h"
#include "hg/backwardinstdata.h"

namespace ito33
{

namespace hg
{

  
void BackwardInstData::Alloc(size_t nNbX)
{
  InstData::Alloc(nNbX);

  // Allocate the fugit arrays
  if (m_bComputeFugit)
  {
    m_pdFugits = Array<double>(nNbX);
    m_pdOldFugits = Array<double>(nNbX);
    m_pdOldOldFugits = Array<double>(nNbX);
  }

  // Clear the sensitivity datas
  m_pSensitivityData.clear();

  // No sensitivity is required
  if ( m_pbComputeSensitivities.empty() )
    return;

  // loop over the parameters. Note that this should use the same order as
  // in Translator
  size_t nCounterF = 0;

  // Sensitivity on volatilities
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    if ( m_pbComputeSensitivities[nCounterF++] )
    {
      SensitivityData data;
      data.m_sensitivityType = SensitivityType_Volatility;
      data.m_nRegime1 = nIdxR;

      m_pSensitivityData.push_back(data);
    }
  }

  // Sensitivity on default intensities
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    if ( m_pbComputeSensitivities[nCounterF++] )
    {
      SensitivityData data;
      data.m_sensitivityType = SensitivityType_DefaultIntensity;
      data.m_nRegime1 = nIdxR;

      m_pSensitivityData.push_back(data);
    }
  }

  // Sensitivities on no default jumps
  for (size_t nIdxR1 = 0; nIdxR1 < m_nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < m_nNbRegimes; nIdxR2++)
    {
      Jumps jumps = m_model.GetJumps(nIdxR1, nIdxR2);

      Jumps::const_iterator iterJumps;
      
      for (iterJumps = jumps.begin(); iterJumps != jumps.end(); ++iterJumps)
      {
        if ( m_pbComputeSensitivities[nCounterF++] )
        {
          SensitivityData data2;
          data2.m_sensitivityType = SensitivityType_JumpIntensity;
          data2.m_nRegime1 = nIdxR1;
          data2.m_nRegime2 = nIdxR2;
          data2.m_dAmplitude = iterJumps->GetAmplitude();
          data2.m_dIntensity = iterJumps->GetIntensity();

          m_pSensitivityData.push_back(data2);
        }

        if ( m_pbComputeSensitivities[nCounterF++] )
        {
          SensitivityData data1;
          data1.m_sensitivityType = SensitivityType_JumpAmplitude;
          data1.m_nRegime1 = nIdxR1;
          data1.m_nRegime2 = nIdxR2;
          data1.m_dAmplitude = iterJumps->GetAmplitude();
          data1.m_dIntensity = iterJumps->GetIntensity();

          m_pSensitivityData.push_back(data1);
        }
      }

    } // loop over regime 2
  } // loop over regime 1

  // Compute sensitivity by PDE requires additional arrays
  if ( m_sensitivityMethod == SensitivityMethod_PDE )
  {
    // Merge into the sensitivity datas?
    size_t nNbSensitivities = m_pSensitivityData.size();

    // allocate the storage for the solution arrays
    m_ppdSensitivities.resize(nNbSensitivities);
    m_ppdOldSensitivities.resize(nNbSensitivities);
    m_ppdOldOldSensitivities.resize(nNbSensitivities);
    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
    {
      m_ppdSensitivities[nIdxD] = Array<double>(nNbX);
      m_ppdOldSensitivities[nIdxD] = Array<double>(nNbX);
      m_ppdOldOldSensitivities[nIdxD] = Array<double>(nNbX);
    }
  }
}

void BackwardInstData::SetInitialSensitivityValue()
{
  m_aData.m_bIsValid = false;

  if ( m_sensitivityMethod == SensitivityMethod_PDE )
    for (size_t nIdxD = 0; nIdxD < m_pSensitivityData.size(); nIdxD++)
    {
      for (size_t nIdxS = 0; nIdxS < m_nNbX; nIdxS++)
      {
        m_ppdSensitivities[nIdxD][nIdxS] = 0.;
        m_ppdOldSensitivities[nIdxD][nIdxS] = 0.;
        m_ppdOldOldSensitivities[nIdxD][nIdxS] = 0.;
      }
    }

  if ( m_sensitivityMethod == SensitivityMethod_Adjoint )
    m_bDualSystemRequired = true;
}

void BackwardInstData::Swap()
{
  InstData::Swap();

  if (m_bComputeFugit)
  {
    swap(m_pdOldOldFugits, m_pdFugits);
    swap(m_pdOldOldFugits, m_pdOldFugits);
  }

  if ( m_sensitivityMethod == SensitivityMethod_PDE )
  {
    for (size_t nIdxD = 0; nIdxD < m_pSensitivityData.size(); nIdxD++)
    {
      swap(m_ppdOldOldSensitivities[nIdxD], m_ppdSensitivities[nIdxD]);
      swap(m_ppdOldOldSensitivities[nIdxD], m_ppdOldSensitivities[nIdxD]);
    }
  }
}

void BackwardInstData::ApplyEvent(const pricing::Event *pEvent)
{
  InstData::ApplyEvent(pEvent); 

  if (m_bComputeFugit)
  {
    double* pdFugits = m_pdFugits.Get();

    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      pEvent->ApplyToGreek(m_pdS, pdFugits + nIdxR * m_nNbS, m_nNbS);
  }

  // Don't store the current data for sensitivity by adjoint, it is not 
  // solution of a non trivial system.
  if ( m_sensitivityMethod == SensitivityMethod_Adjoint )
    m_aData.m_bIsValid = false;

  if ( m_sensitivityMethod == SensitivityMethod_PDE )
    for (size_t nIdxD = 0; nIdxD < m_pSensitivityData.size(); nIdxD++)
    {
      double* pdSensitivities = m_ppdSensitivities[nIdxD].Get();

      for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
        pEvent->ApplyToGreek(m_pdS, pdSensitivities + nIdxR * m_nNbS, m_nNbS);
    }
  
  if ( m_bDualSystemRequired )
  {
    if ( pEvent->GetType() == pricing::ET_Dividend )
    {
      const pricing::DividendEvent* 
        pDividendEvent = (pricing::DividendEvent *)pEvent;

      m_aData.m_pInterpMatrix = AutoPtr<numeric::InterpolationMatrix>
          ( pDividendEvent->GetInterpolationMatrix
                            (m_pdS, m_nNbS, m_nNbRegimes) );
    }
  }
}

double BackwardInstData::GetInitialSpot() const
{
  return m_params.GetSpotSharePrice(); 
}

void BackwardInstData::SetSensitivityBoundary()
{
  if (m_BoundaryCondition.GetLeftType() == numeric::BCType_Dirichlet)
    m_BoundaryCondition.SetLeft(numeric::BCType_Dirichlet, 0.0);

  if (m_BoundaryCondition.GetRightType() == numeric::BCType_Dirichlet)
    m_BoundaryCondition.SetRight(numeric::BCType_Dirichlet, 0.0);
}

void BackwardInstData::SetupFlags(const finance::ComputationalFlags& flags)
{
  // The usual Greek flags
  m_bComputeFugit = flags.GetComputeFugit();

  // If the flag for all the senstivities is set in the computational flags, compute
  // all sensitivities. If not set, then only compute the sensitivities
  // specified in the SensitivityFlag array
  if ( flags.AreAllSensitivitiesActivated() )
    m_pbComputeSensitivities.resize(m_model.GetNbParameters(), true);
  else
    m_pbComputeSensitivities = flags.GetSensitivityFlags();
  
  size_t nNbSensitivities = m_pbComputeSensitivities.size();

  bool bComputeSensitivity = false;
  // Don't count the sensitivity flag on post default volatility since it
  // will be always 0
  if ( nNbSensitivities )
    for (size_t nIdxD = 0; nIdxD < nNbSensitivities - 1; nIdxD++)
      if ( m_pbComputeSensitivities[nIdxD] )
      {
        bComputeSensitivity = true;
        break;
      }
  
  m_sensitivityMethod = SensitivityMethod_None;

  // Default to adjoint if not set to PDE in flags
  if ( bComputeSensitivity )
    m_sensitivityMethod = SensitivityMethod_Adjoint;

  if ( bComputeSensitivity && flags.GetSensitivityMethod() == 1 )
    m_sensitivityMethod = SensitivityMethod_PDE;
}

} // namespace hg

} // namespace ito33
