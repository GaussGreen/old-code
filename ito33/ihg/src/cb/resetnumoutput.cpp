/////////////////////////////////////////////////////////////////////////////
// Name:        cb/resetnumoutput.cpp
// Purpose:     implementation of ResetNumOutput class 
// Author:      David and Yann
// Created:     2004/11/08
// RCS-ID:      $Id: resetnumoutput.cpp,v 1.10 2006/05/01 20:26:10 dave Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/resetparams.h"
#include "ito33/pricing/cbmeshmanager.h"

#include "ito33/numeric/domain_general.h"
#include "ito33/numeric/predicatetime.h"

#include "ihg/resetnumoutput.h"
#include "ihg/backwardinstdata.h"

using namespace ito33;
using namespace ito33::numeric;
using namespace ito33::ihg;


// implement the AutoPtrDeleter for ResetNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::ResetNumOutput);
}

ResetNumOutput::ResetNumOutput(pricing::ResetParams& params)               
  : CBNumOutput(params), 
    m_resetParams(params)
{
  double dCurrentRatio = m_resetParams.GetReset().GetNominal()
    / m_resetParams.GetReset().GetCurrentConversionPrice();

  m_dInverseCurrentRatio = 1.0 / dCurrentRatio;

  // If we get this far, there must be at least one reset time
  m_dFirstResetTime = m_resetParams.GetFirstActiveResetTime();

}

void ResetNumOutput::UpdateMe(BackwardInstData& instdata, double dTime)
{
  // Only save the data if the current time is equal or before the
  // first reset time.
  if (IsEqualOrBefore(dTime, m_dFirstResetTime))
    BackwardNumOutput::UpdateMe(instdata, dTime);
}

void ResetNumOutput::UpdateMeAtEndOfGrid
                     (BackwardInstData& instdata, double dTime)
{
  // Only save the data if the current time is equal or before the
  // first reset time.
  if (IsEqualOrBefore(dTime, m_dFirstResetTime))
    BackwardNumOutput::UpdateMeAtEndOfGrid(instdata, dTime);
}

void ResetNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  DomainGeneral& domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for (size_t nIdx = 0; nIdx < instdata.m_nNbS; nIdx++)
    pdSpots[nIdx] = instdata.m_pdS[nIdx] * m_dInverseCurrentRatio;

  domain.AddSpotsAtTime(pdSpots, dTime);
  
  SaveSurfaceDataFrom(instdata);
}

void ResetNumOutput::SaveSurfaceAtEndOfGrid
     (BackwardInstData& instdata, double dTime)
{
  DomainGeneral& domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for (size_t nIdx = 0; nIdx < instdata.m_nNbS; nIdx++)
    pdSpots[nIdx] = instdata.m_pdS[nIdx] * m_dInverseCurrentRatio;

  domain.AddSpotsAtTime(pdSpots, dTime, true);
  
  SaveSurfaceDataFrom(instdata, true);
}

void ResetNumOutput::Finalize(BackwardInstData& instdata)
{
  // Save the final data, as usual
  CBNumOutput::Finalize(instdata);

  // Scale the final grid
  size_t nNbSpots = m_pdFinalSpots.size();
  for (size_t nIdx = 0; nIdx < nNbSpots; nIdx++)
    m_pdFinalSpots[nIdx] *= m_dInverseCurrentRatio;
}

void ResetNumOutput::CalculateFinalScalarResult(BackwardInstData& instdata)
{
  // create a temporary scaled grid
  std::vector<double> pdTmpGrid(m_nNbS);
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    pdTmpGrid[nIdx] = m_pdS[nIdx] * m_dInverseCurrentRatio;

  // Save pointer to old grid and update
  const double* pdOldGrid = m_pdS;
  m_pdS = &pdTmpGrid[0];

  // now save data as usual. This function uses m_pdS for the grid
  CBNumOutput::CalculateFinalScalarResult(instdata);

  // undo the grid change
  m_pdS = pdOldGrid;
}

void ResetNumOutput::SaveAnalysisData(BackwardInstData& instdata)
{
  // Save as usual
  BackwardNumOutput::SaveAnalysisData(instdata);

  // Scale the grid
  size_t nNbSpots = m_pdAnalysisSpots.size();
  for (size_t nIdx = 0; nIdx < nNbSpots; nIdx++)
    m_pdAnalysisSpots[nIdx] *= m_dInverseCurrentRatio;
}
