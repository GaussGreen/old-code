/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/numoutputanalytical.cpp
// Purpose:     implementation of NumOutput class using analytical formula
// Created:     2006/08/02
// RCS-ID:      $Id: numoutputanalytical.cpp,v 1.4 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/interpolationmatrix.h"
#include "ito33/numeric/deltagamma.h"

#include "ito33/pricing/params.h"

#include "ito33/hg/modeloutput.h"
#include "hg/numoutputanalytical.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::NumOutputAnalytical);

namespace hg
{

  using namespace numeric;

NumOutputAnalytical::NumOutputAnalytical(pricing::Params& params) 
                                       : BackwardNumOutput(params)
{
}

void NumOutputAnalytical::SetFinalValues(const std::vector<double>& pdS,
                                         const std::vector<double>& pdPrices)
{
  m_pdFinalSpots = pdS;
  m_pdFinalValues = pdPrices;

  // Just copy the data to avoid complicated code
  if ( AreTimesEqual(m_params.GetAnalysisTime(), m_params.GetValuationTime()) )
  {
    m_pdAnalysisSpots = pdS;
    m_pdAnalysisValues = pdPrices;
  }
}

shared_ptr<ModelOutput> NumOutputAnalytical::GetModelOutput()
{
  // Construct the output class
  shared_ptr<ModelOutput> pOutput (new ModelOutput);
  
  // Scalar value
  pOutput->SetValueAfterDefault(m_dValueAfterDefault);

  pOutput->SetPrice(m_dPrice);
  pOutput->SetDelta(m_dDelta);
  pOutput->SetGamma(m_dGamma);

  // Analysis data is only supported at valuation date
  if ( !m_pdFinalSpots.empty() )
  {
    pOutput->SetSpotsAtAnalysisDate(m_pdFinalSpots);
    pOutput->SetPricesAtAnalysisDate(m_pdFinalValues);

    m_nNbS = m_pdFinalSpots.size();
    m_pdS = &m_pdFinalSpots[0];

    // Temporary vector for delta/gamma etc
    std::vector<double> pdResults(m_nNbS);
    
    // Delta
    ComputeDelta(m_pdS, &m_pdFinalValues[0], m_nNbS, &pdResults[0]);

    pOutput->SetDeltasAtAnalysisDate(pdResults);

    ComputeGammaFD(m_pdS, &m_pdFinalValues[0], m_nNbS, &pdResults[0]);

    pOutput->SetGammasAtAnalysisDate(pdResults);
  }

  // No surface data is available

  return pOutput;
}

} // namespace hg

} // namespace ito33
