/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/src/analysisdate_surface.cpp
// Purpose:     Tests on Analysis date and ComputeSurface flags
// Created:     2004/11/22
// RCS-ID:      $Id: analysisdate_surface.cpp,v 1.7 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/surfacedouble.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/tests/analysisdate_surface.h"

namespace ito33
{

namespace ihg
{


bool AnalysisDateTest::NoDataForAnalysisDateBeforeValuationDate()
{
  Date analysisDate = m_pDerivative->GetSessionData()->GetValuationDate();
  analysisDate.AddMonths(-1);

  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);
  pFlags->SetAnalysisDate(analysisDate);

  m_model.SetExternalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMo = m_model.Compute(*m_pDerivative);
  
  return !pMo->HasSpotAtAnalysisDate();
}

bool AnalysisDateTest::NoDataForAnalysisDateAtMaturityDate()
{
  Date analysisDate = m_pDerivative->GetMaturityDate();

  shared_ptr<finance::ComputationalFlags> pFlags(new finance::ComputationalFlags);
  pFlags->SetAnalysisDate(analysisDate);

  m_model.SetExternalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMo = m_model.Compute(*m_pDerivative);
  
  return !pMo->HasSpotAtAnalysisDate();
}

void AnalysisDateTest::DataAvailability()
{
  // Check analysis date = valuation date
  Date analysisDate = m_pDerivative->GetSessionData()->GetValuationDate();

  shared_ptr<finance::ComputationalFlags> pFlags(new finance::ComputationalFlags);
  pFlags->SetComputeFugit(true);
  pFlags->SetComputeRho(true);
  pFlags->SetComputeVega(true);
  pFlags->SetComputeFugit(true);

  pFlags->SetAnalysisDate(analysisDate);

  m_model.SetExternalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMo = m_model.Compute(*m_pDerivative);

  if ( !pMo->HasSpotAtAnalysisDate() )
    return;

  const size_t nNbSpots = pMo->GetSpotsAtAnalysisDate().size();

  // Price
  m_pbHasDatas[0] =  pMo->HasPriceAtAnalysisDate()
                  && pMo->GetPricesAtAnalysisDate().size() == nNbSpots;

  // Delta
  m_pbHasDatas[1] =  pMo->HasDeltaAtAnalysisDate()
                  && pMo->GetDeltasAtAnalysisDate().size() == nNbSpots;

  // Gamma
  m_pbHasDatas[2] =  pMo->HasGammaAtAnalysisDate()
                  && pMo->GetGammasAtAnalysisDate().size() == nNbSpots;

  // Theta
  m_pbHasDatas[3] =  pMo->HasThetaAtAnalysisDate()
                  && pMo->GetThetasAtAnalysisDate().size() == nNbSpots;

  // Rho
  m_pbHasDatas[4] =  pMo->HasRhoAtAnalysisDate()
                  && pMo->GetRhosAtAnalysisDate().size() == nNbSpots;
 
  // Vega
  m_pbHasDatas[5] =  pMo->HasVegaAtAnalysisDate()
                  && pMo->GetVegasAtAnalysisDate().size() == nNbSpots;

  // fugit
  m_pbHasDatas[6] =  pMo->HasFugitAtAnalysisDate()
                  && pMo->GetFugitsAtAnalysisDate().size() == nNbSpots;
}

bool AnalysisDateTest::ConsistencyAtValuationDate()
{
  Date analysisDate = m_pDerivative->GetSessionData()->GetValuationDate();
  shared_ptr<finance::ComputationalFlags> pFlags(new finance::ComputationalFlags);

  pFlags->SetAnalysisDate(analysisDate);

  m_model.SetExternalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMo = m_model.Compute(*m_pDerivative);

  double dSpot = m_pDerivative->GetSessionData()->GetSpotSharePrice();
  double dPrice;

  if (!pMo->HasSpotAtAnalysisDate())
    return false;

  // Interpolation on the analysis data
  numeric::Interpolate(&pMo->GetSpotsAtAnalysisDate()[0],
                       &pMo->GetPricesAtAnalysisDate()[0],
                       pMo->GetSpotsAtAnalysisDate().size(),
                       &dSpot,
                       &dPrice,
                       1,
                       numeric::ExtrapolationMode_Linear,
                       numeric::ExtrapolationMode_Linear);

  return fabs(dPrice - pMo->GetPrice()) < 1.e-3;
}

void AnalysisDateTest::Report()
{
  if (!NoDataForAnalysisDateBeforeValuationDate())
    std::cout << "Having datas at an analysis date before valuation date\n";

  if (!NoDataForAnalysisDateAtMaturityDate())
    std::cout << "Having analysis data at maturity date!\n";

  if (!ConsistencyAtValuationDate())
    std::cout << "Analysis data at valuation date not consistent with "
                 "scalar values\n";

  DataAvailability();

  for (size_t nIdx = 0; nIdx < NUMBEROFGREEK; nIdx++)
    if (m_pbShouldHaveDatas[nIdx] != m_pbHasDatas[nIdx])
      std::cout << CanCompute[nIdx] << " expected  but not computed\n";
}

void SurfaceTest::DataAvailability()
{
  shared_ptr<finance::ComputationalFlags> pFlags(new finance::ComputationalFlags);
  pFlags->SetComputeSurface(true);

  pFlags->SetComputeFugit(true);
  pFlags->SetComputeRho(true);
  pFlags->SetComputeVega(true);
  pFlags->SetComputeFugit(true);

  m_model.SetExternalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMo = m_model.Compute(*m_pDerivative);

  if ( !pMo->HasDomain() )
    return;

  // Price
  m_pbHasDatas[0] = pMo->HasPriceSurface();

  // Delta
  m_pbHasDatas[1] = pMo->HasDeltaSurface();

  // Gamma
  m_pbHasDatas[2] = pMo->HasGammaSurface();

  // Theta
  m_pbHasDatas[3] = pMo->HasThetaSurface();

  // Rho
  m_pbHasDatas[4] = pMo->HasRhoSurface();
 
  // Vega
  m_pbHasDatas[5] = pMo->HasVegaSurface();

  // Fugit
  m_pbHasDatas[6] = pMo->HasFugitSurface();
}

bool SurfaceTest::ConsistencyAtValuationDate()
{
  Date analysisDate = m_pDerivative->GetSessionData()->GetValuationDate();
  shared_ptr<finance::ComputationalFlags> pFlags(new finance::ComputationalFlags);

  pFlags->SetAnalysisDate(analysisDate);
  
  pFlags->SetComputeSurface(true);

  m_model.SetExternalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMo = m_model.Compute(*m_pDerivative);
  const std::vector<double>& 
    pdPricesAtAnalysisDate = pMo->GetPricesAtAnalysisDate();
  const std::vector<double>& 
    pdSpotsAtAnalysisDate = pMo->GetSpotsAtAnalysisDate();

  shared_ptr<finance::SurfaceDouble> pSurface = pMo->GetPriceSurface();
  shared_ptr<finance::Domain> pDomain = pMo->GetDomain();

  pDomain->SetUnderlyingSharePrices(pdSpotsAtAnalysisDate);

  const std::vector<double>&
    pPrices = pSurface->GetValuesAt(pDomain->GetDates().size() - 1);

  for (size_t nIdx = 0; nIdx < pdSpotsAtAnalysisDate.size(); nIdx++)
    if (fabs(pPrices[nIdx] - pdPricesAtAnalysisDate[nIdx]) > 1.e-4)
      return false;

  return true;
}

void SurfaceTest::Report()
{
  if (!ConsistencyAtValuationDate())
    std::cout << "surface data not consistent with data at valuation date.\n"; 

  DataAvailability();

  for (size_t nIdx = 0; nIdx < NUMBEROFGREEK; nIdx++)
    if (m_pbShouldHaveDatas[nIdx] != m_pbHasDatas[nIdx])
      std::cout << CanCompute[nIdx] << " expected  but not computed\n";
}


} // namespace ihg

} // namespace ito33
