/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/numoutputtimeonly.cpp
// Purpose:     implementation of time only NumOutput class using HG model
// Created:     2005/06/09
// RCS-ID:      $Id: numoutputtimeonly.cpp,v 1.15 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/hg/modeloutput.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfacezero.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"

#include "hg/instdatatimeonly.h"
#include "hg/numoutputtimeonly.h"
#include "hg/numoutput.h"

// implement the AutoPtrDeleter for CDSNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::NumOutputTimeOnly);
}

namespace ito33
{

namespace hg
{

  using namespace numeric;

NumOutputTimeOnly::NumOutputTimeOnly(pricing::Params& params) 
                         : BackwardNumOutput(params)
{
}

void NumOutputTimeOnly::Init(InstDataTimeOnly& instdata)
{
  m_nNbRegimes = instdata.m_nNbRegimes;

  if ( m_computationalFlags.GetComputeSurface() )
  {
    // Constructs a trivial domain
    m_pFixedDomain = make_ptr( new DomainFixedSpaceMesh
                                   ( m_params.GetSpotSharePrice() ) );

    // price surfaces
    m_ppPriceSurfaces.resize(m_nNbRegimes);
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_ppPriceSurfaces[nIdxR] = make_ptr( new SurfaceGeneral(m_pFixedDomain) );

  }

  // set the analysis time
  m_dAnalysisTime = m_params.GetAnalysisTime();

  // Set sensitivity flags
  m_computationalFlags.SetSensitivityFlags( instdata.m_pbComputeSensitivities );

  // Clear the recovery value vector so we can safely push_back
  m_pdRecoveryValues.clear();

  // Clear the times so we can push_back
  m_pdTimes.clear();

}

void NumOutputTimeOnly::UpdateMe(InstDataTimeOnly& instdata, double dTime)
{
  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pFixedDomain->AddTime(dTime);

    finance::Values pdPricesTmp(1);
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {      
      pdPricesTmp[0] = instdata.m_pdPrices[nIdxR];

      m_ppPriceSurfaces[nIdxR]->Add(pdPricesTmp);
    }
  }

  if ( AreTimesEqual(m_dAnalysisTime, dTime) )
  {
    m_pdAnalysisSpots.clear();
    m_pdAnalysisSpots.push_back( m_params.GetSpotSharePrice() );

    m_dPriceAtAnalysisDate = instdata.m_pdPrices[0];
    m_dThetaAtanalysisDate = - instdata.m_dInverseTimeStep
                           * (instdata.m_pdPrices[0] - instdata.m_pdOldPrices[0]);

    m_pdAnalysisValues.clear();
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_pdAnalysisValues.push_back(instdata.m_pdPrices[nIdxR]);   
  }

  // Save the recovery values, and the time mesh
  size_t nNbTimes = m_pdTimes.size();
  if ( nNbTimes > 0 && AreTimesEqual(dTime, m_pdTimes[nNbTimes-1]) )
  {
    // an event must have occured.  Save the new recovery value    
    m_pdRecoveryValues[nNbTimes-1] = instdata.m_dRecoveryValue;
  }
  else
  {
    m_pdRecoveryValues.push_back( instdata.m_dRecoveryValue );
    m_pdTimes.push_back( dTime );
  }

}

void NumOutputTimeOnly::Finalize(InstDataTimeOnly& instdata)
{
  m_dTheta = - instdata.m_dInverseTimeStep
           * (instdata.m_pdPrices[0] - instdata.m_pdOldPrices[0]);

  if ( m_computationalFlags.HasSensitivityFlags() )
  {
    m_pdSensitivities = instdata.GetSensitivities();

    // Sensitivity w.r.t post default volatility will just be zero
    if ( m_computationalFlags.GetSensitivityFlags().back() )
      m_pdSensitivities.push_back(0);
  }

  m_pdFinalValues.resize(m_nNbRegimes);
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    m_pdFinalValues[nIdxR] = instdata.m_pdPrices[nIdxR];

  m_dValueAfterDefault = instdata.m_dRecoveryValue;
}

void NumOutputTimeOnly::SetPrices(const std::vector<double>& pdPrices)
{
  m_pdFinalValues = pdPrices;
  m_pdFinalSpots.resize(1, m_params.GetSpotSharePrice());
}

// Return the requested data to the user
shared_ptr<ModelOutput> NumOutputTimeOnly::GetModelOutput()
{
  // Construct the output class
  shared_ptr<ModelOutput> pOutput (new ModelOutput);

  pOutput->SetPrice(m_pdFinalValues[0]);
  pOutput->SetDelta(0);
  pOutput->SetGamma(0);
  pOutput->SetTheta(m_dTheta);

  if ( m_computationalFlags.HasSensitivityFlags() )
  {
    shared_ptr<hg::NumOutput> pNumOutput( new hg::NumOutput() );
    pNumOutput->SetSensitivities(m_pdSensitivities);
    pOutput->SetNumOutput( pNumOutput );
  }

  pOutput->SetValueAfterDefault(m_dValueAfterDefault);
  
  // A temporary work around for analysis data at valuation date since
  // sometimes we computed analytically and only values at valuation date 
  // are set
  if ( AreTimesEqual(m_params.GetValuationTime(), m_params.GetAnalysisTime()) )
  {
    m_pdAnalysisValues = m_pdFinalValues;
    m_pdAnalysisSpots.resize(1, m_params.GetSpotSharePrice());
  }

  if (m_dAnalysisTime > 0)
  {
    finance::Domain::Spots spots(1);
    spots[0] = m_params.GetSpotSharePrice();

    pOutput->SetSpotsAtAnalysisDate(spots);

    finance::SurfaceDouble::Doubles valuesTmp(1);
    
    pOutput->SetDeltasAtAnalysisDate(valuesTmp);
    pOutput->SetGammasAtAnalysisDate(valuesTmp);
    
    pOutput->SetPricesAtAnalysisDate(m_pdAnalysisValues);

    valuesTmp[0] = m_dThetaAtanalysisDate;
    pOutput->SetThetasAtAnalysisDate(valuesTmp);
  }

  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pFixedDomain->GenerateOutputDates();
    
    pOutput->SetDomain(m_pFixedDomain);

    pOutput->SetPriceSurface
                 ( shared_ptr<finance::SurfaceDouble>
                   ( new finance::SurfaceDouble(m_ppPriceSurfaces[0]) ) );

    shared_ptr<finance::SurfaceDouble>
      pSurfaceZero( new finance::SurfaceDouble
          (shared_ptr<SurfaceZero>(new SurfaceZero(m_pFixedDomain))) );

    pOutput->SetDeltaSurface(pSurfaceZero);
    pOutput->SetGammaSurface(pSurfaceZero);

    shared_ptr<SurfaceDouble> pThetaSurface;
    m_ppPriceSurfaces[0]->GetThetaBackwardOnly(pThetaSurface);

    pOutput->SetThetaSurface(shared_ptr<finance::SurfaceDouble>
                              ( new finance::SurfaceDouble(pThetaSurface) )
                            );
  }

  return pOutput;
}

void NumOutputTimeOnly::SaveSurface
     (BackwardInstData& /* instdata */, double /* dTime */)
{
}


} // namespace hg

} // namespace ito33
