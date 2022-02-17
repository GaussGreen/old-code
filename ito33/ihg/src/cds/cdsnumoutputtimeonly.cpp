/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cds/cdsnumoutputtimeonly.cpp
// Purpose:     implementation of time only cds NumOutput class 
// Created:     2004/03/18
// RCS-ID:      $Id: cdsnumoutputtimeonly.cpp,v 1.42 2006/08/21 15:00:17 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfacezero.h"

#include "ito33/pricing/cdsparams.h"
#include "ito33/pricing/cdsmeshmanager.h"

#include "ihg/cdsnumoutputtimeonly.h"
#include "ihg/cdsinstdatatimeonly.h"

#include "ito33/ihg/modeloutput.h"

using ito33::numeric::AreTimesEqual;

// implement the AutoPtrDeleter for CDSNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::CDSNumOutputTimeOnly);
}

namespace ito33
{

namespace ihg
{


CDSNumOutputTimeOnly::CDSNumOutputTimeOnly(pricing::CDSParams& params) 
                                         : BackwardNumOutput(), 
                                           m_params(params)
{
}

void CDSNumOutputTimeOnly::Init(CDSInstDataTimeOnly& /* instdata */)
{
  if ( m_computationalFlags.GetComputeSurface() )
  {
    // Constructs a trivial domain
    m_pFixedDomain = make_ptr( new numeric::DomainFixedSpaceMesh
                                   ( m_params.GetSpotSharePrice() ) );

    m_pPriceSurface = make_ptr( new numeric::SurfaceGeneral(m_pFixedDomain) );
  }

  // set the analysis time
  m_dAnalysisTime = m_params.GetAnalysisTime();
}


void CDSNumOutputTimeOnly::UpdateMe(CDSInstDataTimeOnly& instdata, double dTime)
{
  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pFixedDomain->AddTime(dTime);

    finance::Values pdPricesTmp(1);
    pdPricesTmp[0] = instdata.m_dPrice;

    m_pPriceSurface->Add(pdPricesTmp);
  }

  if ( AreTimesEqual(m_dAnalysisTime, dTime) )
  {
    m_dPriceAtAnalysisDate = instdata.m_dPrice;
    m_dThetaAtanalysisDate = - instdata.m_dInverseTimeStep
                           * (instdata.m_dPrice - instdata.m_dOldPrice);

    m_pdAnalysisValues.clear();
    m_pdAnalysisValues.push_back(m_dPriceAtAnalysisDate);
  }

  // Always update from the current price 
  m_dPrice = instdata.m_dPrice;
}

void CDSNumOutputTimeOnly::Finalize(CDSInstDataTimeOnly& instdata)
{
  m_dTheta = - instdata.m_dInverseTimeStep
           * (instdata.m_dPrice - instdata.m_dOldPrice);

  m_dValueAfterDefault = instdata.m_dRecoveryValue;
}

// Return the requested data to the user
shared_ptr<ito33::ihg::ModelOutput> CDSNumOutputTimeOnly::GetOutput()
{
  // Construct the output class
  shared_ptr<ihg::ModelOutput> pOutput (new ihg::ModelOutput);

  pOutput->SetPrice(m_dPrice);

  pOutput->SetDelta(0);

  pOutput->SetGamma(0);

  pOutput->SetTheta(m_dTheta);
  
  pOutput->SetValueAfterDefault(m_dValueAfterDefault);

  if ( m_computationalFlags.GetComputeVega() )
    pOutput->SetVega(0.);

  if (m_dAnalysisTime > 0)
  {
    finance::Domain::Spots spots(1);
    spots[0] = m_params.GetSpotSharePrice();

    pOutput->SetSpotsAtAnalysisDate(spots);

    finance::SurfaceDouble::Doubles valuesTmp(1);
    
    pOutput->SetDeltasAtAnalysisDate(valuesTmp);
    pOutput->SetGammasAtAnalysisDate(valuesTmp);
    
    if ( m_computationalFlags.GetComputeVega() )
      pOutput->SetVegasAtAnalysisDate(valuesTmp);

    pOutput->SetPricesAtAnalysisDate(m_pdAnalysisValues);

    valuesTmp[0] = m_dThetaAtanalysisDate;
    pOutput->SetThetasAtAnalysisDate(valuesTmp);
  }

  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pFixedDomain->GenerateOutputDates();
    
    pOutput->SetDomain(m_pFixedDomain);

    pOutput->SetPriceSurface( shared_ptr<finance::SurfaceDouble>
                              ( new finance::SurfaceDouble(m_pPriceSurface) )
                            );

    shared_ptr<finance::SurfaceDouble>
      pSurfaceZero( new finance::SurfaceDouble(
                         shared_ptr<numeric::SurfaceZero>(
                        (new numeric::SurfaceZero(m_pFixedDomain)))) );

    pOutput->SetDeltaSurface(pSurfaceZero);
    pOutput->SetGammaSurface(pSurfaceZero);

    shared_ptr<numeric::SurfaceDouble> pThetaSurface;
    m_pPriceSurface->GetThetaBackwardOnly(pThetaSurface);

    pOutput->SetThetaSurface(shared_ptr<finance::SurfaceDouble>
                              ( new finance::SurfaceDouble(pThetaSurface) )
                            );

    if ( m_computationalFlags.GetComputeVega() )
      pOutput->SetVegaSurface( pSurfaceZero  );
  }
  
  return pOutput;
}

void CDSNumOutputTimeOnly::SaveSurface
     (BackwardInstData& /* instdata */, double /* dTime */)
{
}

} // namespace ihg

} // namespace ito33
