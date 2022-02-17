/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/parbond/parbondnumoutputtimeonly.cpp
// Purpose:     implementation of time only parbond NumOutput class 
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondnumoutputtimeonly.cpp,v 1.6 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfacezero.h"

#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

#include "ihg/parbondnumoutputtimeonly.h"
#include "ihg/parbondinstdatatimeonly.h"

#include "ito33/ihg/modeloutput.h"

using ito33::numeric::AreTimesEqual;

// implement the AutoPtrDeleter for ParBondNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::ParBondNumOutputTimeOnly);
}

namespace ito33
{

namespace ihg
{


ParBondNumOutputTimeOnly::ParBondNumOutputTimeOnly(pricing::ParBondParams& params) 
                                         : BackwardNumOutput(), 
                                           m_params(params)
{
}

void ParBondNumOutputTimeOnly::Init(ParBondInstDataTimeOnly& /* instdata */)
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


void ParBondNumOutputTimeOnly::UpdateMe(ParBondInstDataTimeOnly& instdata, double dTime)
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

void ParBondNumOutputTimeOnly::Finalize(ParBondInstDataTimeOnly& instdata)
{
  m_dTheta = - instdata.m_dInverseTimeStep
           * (instdata.m_dPrice - instdata.m_dOldPrice);

  m_dValueAfterDefault = instdata.m_dRecoveryValue;
}

// Return the requested data to the user
shared_ptr<ModelOutput> ParBondNumOutputTimeOnly::GetOutput()
{
  // Construct the output class
  shared_ptr<ModelOutput> pOutput (new ModelOutput);

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
      pSurfaceZero( new finance::SurfaceDouble
                        (shared_ptr<numeric::SurfaceZero>
                          (new numeric::SurfaceZero(m_pFixedDomain)) ) );

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

void ParBondNumOutputTimeOnly::SaveSurface
     (BackwardInstData& /* instdata */, double /* dTime */)
{
}

} // namespace ihg

} // namespace ito33
