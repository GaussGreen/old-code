/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cbnumoutput.cpp
// Purpose:     implementation of CBNumOutput class 
// Created:     2005/04/11
// RCS-ID:      $Id: cbnumoutput.cpp,v 1.10 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/domain_general.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/hg/bondlikeoutput.h"

#include "hg/cbinstdata.h"
#include "hg/cbnumoutput.h"
#include "hg/backwardinstdata.h"

// implement the AutoPtrDeleter for CBNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::CBNumOutput);
}

namespace ito33
{

namespace hg
{

  using namespace numeric;

CBNumOutput::CBNumOutput(pricing::CBLikeParams& params) 
                       : BackwardNumOutput(params), m_params(params)
{
  m_bHasConstraintFlags = true;
}

void CBNumOutput::Init(CBInstData& cbinstdata)
{
  m_nNbRegimes = cbinstdata.m_nNbRegimes;

  /// Flag indicates if we compute fugit
  m_computationalFlags.SetComputeFugit(cbinstdata.m_bComputeFugit);
  
  // set the analysis time
  m_dAnalysisTime = m_params.GetAnalysisTime();
  
  m_bHasConstraintFlags = true;

  // If surfaces are requested, construct them. The first step is constructing
  // the underlying domain which is shared between the surfaces.
  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pDomain = make_ptr( new DomainGeneral() );
    
    // price surfaces
    m_ppPriceSurfaces.resize(m_nNbRegimes);
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_ppPriceSurfaces[nIdxR] = make_ptr( new SurfaceGeneral(m_pDomain) );

    m_pFlagSurface = make_ptr( new SurfaceFlag(m_pDomain) );

    // fugit surface (fugit is computed by PDE)
    if ( m_computationalFlags.GetComputeFugit() )
      m_pFugitSurface = make_ptr( new SurfaceGeneral(m_pDomain) );
  }

}

void CBNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainGeneral*>(m_pDomain.get()), 
    "The type of Domain in CBNumOutput should be DomainGeneral."
  );

  DomainGeneral& 
    domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for(size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime);

  BackwardNumOutput::SaveSurfaceDataFrom(instdata);
}

void CBNumOutput::SaveSurfaceAtEndOfGrid
     (BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainGeneral*>(m_pDomain.get()), 
    "The type of Domain in CBNumOutput should be DomainGeneral."
  );

  DomainGeneral& 
    domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for(size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime, true);
  
  BackwardNumOutput::SaveSurfaceDataFrom(instdata, true);
}

void CBNumOutput::CalculateFinalScalarResult(BackwardInstData& instdata)
{
  BackwardNumOutput::CalculateFinalScalarResult(instdata);
 
  double dSpot = 0.;
 
  // Interpolate for price at the spot
  Interpolate(m_pdS, instdata.m_pdPrices.Get(),
              m_nNbS, &dSpot, &m_dBondFloor, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  if ( m_dBondFloor < 0. )
    m_dBondFloor = 0.;
}

// Return the requested data to the user
shared_ptr<BondLikeOutput> CBNumOutput::GetBondLikeOutput()
{
  // Construct the output class
  shared_ptr<BondLikeOutput> pOutput (new BondLikeOutput());

  GetOutput(*pOutput);

  pOutput->SetBondFloor(m_dBondFloor);

  return pOutput;
}


} // namespace hg

} // namespace ito33
