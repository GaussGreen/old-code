/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/cbnumoutput.cpp
// Purpose:     implementation of CBNumOutput class 
// Author:      Nabil
// Created:     2004/04/13
// RCS-ID:      $Id: cbnumoutput.cpp,v 1.32 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/domain_general.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/numeric/mesh/spacemesh.h"

#include "ito33/ihg/bondlikeoutput.h"

#include "ihg/cbinstdata.h"
#include "ihg/cbnumoutput.h"
#include "ihg/backwardinstdata.h"

using ito33::AutoPtr;
using ito33::shared_ptr;

using namespace ito33::numeric;
using namespace ito33::numeric::mesh;

using namespace ito33::ihg;


// implement the AutoPtrDeleter for CBNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::CBNumOutput);
}

void CBNumOutput::Init(CBInstData& cbinstdata)
{
  /// Flag indicates if we compute fugit
  m_computationalFlags.SetComputeFugit(cbinstdata.m_bComputeFugit);
  
  // Vega if computed is by PDE
  m_computationalFlags.SetComputeVega(cbinstdata.m_bComputeVega);
  
  // set the analysis time
  m_dAnalysisTime = m_params.GetAnalysisTime();
  
  m_bHasConstraintFlags = true;

  // If surfaces are requested, construct them. The first step is constructing
  // the underlying domain which is shared between the surfaces.
  if( m_computationalFlags.GetComputeSurface() )
  {
    m_pDomain = shared_ptr<Domain>( new DomainGeneral() );
    
    // price surface
    m_pPriceSurface = make_ptr( new SurfaceGeneral(m_pDomain) );

    m_pFlagSurface = make_ptr( new SurfaceFlag(m_pDomain) );

    // vega surface (vega is computed by PDE)
    if ( m_computationalFlags.GetComputeVega() )
      m_pVegaSurface = make_ptr( new SurfaceGeneral(m_pDomain) );

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
  for (size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime);

  SaveSurfaceDataFrom(instdata);
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
  for (size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime, true);
  
  SaveSurfaceDataFrom(instdata, true);
}

void CBNumOutput::CalculateFinalScalarResult(BackwardInstData& instdata)
{
  BackwardNumOutput::CalculateFinalScalarResult(instdata);
  
  double dSpot = 0.;
 
  // Interpolate for bond floor at the spot
  Interpolate(m_pdS, instdata.m_pdPrices.Get(),
              m_nNbS, &dSpot, &m_dBondFloor, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  if (m_dBondFloor < 0.)
    m_dBondFloor = 0.;
}

// Return the requested data to the user
shared_ptr<BondLikeOutput> CBNumOutput::GetOutput()
{
  // Construct the output class
  shared_ptr<BondLikeOutput> pOutput (new BondLikeOutput());

  BackwardNumOutput::GetOutput(pOutput.get());

  pOutput->SetBondFloor(m_dBondFloor);

  return pOutput;
}
