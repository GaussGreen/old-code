/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cboptionnumoutput.cpp
// Purpose:     implementation of CBOptionNumOutput class 
// Created:     2006/01/19
// RCS-ID:      $Id: cboptionnumoutput.cpp,v 1.3 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/domain_general.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/numeric/mesh/spacemesh.h"

#include "hg/cbnumoutput.h"
#include "hg/cboptioninstdata.h"
#include "hg/cboptionnumoutput.h"

#include "ito33/hg/bondlikeoutput.h"
#include "ito33/hg/cboptionoutput.h"

// implement the AutoPtrDeleter for CBOptionNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::CBOptionNumOutput);
}

namespace ito33
{

namespace hg
{

  using namespace numeric;

CBOptionNumOutput::CBOptionNumOutput(pricing::CBLikeParams& params) 
                  : BackwardNumOutput(params), m_params(params),
                    m_pCBNumOutput( new CBNumOutput(params) )
{
  m_bHasConstraintFlags = true;
}

void CBOptionNumOutput::Init(CBOptionInstData& cboptioninstdata)
{
  m_nNbRegimes = cboptioninstdata.m_nNbRegimes;

  // Init the CBNumOutput
  m_pCBNumOutput->Init( cboptioninstdata.GetCBInstData() );
  
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

void CBOptionNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainGeneral*>( m_pDomain.get() ), 
    "The type of Domain in CBOptionNumOutput should be DomainGeneral."
  );

  DomainGeneral& domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for (size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime);

  SaveSurfaceDataFrom(instdata);
}

void CBOptionNumOutput::SaveSurfaceAtEndOfGrid
     (BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainGeneral*>( m_pDomain.get() ), 
    "The type of Domain in CBOptionNumOutput should be DomainGeneral."
  );

  DomainGeneral& domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for (size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime, true);
  
  SaveSurfaceDataFrom(instdata, true);
}

void CBOptionNumOutput::UpdateMe(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<CBOptionInstData*>(&instdata), 
    "The type of instdata in CBOptionNumOutput should be CBOptionInstData."
  );

  // Update for the CBNumOutput
  CBOptionInstData& 
    cboptioninstdata = static_cast<CBOptionInstData&>(instdata);
  m_pCBNumOutput->UpdateMe(cboptioninstdata.GetCBInstData(), dTime);

  BackwardNumOutput::UpdateMe(instdata, dTime);
}

void CBOptionNumOutput::UpdateMeAtEndOfGrid(BackwardInstData& instdata, 
                                            double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<CBOptionInstData*>(&instdata), 
    "The type of instdata in CBOptionNumOutput should be CBOptionInstData."
  );

  // Update for the CBNumOutput
  CBOptionInstData& 
    cboptioninstdata = static_cast<CBOptionInstData&>(instdata);
  m_pCBNumOutput->UpdateMeAtEndOfGrid(cboptioninstdata.GetCBInstData(), dTime);

  BackwardNumOutput::UpdateMeAtEndOfGrid(instdata, dTime);
}

void CBOptionNumOutput::Finalize(BackwardInstData& instdata)
{
  ASSERT_MSG
  (
    dynamic_cast<CBOptionInstData*>(&instdata), 
    "The type of instdata in CBOptionNumOutput should be CBOptionInstData."
  );

  // Finalize for the CBNumOutput
  CBOptionInstData& 
    cboptioninstdata = static_cast<CBOptionInstData&>(instdata);
  m_pCBNumOutput->Finalize( cboptioninstdata.GetCBInstData() );

  BackwardNumOutput::Finalize(instdata);
}

// Return the requested data to the user
shared_ptr<CBOptionOutput> CBOptionNumOutput::GetCBOptionOutput()
{
  // Construct the output class
  shared_ptr<CBOptionOutput> pOutput( new CBOptionOutput() );

  // CBOutput treatment
  pOutput->SetCBOutput( m_pCBNumOutput->GetBondLikeOutput() );

  GetOutput(*pOutput);

  return pOutput;
}

} // namespace hg

} // namespace ito33
