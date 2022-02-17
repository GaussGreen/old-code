/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/cboptionnumoutput.cpp
// Purpose:     implementation of CBOptionNumOutput class 
// Author:      Nabil
// Created:     2005/10/18
// RCS-ID:      $Id: cboptionnumoutput.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/domain_general.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/numeric/mesh/spacemesh.h"

#include "ito33/ihg/bondlikeoutput.h"
#include "ito33/ihg/cboptionoutput.h"

#include "ihg/cboptioninstdata.h"
#include "ihg/cboptionnumoutput.h"

using ito33::AutoPtr;
using ito33::shared_ptr;

using namespace ito33::pricing;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;

using namespace ito33::ihg;


// implement the AutoPtrDeleter for CBOptionNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::CBOptionNumOutput);
}


CBOptionNumOutput::CBOptionNumOutput(ito33::pricing::CBLikeParams& params) 
                  : BackwardNumOutput(), m_params(params)
{
  m_pCBNumOutput = shared_ptr<CBNumOutput>( new CBNumOutput(params) );
  m_bHasConstraintFlags = true;
}

void CBOptionNumOutput::Init(CBOptionInstData& cboptioninstdata)
{
  // Init the CBNumOutput
  m_pCBNumOutput->Init(*cboptioninstdata.GetCBInstData());
  
  // Vega if computed is by PDE
  m_computationalFlags.SetComputeVega(cboptioninstdata.m_bComputeVega);

  /// Flag indicates if we compute fugit
  m_computationalFlags.SetComputeFugit(cboptioninstdata.m_bComputeFugit);
  
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
    if( m_computationalFlags.GetComputeVega() )
      m_pVegaSurface = make_ptr( new SurfaceGeneral(m_pDomain) );

    // fugit surface (fugit is computed by PDE)
    if ( m_computationalFlags.GetComputeFugit() )
      m_pFugitSurface = make_ptr( new SurfaceGeneral(m_pDomain) );
  }
}

void CBOptionNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainGeneral*>(m_pDomain.get()), 
    "The type of Domain in CBOptionNumOutput should be DomainGeneral."
  );

  DomainGeneral& 
    domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for(size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime);

  SaveSurfaceDataFrom(instdata);
}

void CBOptionNumOutput::SaveSurfaceAtEndOfGrid
     (BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainGeneral*>(m_pDomain.get()), 
    "The type of Domain in CBOptionNumOutput should be DomainGeneral."
  );

  DomainGeneral& 
    domain = static_cast<DomainGeneral&>(*m_pDomain);
  finance::Domain::Spots pdSpots(instdata.m_nNbS);
  for(size_t n = 0; n < instdata.m_nNbS; n++)
    pdSpots[n] = instdata.m_pdS[n];
  
  domain.AddSpotsAtTime(pdSpots, dTime, true);
  
  SaveSurfaceDataFrom(instdata, true);
}

void CBOptionNumOutput::UpdateMe(BackwardInstData& instdata, double dTime)
{
  // Update for the CBNumOutput
  CBOptionInstData& 
    cboptioninstdata = dynamic_cast<CBOptionInstData&>(instdata);
  m_pCBNumOutput->UpdateMe(*cboptioninstdata.GetCBInstData(), dTime);

  BackwardNumOutput::UpdateMe(instdata, dTime);
}

void CBOptionNumOutput::UpdateMeAtEndOfGrid(BackwardInstData& instdata, 
                                            double dTime)
{
  // Update for the CBNumOutput
  CBOptionInstData& 
    cboptioninstdata = dynamic_cast<CBOptionInstData&>(instdata);
  m_pCBNumOutput->UpdateMeAtEndOfGrid(*cboptioninstdata.GetCBInstData(), 
                                      dTime);

  BackwardNumOutput::UpdateMeAtEndOfGrid(instdata, dTime);
}

void CBOptionNumOutput::Finalize(BackwardInstData& instdata)
{
  // Finalize for the CBNumOutput
  CBOptionInstData& 
    cboptioninstdata = dynamic_cast<CBOptionInstData&>(instdata);
  m_pCBNumOutput->Finalize(*cboptioninstdata.GetCBInstData());

  BackwardNumOutput::Finalize(instdata);
}

// Return the requested data to the user
shared_ptr<ito33::ihg::CBOptionOutput> CBOptionNumOutput::GetOutput()
{
  // Construct the output class
  shared_ptr<ihg::CBOptionOutput> pOutput (new ihg::CBOptionOutput());

  // CBOutput treatment
  pOutput->SetCBOutput( m_pCBNumOutput->GetOutput() );

  BackwardNumOutput::GetOutput(pOutput.get());

  return pOutput;
}
