/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/heronumoutput.cpp
// Purpose:     implementation of NumOutput for Hero
// Created:     2005/09/26
// RCS-ID:      $Id: heronumoutput.cpp,v 1.7 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/hg/modeloutput.h"

#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"

#include "hg/backwardnumoutput.h"
#include "hg/heronumoutput.h"

// implement the AutoPtrDeleter for HeroNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::HeroNumOutput);
}

namespace ito33
{

  using namespace numeric;

namespace hg
{


void HeroNumOutput::Init(HeroInstData& instdata)
{  
  BackwardNumOutput::Init(instdata);
}


shared_ptr<ModelOutput> HeroNumOutput::GetModelOutput()
{  
  // Construct the output class
  shared_ptr<ModelOutput> pOutput( new ModelOutput() );

  // Only output the hero price, analysis date values, and surface. In
  // each case, output the square root of the values.
  HeroSQRT heroSQRT;
  pOutput->SetPrice( heroSQRT(m_dPrice) );    
    
  // Set the analysis date data, if available
  if ( m_dAnalysisTime > 0 )
  {
    size_t nNbS = m_pdAnalysisSpots.size();

    pOutput->SetSpotsAtAnalysisDate(m_pdAnalysisSpots);
    
    // Only save one regime in the finance output
    std::vector<double> pdAnalysisValues(nNbS);

    for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
      pdAnalysisValues[nIdx] = heroSQRT(m_pdAnalysisValues[nIdx]);

    pOutput->SetPricesAtAnalysisDate(pdAnalysisValues);
  }

  // If surfaces are needed, set them.
  if ( m_computationalFlags.GetComputeSurface() )
  {
    // Update the domain before creating the user output
    m_pDomain->GenerateOutputDates();

    pOutput->SetDomain(m_pDomain);

    // Need to convert from numeric price surface to finance price surface 
    // before setting. The former is not derived from the later
    m_ppPriceSurfaces[0]->ApplyFunctor( heroSQRT );
    shared_ptr<finance::SurfaceDouble> 
      pTmpSurface(new finance::SurfaceDouble(m_ppPriceSurfaces[0]));
    pOutput->SetPriceSurface(pTmpSurface);
  }

  return pOutput;

}


} // namespace hg

} // namespace ito33
