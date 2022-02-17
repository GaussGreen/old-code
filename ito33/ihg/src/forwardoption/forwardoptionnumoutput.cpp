/////////////////////////////////////////////////////////////////////////////
// Name:        forwardoption/forwardoptionnumoutput.cpp
// Purpose:     implementation of Option NumOutput class using forward PDE 
// Author:      Wang
// RCS-ID:      $Id: forwardoptionnumoutput.cpp,v 1.14 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ihg/forwardoptioninstdata.h"
#include "ihg/forwardoptionnumoutput.h"

using ito33::shared_ptr;

using namespace ito33::numeric;
using namespace ito33;

using ito33::ihg::ForwardOptionNumOutput;
using ito33::ihg::ForwardOptionInstData;

// implement the AutoPtrDeleter for OptionTotalData
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::ForwardOptionNumOutput);
}

void ForwardOptionNumOutput::Init(ForwardOptionInstData& instdata)
{
  // The number of space points
  m_nNbS = instdata.m_nNbS;
  const double *pdS = instdata.m_pdS;

  // We are using fix mesh so we can save the space mesh here.
  // We always save m_pdS as it might be needed for interpolation, output etc. 
  m_pdS.resize(m_nNbS);
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdS[nIdx] = pdS[nIdx];

  // If surfaces are requested, construct them. The first step is constructing
  // the underlying domain which is shared between the surfaces.
  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pDomain = make_ptr( new DomainFixedSpaceMesh(pdS, m_nNbS) );
    
    // price surface
    m_pPriceSurface = make_ptr( new SurfaceGeneral(m_pDomain) );
  }
}

void ForwardOptionNumOutput::UpdateMe
     (ForwardOptionInstData &instdata, double dTime)
{
  size_t nIdx;

  // Always save the price data pointer
  m_pdTmpPrices = instdata.m_pdPrices.Get();

  // If surfaces are being computed, then save the data at this step
  // This function is called by init, so initial data is automatically stored.
  // Also remember to update the shared underlying domain!
  if ( m_computationalFlags.GetComputeSurface() )
  {
    m_pDomain->AddTime(dTime);

    // In all cases, convert the data to the right format before adding
    // to the surface.  Do the price first    
    numeric::SurfaceDouble::Doubles pdTmpData;
    pdTmpData.resize(m_nNbS);
    for (nIdx = 0; nIdx < m_nNbS; nIdx++)
      pdTmpData[nIdx] = m_pdTmpPrices[nIdx];

    m_pPriceSurface->Add(pdTmpData);

  } // if computing surfaces
}

void ForwardOptionNumOutput::Finalize(ForwardOptionInstData& /*instdata*/)
{
  size_t nIdx;

  // Copy the final values stored in instdata, since the life of instdata may be
  // less than the life of numoutput
  m_pdValues = CountedArray<double>(m_nNbS);
  for (nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdValues[nIdx] = m_pdTmpPrices[nIdx];
}

shared_ptr<finance::ModelOutput> ForwardOptionNumOutput::GetModelOutput()
{ 
  // Construct the output class
  shared_ptr<finance::ModelOutput> pOutput (new finance::ModelOutput());

  // Scalar data is found at the spot
  double dSpot = m_params.GetSpotSharePrice();

  // Save the price data
  double dPrice;
  Interpolate(&m_pdS[0], m_pdValues.Get(),
              m_nNbS, &dSpot, &dPrice, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  pOutput->SetPrice(dPrice);

  return pOutput;
}

