/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/option/numoutput_backward.cpp
// Purpose:     implementation of BackwardNumOutput class 
// Author:      Based on ICARE version
// RCS-ID:      $Id: backwardnumoutput.cpp,v 1.27 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/pricing/backwardmeshmanager_fix.h"
#include "ito33/pricing/params.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/computationalflags.h"

#include "ihg/backwardnumoutput.h"
#include "ihg/backwardinstdata.h"

using ito33::shared_ptr;
using ito33::AutoPtr;

using namespace ito33::numeric;

using ito33::ihg::BackwardNumOutput;
using ito33::ihg::BackwardInstData;

void BackwardNumOutput::InitWithFixMesh
                        (
                          BackwardInstData& instdata,
                          const ito33::pricing::Params& params,
                          bool bHasConstraintFlags
                        )
{
  // Get the space mesh
  m_pdS = instdata.GetSpaceMesh(m_nNbS);

  // Vega if computed is by PDE
  m_computationalFlags.SetComputeVega(instdata.m_bComputeVega);

  // Check if instdata has constraint flags
  m_bHasConstraintFlags = bHasConstraintFlags;

  // If surfaces are requested, construct them. The first step is constructing
  // the underlying domain which is shared between the surfaces.
  if ( m_computationalFlags.GetComputeSurface() )
  {
    // copy the mesh in case extra boundary points need to be added
    std::vector<double> pdS = CopySpaceMesh(instdata);

    m_pDomain = make_ptr( new DomainFixedSpaceMesh( &pdS[0], pdS.size() ) );
  
    // price surface
    m_pPriceSurface = make_ptr( new SurfaceGeneral(m_pDomain) );

    // constraint flag surface
    if ( m_bHasConstraintFlags )
      m_pFlagSurface = make_ptr( new SurfaceFlag(m_pDomain) );

    // vega surface (vega is computed by PDE)
    if ( m_computationalFlags.GetComputeVega() )
      m_pVegaSurface = make_ptr( new SurfaceGeneral(m_pDomain) );
  }

  // set the analysis time
  m_dAnalysisTime = params.GetAnalysisTime();
}

void BackwardNumOutput::SaveAnalysisData(BackwardInstData& instdata)
{    
  // Extra points are added for each Dirichlet boundary.  Linear 
  // extrapolation of analysis date prices will then work in all cases.
  // The new points will be just outside the grid. Adding
  // outside the grid could cause problems for plotting.  It may
  // also add negative spots if the barrier is near zero.  Adding
  // inside the grid is technically wrong.  Neither choice (inside or 
  // outside) is perfect.
  //
  // See also the surface code
  size_t nIdx;

  size_t nNbS = m_nNbS;

  // Check if extra points are being added at boundaries
  size_t nNbExtraLeft = 0;
  size_t nNbExtraRight = 0;
  GetNbExtraPoints(instdata, nNbExtraLeft, nNbExtraRight);

  // New vectors will have the adjusted size
  nNbS += nNbExtraLeft + nNbExtraRight;

  // The spots at analysis date, copied in case extra points are added
  m_pdAnalysisSpots = CopySpaceMesh(instdata);

  // Save the data.  Cannot just set in the output class since
  // the output class has not been created yet
  m_pdAnalysisValues.resize(nNbS);
  CopyData(instdata, instdata.m_pdPrices.Get(), &m_pdAnalysisValues[0]);

  // Use finite differences for theta.  For options, the old solution
  // should be valid, even with events, since the old solution is set
  // at the start of the current timestep.  
  m_pdAnalysisThetas.resize(nNbS);
  if ( nNbExtraLeft == 0 && nNbExtraRight == 0)
  {
    // Avoid copy if possible
    for (nIdx = 0; nIdx < m_nNbS; nIdx++)
      m_pdAnalysisThetas[nIdx] = - instdata.m_dInverseTimeStep
                  * (instdata.m_pdPrices[nIdx] - instdata.m_pdOldPrices[nIdx]);
  }
  else
  {
    // Copying the data is not the most efficient, but avoids duplicating
    // the CopyData code
    std::vector<double> pdTmp(nNbS);

    for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
      pdTmp[nIdx] = - instdata.m_dInverseTimeStep
          * (instdata.m_pdPrices[nIdx] - instdata.m_pdOldPrices[nIdx]);

    CopyData(instdata, &pdTmp[0], &m_pdAnalysisThetas[0]);
  }


  if ( m_computationalFlags.GetComputeVega() )
  {
    m_pdAnalysisVegas.resize(nNbS);
    CopyData(instdata, instdata.m_pdVegas.Get(), &m_pdAnalysisVegas[0]);
  }

  if( m_computationalFlags.GetComputeFugit() )
  {
    m_pdAnalysisFugits.resize(nNbS);
    CopyData(instdata, instdata.m_pdFugits.Get(), &m_pdAnalysisFugits[0]);
  }

  // if there has constraint flags
  if ( m_bHasConstraintFlags )
  {
    m_pAnalysisFlags.resize(nNbS);
    const int *piFlags = instdata.GetConstraintFlags();

    CopyData(instdata, piFlags, &m_pAnalysisFlags[0]);
  }

} 

void BackwardNumOutput::SaveSurfaceDataFrom
     (BackwardInstData& instdata, bool bEndOfGrid)
{
  // Extra points are added for each Dirichlet boundary.  Linear 
  // extrapolation of the surface prices should then work in all cases.
  // The new points will be just outside the grid. Adding
  // outside the grid could cause problems for plotting.  It may
  // also add negative spots if the barrier is near zero.  Adding
  // inside the grid is technically wrong.  Neither choice (inside or 
  // outside) is perfect.
  //
  // See also the analysis date code

  // check if extra points are needed
  size_t nNbExtraLeft = 0;
  size_t nNbExtraRight = 0;
  GetNbExtraPoints(instdata, nNbExtraLeft, nNbExtraRight);

  // The grid may have an adjusted size
  size_t nNbS = m_nNbS;
  nNbS += nNbExtraLeft + nNbExtraRight;

  // In all cases, convert the data to the right format before adding
  // to the surface.  Do the price first    
  numeric::SurfaceDouble::Doubles pdTmpData;
  pdTmpData.resize(nNbS);
  CopyData( instdata, &instdata.m_pdPrices[0], &pdTmpData[0] );
  m_pPriceSurface->Add(pdTmpData, bEndOfGrid);

  // if there has constraint flags surface
  if ( m_bHasConstraintFlags )
  {
    numeric::SurfaceFlag::Flags piTmpFlags(nNbS);
    const int *piFlags = instdata.GetConstraintFlags();
    CopyData( instdata, piFlags, &piTmpFlags[0] );
    m_pFlagSurface->Add(piTmpFlags, bEndOfGrid);
  }

  // Vega, if requested
  if ( m_computationalFlags.GetComputeVega() )
  {
    CopyData( instdata, instdata.m_pdVegas.Get(), &pdTmpData[0] ); 
    m_pVegaSurface->Add(pdTmpData, bEndOfGrid);
  }

  // fugit, if requested
  if( m_computationalFlags.GetComputeFugit() )
  {
    CopyData( instdata, instdata.m_pdFugits.Get(), &pdTmpData[0] ); 
    m_pFugitSurface->Add(pdTmpData, bEndOfGrid);
  }

}

void BackwardNumOutput::SaveSurfaceWithFixMesh
     (BackwardInstData& instdata, DomainFixedSpaceMesh& domain, double dTime)
{
  domain.AddTime(dTime);

  SaveSurfaceDataFrom(instdata);
} 

void BackwardNumOutput::SaveSurfaceAtEndOfGrid
     (BackwardInstData& /* instdata */, double /* dTime */)
{
  FAIL("SaveSurfaceAtEndOfGrid() must be implmented for specific NumOutput "
       "when using Special Engine.");

  /* exemple :
  m_pDomain->AddSpotsAtTime(pdS, instdata.m_meshes.GetTime(), true);

  SaveSurfaceFrom(instdata, true);
  */
}

void BackwardNumOutput::UpdateMeAtEndOfGrid
     (BackwardInstData& instdata, double dTime)
{
  m_pdS = instdata.GetSpaceMesh(m_nNbS);

  // do other thing than analysis data save, surface save and final save
  DoSpecialWorkWith(instdata, dTime);

  // don't need to Check if this is the analysis date, as this is end of grid

  // If surfaces are being computed, then save the data at this step
  // This function is called by init, so initial data is automatically stored.
  // Also remember to update the shared underlying domain!
  if ( m_computationalFlags.GetComputeSurface() )
    SaveSurfaceAtEndOfGrid(instdata, dTime);
} 

void BackwardNumOutput::DoSpecialWorkWith
     (BackwardInstData& /* instdata */, double /* dTime */)
{
}

void BackwardNumOutput::UpdateMe(BackwardInstData& instdata, double dTime)
{
  m_pdS = instdata.GetSpaceMesh(m_nNbS);

  // do other thing than analysis data save, surface save and final save
  DoSpecialWorkWith(instdata, dTime);

  // Check if this is the analysis date

  // when an event occurs, UpdateMe will be called twice
  // And since for now we want the the values after having applied the events,
  // we just leave the analysisdate values be updated twice
  if ( AreTimesEqual(m_dAnalysisTime, dTime) )
    SaveAnalysisData(instdata);

  // If surfaces are being computed, then save the data at this step
  // This function is called by init, so initial data is automatically stored.
  // Also remember to update the shared underlying domain!
  if ( m_computationalFlags.GetComputeSurface() )
    SaveSurface(instdata, dTime);
}

// Instead of copying data from tmp arrays, could also pass instdata
// to this function
void BackwardNumOutput::Finalize(BackwardInstData& instdata)
{
  CalculateFinalScalarResult(instdata);

  if (!m_bFinalSave)
    return;

  size_t nIdx;

  // The lifetime of the pricer objects may be less than the lifetime
  // of numoutput. So, we need to copy data that is needed for 
  // returning the model output, or data that is needed by other
  // internal routines

  // Copy the final grid referenced in instdata
  m_pdFinalSpots.resize(m_nNbS);
  for (nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdFinalSpots[nIdx] = m_pdS[nIdx];

  // Copy the final values stored in instdata
  m_pdFinalValues.resize(m_nNbS);
  for (nIdx = 0; nIdx < m_nNbS; nIdx++)
    m_pdFinalValues[nIdx] = instdata.m_pdPrices[nIdx];

  // Copy vega, if requested
  if ( m_computationalFlags.GetComputeVega() )
  {
    m_pdFinalVegas.resize(m_nNbS);
    for (nIdx = 0; nIdx < m_nNbS; nIdx++)
      m_pdFinalVegas[nIdx] = instdata.m_pdVegas[nIdx];
  }

}

void BackwardNumOutput::CalculateFinalScalarResult
     (BackwardInstData& instdata)
{
  double dInitialSpot = instdata.GetInitialSpot();

  double *pdTmpPrice = instdata.m_pdPrices.Get();

  // Interpolate for price at the spot
  Interpolate(m_pdS, instdata.m_pdPrices.Get(),
              m_nNbS, &dInitialSpot, &m_dPrice, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  std::vector<double> pdResults;
  pdResults.resize(m_nNbS);

  // Interpolate for delta at the spot
  // TODO: Only compute and interpolate over a subset of data near the spot
  
  ComputeDelta(m_pdS, instdata.m_pdPrices.Get(), m_nNbS, &pdResults[0]);
  
  Interpolate(m_pdS, &pdResults[0],
              m_nNbS, &dInitialSpot, &m_dDelta, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  // Interpolate for gamma at the spot  
  // TODO: Only compute and interpolate over a subset of data near the spot
  ComputeGammaFD(m_pdS, instdata.m_pdPrices.Get(), m_nNbS, &pdResults[0]);

  Interpolate(m_pdS, &pdResults[0],
              m_nNbS, &dInitialSpot, &m_dGamma, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);


  // Interpolate for vega at the spot (if requested)
  if ( m_computationalFlags.GetComputeVega() )
  {
    Interpolate(m_pdS, instdata.m_pdVegas.Get(),
                m_nNbS, &dInitialSpot, &m_dVega, 1, 
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  }

  // Interpolate for fugit at the spot (if requested)
  if ( m_computationalFlags.GetComputeFugit() )
  {
    Interpolate(m_pdS, instdata.m_pdFugits.Get(),
                m_nNbS, &dInitialSpot, &m_dFugit, 1, 
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  }

  // Interpolate for Theta at the spot

  // Take the negative, since we want deriv wrt current time, not
  // the valuation time
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
    pdResults[nIdx] = -instdata.m_dInverseTimeStep
              * (pdTmpPrice[nIdx] - instdata.m_pdOldPrices[nIdx]);

  Interpolate(m_pdS, &pdResults[0],
              m_nNbS, &dInitialSpot, &m_dTheta, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  
  // Value after default
  m_dValueAfterDefault = instdata.GetValueAfterDefaultAtInitialSpot();
}

void BackwardNumOutput::GetOutput(ito33::finance::ModelOutput* pOutput)
{
  pOutput->SetPrice(m_dPrice);
  pOutput->SetDelta(m_dDelta);
  pOutput->SetGamma(m_dGamma);
  pOutput->SetTheta(m_dTheta);
  pOutput->SetValueAfterDefault(m_dValueAfterDefault);

  if ( m_computationalFlags.GetComputeVega() )
    pOutput->SetVega(m_dVega);

  if ( m_computationalFlags.GetComputeFugit() )
    pOutput->SetFugit(m_dFugit);

  // Set the analysis date data, if available
  if ( m_pdAnalysisValues.size() > 0 )
  {
    pOutput->SetSpotsAtAnalysisDate(m_pdAnalysisSpots);
    pOutput->SetPricesAtAnalysisDate(m_pdAnalysisValues);

    size_t nNbS = m_pdAnalysisSpots.size();

    std::vector<double> pdAnalysisDeltas(nNbS);
    
    ComputeDelta(&m_pdAnalysisSpots[0], &m_pdAnalysisValues[0], nNbS, 
                 &pdAnalysisDeltas[0]);

    pOutput->SetDeltasAtAnalysisDate(pdAnalysisDeltas);

    std::vector<double> pdAnalysisGammas(nNbS);
    
    ComputeGammaFD(&m_pdAnalysisSpots[0], &m_pdAnalysisValues[0], nNbS,
                   &pdAnalysisGammas[0]);

    pOutput->SetGammasAtAnalysisDate(pdAnalysisGammas);

    pOutput->SetThetasAtAnalysisDate(m_pdAnalysisThetas);

    if ( m_computationalFlags.GetComputeVega() )
      pOutput->SetVegasAtAnalysisDate(m_pdAnalysisVegas);

    if ( m_computationalFlags.GetComputeFugit() )
      pOutput->SetFugitsAtAnalysisDate(m_pdAnalysisFugits);
  }

  // If surfaces are needed, set them.
  if ( m_computationalFlags.GetComputeSurface() )
  {
    // Only save a surface with 2 or more time points.
    // (continuous path dep problems will only save at valuation date)  
    size_t nNbTimes = m_pDomain->GetTimes().size();
    if ( nNbTimes > 1 )
    {
      // Update the domain before creating the user output
      m_pDomain->GenerateOutputDates();

      pOutput->SetDomain(m_pDomain);

      // Need to convert from numeric price surface to finance price surface 
      // before setting. The former is not derived from the later
      shared_ptr<finance::SurfaceDouble> 
      pTmpSurface(new finance::SurfaceDouble(m_pPriceSurface));
      pOutput->SetPriceSurface(pTmpSurface);

      shared_ptr<SurfaceDouble> pDeltaSurface( new SurfaceGeneral(m_pDomain) );
      shared_ptr<SurfaceDouble> pGammaSurface( new SurfaceGeneral(m_pDomain) );
      shared_ptr<SurfaceDouble> pThetaSurface( new SurfaceGeneral(m_pDomain) );

      m_pPriceSurface->GetDeltaAndGamma(pDeltaSurface, pGammaSurface);
      m_pPriceSurface->GetThetaBackwardOnly(pThetaSurface);

      pTmpSurface = make_ptr( new finance::SurfaceDouble(pDeltaSurface) );
      pOutput->SetDeltaSurface(pTmpSurface);

      pTmpSurface = make_ptr( new finance::SurfaceDouble(pGammaSurface) );
      pOutput->SetGammaSurface(pTmpSurface);

      pTmpSurface = make_ptr( new finance::SurfaceDouble(pThetaSurface) );
      pOutput->SetThetaSurface(pTmpSurface);

      if ( m_computationalFlags.GetComputeVega() )
      {
        pTmpSurface = make_ptr( new finance::SurfaceDouble(m_pVegaSurface) );
        pOutput->SetVegaSurface(pTmpSurface);
      }

      if ( m_computationalFlags.GetComputeFugit() )
      {
        pTmpSurface = make_ptr( new finance::SurfaceDouble(m_pFugitSurface) );
        pOutput->SetFugitSurface(pTmpSurface);
      }
    } // if more than one point
  } // if saving surface
}

std::vector<double> 
BackwardNumOutput::CopySpaceMesh(const BackwardInstData& instdata)
{
  // check if extra points are needed
  size_t nNbExtraLeft = 0;
  size_t nNbExtraRight = 0;
  GetNbExtraPoints(instdata, nNbExtraLeft, nNbExtraRight);

  // The grid may have an adjusted size
  size_t nNbS = m_nNbS;
  nNbS += nNbExtraLeft + nNbExtraRight;

  // Copy the spots. Create return vector with appropriate size.
  std::vector<double> pdNewSpots(nNbS);
  for (size_t nIdxS = nNbExtraLeft; nIdxS < m_nNbS + nNbExtraLeft; nIdxS++)
    pdNewSpots[nIdxS] = instdata.m_pdS[nIdxS - nNbExtraLeft];

  // Add the extra points, if needed
  if (nNbExtraLeft > 0)
  {
    // try to shift by some percentage of the grid value. Do not
    // shift less than the double tolerance
    double dShift = 100.0 * DOUBLETOLERANCE * instdata.m_pdS[0];
    if ( dShift < DOUBLETOLERANCE )
      dShift = 10.0 * DOUBLETOLERANCE;

    // Add points just to left of the Dirichlet node
    for (size_t nIdx = 0; nIdx < nNbExtraLeft; nIdx++)
      pdNewSpots[nIdx] = instdata.m_pdS[0] - (nNbExtraLeft - nIdx) * dShift;
  }

  if (nNbExtraRight > 0)
  {
    // try to shift by some percentage of the grid value. Do not
    // shift less than the double tolerance
    double dShift = 100.0 * DOUBLETOLERANCE * instdata.m_pdS[m_nNbS - 1];
    if ( dShift < DOUBLETOLERANCE )
      dShift = 10.0 * DOUBLETOLERANCE;

    // Add points just to right of the Dirichlet node
    for (size_t nIdx = 0; nIdx < nNbExtraRight; nIdx++)
      pdNewSpots[nNbS - nIdx - 1] = instdata.m_pdS[m_nNbS - 1] 
                                  + (nNbExtraRight - nIdx) * dShift;
  }

  return pdNewSpots;
}


template <typename T>
void BackwardNumOutput::CopyData(const BackwardInstData& instdata,
                                 const T* pdSource, 
                                 T* pdTarget)
{

  // check if extra points are needed
  size_t nNbExtraLeft = 0;
  size_t nNbExtraRight = 0;
  GetNbExtraPoints(instdata, nNbExtraLeft, nNbExtraRight);

  size_t nIdxTarget = 0;

  // add extra points, if needed
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbExtraLeft; nIdx++)
    pdTarget[nIdxTarget++] = pdSource[0];

  // copy the original points
  for (nIdx = 0; nIdx < m_nNbS; nIdx++)
    pdTarget[nIdxTarget++] = pdSource[nIdx];

  // add extra points, if needed
  for (nIdx = 0; nIdx < nNbExtraRight; nIdx++)
    pdTarget[nIdxTarget++] = pdSource[m_nNbS - 1];

}

void BackwardNumOutput::GetNbExtraPoints(const BackwardInstData& instdata,
                                         size_t& nNbLeft, size_t& nNbRight)
{
  // Add two extra points so that delta and gamma extrapolation will be
  // correct.  Only one point is needed for the price.  However, delta
  // and gamma values at barriers are taken to be interior limits, so we
  // need two extra points with zero values to ensure extrapolated values
  // are zero.
  if ( instdata.GetBoundaryCondition().GetLeftType() == BCType_Dirichlet )
    nNbLeft = 2;
  else
    nNbLeft = 0;

  if ( instdata.GetBoundaryCondition().GetRightType() == BCType_Dirichlet )
    nNbRight = 2;
  else
    nNbRight = 0;
}
