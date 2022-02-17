/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/backwardnumoutput.cpp
// Purpose:     implementation of BackwardNumOutput class 
// Created:     2005/01/13
// RCS-ID:      $Id: backwardnumoutput.cpp,v 1.31 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"
#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/params.h"

#include "ito33/hg/modeloutput.h"
#include "ito33/hg/multioutput.h"

#include "hg/model.h"
#include "hg/backwardnumoutput.h"
#include "hg/backwardinstdata.h"
#include "hg/sensitivitymethod.h"

// implement the AutoPtrDeleter for BackwardNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::BackwardNumOutput);
}

namespace ito33
{

namespace hg
{
  using namespace numeric;


void BackwardNumOutput::Init(BackwardInstData& instdata)
{
  m_nNbRegimes = instdata.m_nNbRegimes;
  
  // Get the space mesh. Also sets m_nNbS.
  m_pdS = instdata.GetSpaceMesh(m_nNbS);
  m_nNbX = m_nNbS * m_nNbRegimes;

  // If surfaces are requested, construct them. The first step is constructing
  // the underlying domain which is shared between the surfaces.
  if ( m_computationalFlags.GetComputeSurface() )
  {

    // copy the mesh in case extra boundary points need to be added
    std::vector<double> pdS = CopySpaceMesh(instdata);

    m_pDomain = make_ptr( new DomainFixedSpaceMesh( &pdS[0], pdS.size() ) );
    
    // price surfaces
    m_ppPriceSurfaces.resize(m_nNbRegimes);
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_ppPriceSurfaces[nIdxR] = make_ptr( new SurfaceGeneral(m_pDomain) );

    // constraint flag surface
    if ( m_bHasConstraintFlags )
      m_pFlagSurface = make_ptr( new SurfaceFlag(m_pDomain) );

  }

  // set the analysis time
  m_dAnalysisTime = m_params.GetAnalysisTime();

  // Set sensitivity flags
  m_computationalFlags.SetSensitivityFlags( instdata.m_pbComputeSensitivities );

  // Clear the recovery value vector so we can safely push_back
  m_pdRecoveryValues.clear();

  // Clear the times so we can push_back
  m_pdTimes.clear();

  // initialize m_dObjectif to negative value which means that it's not
  // computed by default.
  m_dObjectif = -1;
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
  size_t nNbX = m_nNbX;

  // Check if extra points are being added at boundaries
  size_t nNbExtraLeft = 0;
  size_t nNbExtraRight = 0;
  GetNbExtraPoints(instdata, nNbExtraLeft, nNbExtraRight);

  // New vectors will have the adjusted size
  nNbX += nNbExtraLeft * m_nNbRegimes + nNbExtraRight * m_nNbRegimes;

  // The spots at analysis date, copied in case extra points are added
  m_pdAnalysisSpots = CopySpaceMesh(instdata);

  // Save the data.  Cannot just set in the output class since
  // the output class has not been created yet. And we may store only the
  // price at the first regime in the output class  
  m_pdAnalysisValues.resize(nNbX);
  CopyData(instdata, instdata.m_pdPrices.Get(), &m_pdAnalysisValues[0], 
           m_nNbRegimes); 

  // Use finite differences for theta.  The old solution should be valid, 
  // even with events, since the old solution is set at the start of the 
  // current timestep. 
  m_pdAnalysisThetas.resize(nNbX);
  if ( nNbExtraLeft == 0 && nNbExtraRight == 0)
  {
    // Avoid copy if possible
    for (size_t nIdx = 0; nIdx < nNbX; nIdx++)
      m_pdAnalysisThetas[nIdx] = - instdata.m_dInverseTimeStep
          * (instdata.m_pdPrices[nIdx] - instdata.m_pdOldPrices[nIdx]);
  }
  else
  {
    // Copying the data is not the most efficient, but avoids duplicating
    // the CopyData code
    std::vector<double> pdTmp(m_nNbX);

    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      pdTmp[nIdx] = - instdata.m_dInverseTimeStep
          * (instdata.m_pdPrices[nIdx] - instdata.m_pdOldPrices[nIdx]);

    CopyData(instdata, &pdTmp[0], &m_pdAnalysisThetas[0], m_nNbRegimes);
  }

  // if fugit is computed
  if ( m_computationalFlags.GetComputeFugit() )
  {
    m_pdAnalysisFugits.resize(nNbX);
    CopyData(instdata, instdata.m_pdFugits.Get(), &m_pdAnalysisFugits[0], 
             m_nNbRegimes);
  }

  // if there are constraint flags
  if ( m_bHasConstraintFlags )
  {
    m_pAnalysisFlags.resize(nNbX);
    CopyData(instdata, instdata.GetConstraintFlags(), &m_pAnalysisFlags[0], 
             m_nNbRegimes);
  }

} 

void BackwardNumOutput::SaveSurfaceDataFrom
     (BackwardInstData& instdata, bool bEndOfGrid)
{
  // The space mesh should already be updated

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

  // For now, save the prices in different regimes in different
  // surfaces.  For other values (constraints, fugit), only save 
  // values from the first regime.

  // In all cases, convert the data to the right format before adding
  // to the surface.  Do the price first      
  numeric::SurfaceDouble::Doubles pdTmpData;

  pdTmpData.resize(nNbS);
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    size_t nOffset = nIdxR * m_nNbS;

    CopyData( instdata, &instdata.m_pdPrices[nOffset], &pdTmpData[0], 1);

    m_ppPriceSurfaces[nIdxR]->Add(pdTmpData, bEndOfGrid);
  }

  // if there has constraint flags surface
  if ( m_bHasConstraintFlags )
  {
    numeric::SurfaceFlag::Flags piTmpFlags(nNbS);

    const int *piFlags = instdata.GetConstraintFlags();
    CopyData( instdata, piFlags, &piTmpFlags[0], 1);

    m_pFlagSurface->Add(piTmpFlags, bEndOfGrid);
  }

  // fugit, if requested
  if( m_computationalFlags.GetComputeFugit() )
  {
    CopyData( instdata, instdata.m_pdFugits.Get(), &pdTmpData[0], 1); 

    m_pFugitSurface->Add(pdTmpData, bEndOfGrid);
  }
}

void BackwardNumOutput::SaveSurface
     (BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainFixedSpaceMesh*>(m_pDomain.get()), 
    "The type of Domain in BackwardNumOutput should be DomainFixedSpaceMesh."
  );

  DomainFixedSpaceMesh& 
    domainFixed = static_cast<DomainFixedSpaceMesh&>(*m_pDomain);
  
  domainFixed.AddTime(dTime);

  SaveSurfaceDataFrom(instdata);
}

void BackwardNumOutput::SaveSurfaceAtEndOfGrid
     (BackwardInstData& /* instdata */, double /* dTime */)
{
  FAIL("SaveSurfaceAtEndOfGrid() must be implemented for specific NumOutput "
       "when using Special Engine.");
}

void BackwardNumOutput::UpdateMeAtEndOfGrid
     (BackwardInstData& instdata, double dTime)
{
  m_pdS = instdata.GetSpaceMesh(m_nNbS);
  m_nNbX = m_nNbS * m_nNbRegimes;

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
     (BackwardInstData& instdata, double /* dTime */)
{
  if ( instdata.m_aData.m_bIsValid )
  {
    m_pAdjointDatas.push_back(instdata.m_aData);

    // The current implementation will add the interpolation matrix to
    // the step of the previous step having event. Correct it manually
    size_t nNbDatas = m_pAdjointDatas.size();
    if ( nNbDatas > 1 )
      m_pAdjointDatas[nNbDatas - 2].m_pInterpMatrix 
           = m_pAdjointDatas[nNbDatas - 1].m_pInterpMatrix;
  }
}

void BackwardNumOutput::UpdateMe(BackwardInstData& instdata, double dTime)
{
  m_pdS = instdata.GetSpaceMesh(m_nNbS);
  m_nNbX = m_nNbS * m_nNbRegimes;

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

// Instead of copying data from tmp arrays, could also pass instdata
// to this function
void BackwardNumOutput::Finalize(BackwardInstData& instdata)
{
  m_pdS = instdata.GetSpaceMesh(m_nNbS);
  m_nNbX = m_nNbS * m_nNbRegimes;

  CalculateFinalScalarResult(instdata);

  if ( instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
    SensitivityByAdjoint(instdata);

  // Sensitivity w.r.t post default volatility will just be zero
  if (    m_computationalFlags.HasSensitivityFlags()
       && m_computationalFlags.GetSensitivityFlags().back() )
    m_pdSensitivities.push_back(0);

  if (!m_bFinalSave)
    return;

  // The lifetime of the pricer objects may be less than the lifetime
  // of numoutput. So, we need to copy data that is needed for 
  // returning the model output, or data that is needed by other
  // internal routines

  // Copy the final grid referenced in instdata
  m_pdFinalSpots.resize(m_nNbS);
  for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    m_pdFinalSpots[nIdxS] = m_pdS[nIdxS];

  // Copy the final values stored in instdata
  m_pdFinalValues.resize(m_nNbX);
  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdFinalValues[nIdx] = instdata.m_pdPrices[nIdx];

  // we don't care about fugit now, maybe later

}

void BackwardNumOutput::CalculateFinalScalarResult
     (BackwardInstData& instdata)
{
  // Only the scalar values at the first regime is computed.
  const double dInitialSpot = instdata.GetInitialSpot();
  const double *pdPrices = instdata.m_pdPrices.Get();

  AutoPtr<InterpolationMatrix>
    pMatrix( new InterpolationMatrix(&dInitialSpot, 1, m_pdS, m_nNbS, 1) );
  
  // Interpolate for price at the spot
  pMatrix->ProductMatrixVector(pdPrices, &m_dPrice);

  m_dValueAfterDefault = instdata.GetValueAfterDefaultAtInitialSpot();

  std::vector<double> pdResults(m_nNbS);
  
  ComputeDelta(m_pdS, pdPrices, m_nNbS, &pdResults[0]);

  // Interpolate for delta at the spot
  pMatrix->ProductMatrixVector(&pdResults[0], &m_dDelta);

  ComputeGammaFD(m_pdS, pdPrices, m_nNbS, &pdResults[0]);

  // Interpolate for gamma at the spot  
  pMatrix->ProductMatrixVector(&pdResults[0], &m_dGamma);

  // Interpolate for fugit at the spot (if requested)
  if ( m_computationalFlags.GetComputeFugit() )
    pMatrix->ProductMatrixVector(instdata.m_pdFugits.Get(), &m_dFugit);

  // Interpolate for Theta at the spot

  // Take the negative, since we want deriv wrt current time, not
  // the valuation time
  for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    pdResults[nIdxS] = - instdata.m_dInverseTimeStep
                     * (pdPrices[nIdxS] - instdata.m_pdOldPrices[nIdxS]);

  pMatrix->ProductMatrixVector(&pdResults[0], &m_dTheta);

  // Interpolate for sensitivities at the spot (if requested)
  if ( instdata.m_sensitivityMethod == SensitivityMethod_PDE )
  {
    size_t nNbSensitivities = instdata.m_pSensitivityData.size();
    m_pdSensitivities.resize( nNbSensitivities );

    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
      pMatrix->ProductMatrixVector
               ( instdata.m_ppdSensitivities[nIdxD].Get(), 
                 &m_pdSensitivities[nIdxD] );
  }

  const std::vector<double>& pdSpots = m_params.GetSpots();
  size_t nNbSpots = pdSpots.size();
  if ( nNbSpots )
  {
    pMatrix = AutoPtr<InterpolationMatrix>
        ( new InterpolationMatrix(&pdSpots[0], nNbSpots, m_pdS, m_nNbS, 1) );

    m_pdPricesAtObservations.resize(nNbSpots);

    pMatrix->ProductMatrixVector(instdata.m_pdPrices.Get(),
                                 &m_pdPricesAtObservations[0]); 
  }
  
  if ( instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
  {
    ASSERT( !m_pAdjointDatas.empty() );

    Array<double> pdResiduals;

    if ( nNbSpots )
    {
      m_bSensitivityOnObjectif = true;
    
      const std::vector<double>& pdMarketPrices = m_params.GetMarketPrices();
      
      ASSERT( pdMarketPrices.size() == nNbSpots );

      m_dObjectif = 0.;
      pdResiduals = Array<double>(nNbSpots);
      for (size_t nIdxObs = 0; nIdxObs < nNbSpots; nIdxObs++)
      {
        pdResiduals[nIdxObs] = m_pdPricesAtObservations[nIdxObs]
                             - pdMarketPrices[nIdxObs];
        m_dObjectif += 0.5 * pdResiduals[nIdxObs] * pdResiduals[nIdxObs];
      }
    }
    else
    {
      pdResiduals = Array<double>(1);
      pdResiduals[0] = m_dPrice;
    }

    m_pAdjointDatas.back().m_pRMatrix = pMatrix;
    m_pAdjointDatas.back().m_pdResiduals = pdResiduals;
  }
}

void BackwardNumOutput::GetOutput(finance::ModelOutput& output)
{
  output.SetPrice(m_dPrice);
  output.SetValueAfterDefault(m_dValueAfterDefault);

  output.SetDelta(m_dDelta);
  output.SetGamma(m_dGamma);
  output.SetTheta(m_dTheta);

  if ( m_computationalFlags.HasSensitivityFlags() )
    SetSensitivities(m_pdSensitivities);

/*  
  if ( m_computationalFlags.GetComputeFugit() )
    output.SetFugit(m_dFugit);
*/
  // Set the analysis date data, if available
  if ( m_pdAnalysisValues.size() > 0 )
  {
    size_t nNbS = m_pdAnalysisSpots.size();

    output.SetSpotsAtAnalysisDate(m_pdAnalysisSpots);
    
    // Only save one regime in the finance output
    if (  m_pdAnalysisValues.size() > nNbS )
    {
      std::vector<double> pdAnalysisValues(nNbS);

      for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
        pdAnalysisValues[nIdx] = m_pdAnalysisValues[nIdx];

      output.SetPricesAtAnalysisDate(pdAnalysisValues);
    }
    else
      output.SetPricesAtAnalysisDate(m_pdAnalysisValues);

    
    std::vector<double> pdAnalysisDeltas(nNbS);
    
    ComputeDelta(&m_pdAnalysisSpots[0], &m_pdAnalysisValues[0], nNbS, 
                 &pdAnalysisDeltas[0]);

    output.SetDeltasAtAnalysisDate(pdAnalysisDeltas);

    std::vector<double> pdAnalysisGammas(nNbS);
    
    ComputeGammaFD(&m_pdAnalysisSpots[0], &m_pdAnalysisValues[0], nNbS,
                   &pdAnalysisGammas[0]);

    output.SetGammasAtAnalysisDate(pdAnalysisGammas);

    output.SetThetasAtAnalysisDate(m_pdAnalysisThetas);
/*
    if ( m_computationalFlags.GetComputeFugit() )
      output.SetFugitsAtAnalysisDate(m_pdAnalysisFugits);

    if ( m_bHasConstraintFlags )
      output.SetConstraintFlagsAtAnalysisDate(m_pAnalysisFlags);
*/
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

      output.SetDomain(m_pDomain);

      // Need to convert from numeric price surface to finance price surface 
      // before setting. The former is not derived from the later
      shared_ptr<finance::SurfaceDouble> 
        pTmpSurface(new finance::SurfaceDouble(m_ppPriceSurfaces[0]));
      output.SetPriceSurface(pTmpSurface);

      shared_ptr<SurfaceDouble> pDeltaSurface( new SurfaceGeneral(m_pDomain) );
      shared_ptr<SurfaceDouble> pGammaSurface( new SurfaceGeneral(m_pDomain) );
      shared_ptr<SurfaceDouble> pThetaSurface( new SurfaceGeneral(m_pDomain) );

      m_ppPriceSurfaces[0]->GetDeltaAndGamma(pDeltaSurface, pGammaSurface);
      m_ppPriceSurfaces[0]->GetThetaBackwardOnly(pThetaSurface);

      pTmpSurface = make_ptr( new finance::SurfaceDouble(pDeltaSurface) );
      output.SetDeltaSurface(pTmpSurface);

      pTmpSurface = make_ptr( new finance::SurfaceDouble(pGammaSurface) );
      output.SetGammaSurface(pTmpSurface);

      pTmpSurface = make_ptr( new finance::SurfaceDouble(pThetaSurface) );
      output.SetThetaSurface(pTmpSurface);
/*    
      if ( m_bHasConstraintFlags )
      { 
        shared_ptr<finance::SurfaceFlag> 
          pTmpFlagSurface(new finance::SurfaceFlag(m_pFlagSurface));
        output.SetConstraintFlagSurface(pTmpFlagSurface);
      }

      if ( m_computationalFlags.GetComputeFugit() )
      {
        pTmpSurface = new finance::SurfaceDouble(m_pFugitSurface);
        output.SetFugitSurface(pTmpSurface);
      }
  */
    } // if more than one surface point
  } // saving surfaces
  
}

shared_ptr<MultiOutput> BackwardNumOutput::GetMultiOutput()
{
  shared_ptr<MultiOutput> pOutput(new MultiOutput);
  
  GetOutput( *pOutput );

  pOutput->SetPrices(m_pdPricesAtObservations);

  if ( m_dObjectif >= 0 )
  {
    // Update the objective function appropriately    
    pOutput->SetObjectif(m_dObjectif);
  }

  return pOutput;
}

shared_ptr<ModelOutput> BackwardNumOutput::GetModelOutput()
{
  shared_ptr<ModelOutput> pOutput(new ModelOutput);

  GetOutput(*pOutput);

  return pOutput;
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
                                 T* pdTarget, 
                                 size_t nNbRegimes)
{

  ASSERT_MSG(nNbRegimes <= m_nNbRegimes,
    "Invalid number of regimes in call to CopyData");

  // check if extra points are needed
  size_t nNbExtraLeft = 0;
  size_t nNbExtraRight = 0;
  GetNbExtraPoints(instdata, nNbExtraLeft, nNbExtraRight);

  size_t nIdxTarget = 0;
  for (size_t nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    size_t nOffset = nIdxR * m_nNbS;

    // add extra points, if needed
    size_t nIdx;
    for (nIdx = 0; nIdx < nNbExtraLeft; nIdx++)
      pdTarget[nIdxTarget++] = pdSource[nOffset];

    // copy the original points
    for (nIdx = 0; nIdx < m_nNbS; nIdx++)
      pdTarget[nIdxTarget++] = pdSource[nOffset + nIdx];

    // add extra points, if needed
    for (nIdx = 0; nIdx < nNbExtraRight; nIdx++)
      pdTarget[nIdxTarget++] = pdSource[nOffset + m_nNbS - 1];

  } // loop over regimes

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


void BackwardNumOutput::SubractNumOutput(const BackwardNumOutput& otherNumOutput)
{
  // Subtract otherNumOutput from this BackwardNumOutput
  size_t nIdx;

  // scalar data
  m_dPrice -= otherNumOutput.m_dPrice;
  m_dDelta -= otherNumOutput.m_dDelta;
  m_dGamma -= otherNumOutput.m_dGamma;
  m_dTheta -= otherNumOutput.m_dTheta;
  m_dFugit = std::min(m_dFugit, otherNumOutput.m_dFugit);
  for (nIdx = 0; nIdx < m_pdSensitivities.size(); nIdx++)
    m_pdSensitivities[nIdx] -= otherNumOutput.m_pdSensitivities[nIdx];

  // analysis date data
  if (m_pdAnalysisValues.size() > 0)
  {
    size_t nNbOther = otherNumOutput.m_pdAnalysisSpots.size();
    size_t nNbS = m_pdAnalysisSpots.size();
    std::vector<double> pdOtherAnalysisData(nNbOther);
    numeric::QuadraticInterpolate(&otherNumOutput.m_pdAnalysisSpots[0],
                                  &otherNumOutput.m_pdAnalysisValues[0],
                                  nNbOther,
                                  &m_pdAnalysisSpots[0],
                                  &pdOtherAnalysisData[0],
                                  nNbS);
                                  
    for (nIdx = 0; nIdx < nNbS; nIdx++)
      m_pdAnalysisValues[nIdx] -= pdOtherAnalysisData[nIdx];
  
    numeric::QuadraticInterpolate(&otherNumOutput.m_pdAnalysisSpots[0],
                                  &otherNumOutput.m_pdAnalysisThetas[0],
                                  nNbOther,
                                  &m_pdAnalysisSpots[0],
                                  &pdOtherAnalysisData[0],
                                  nNbS);

    for (nIdx = 0; nIdx < nNbS; nIdx++)
      m_pdAnalysisThetas[nIdx] -= pdOtherAnalysisData[nIdx];

  }

  // surface data
  // TODO
}




} // namespace hg

} // namespace ito33
