/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/cbmeshmanager.cpp
// Purpose:     cb multi-grid mesh manager for backward PDE problems
// Author:      Nabil
// Created:     2004/03/17
// RCS-ID:      $Id: cbmeshmanager.cpp,v 1.58 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/constants.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/numparams.h"

#include "ito33/numeric/mesh/cbmesh.h"
#include "ito33/numeric/mesh/genraf.h"
#include "ito33/numeric/mesh/gridwithcv.h"
#include "ito33/numeric/mesh/roots.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cbevent.h"
#include "ito33/pricing/model.h"
#include "ito33/pricing/cashflows.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"

namespace ito33
{

using namespace numeric;
using namespace numeric::mesh;

namespace pricing
{

struct GridInfo
{
  bool
    m_bHasOID,
    m_bHasCall,
    m_bHasConversion;

  size_t
    m_nIdxLastCall,
    m_nIdxLastConversion,
    m_nIdxNextCoupon;
};

void CBMeshManager::ConstructSpaceMesh()
{
  bool
    bMeetGridNode;

  size_t
    nIdxBeginNow,
    nIdxTime,
    nNbEvents;

  double
    dDateEventOld,
    dDateEventNow;
  
  std::vector<GridInfo> pGridInfo;
  GridInfo gridInfoCurrent;

  const std::list< shared_ptr<CBEvent> >
    &cbEvents = m_cbParams.GetAllCBEvents();


  
  //cbEvents doesn't contain the pricing date and the maturity date.
  //Then, we add 2 to the size of cbEvents in order to allocate a sufficient 
  //memory. 
  nNbEvents = cbEvents.size() + 2;
      
  m_pGrids.resize(nNbEvents);
    
  m_nNbGrids = 0;

  nIdxBeginNow = 0;
  
  dDateEventOld = m_cbParams.GetValuationTime();
  
  std::list< shared_ptr<CBEvent> >::const_iterator 
    itercbEvents;
     
  gridInfoCurrent.m_bHasCall = false;
  gridInfoCurrent.m_bHasConversion = false;
  gridInfoCurrent.m_bHasOID = false;
  gridInfoCurrent.m_nIdxLastCall = 0;
  gridInfoCurrent.m_nIdxLastConversion = 0;
  gridInfoCurrent.m_nIdxNextCoupon = 0;

  nIdxTime = 0;
  for(itercbEvents = cbEvents.begin(); itercbEvents != cbEvents.end();
      ++itercbEvents)
  {
    for(; 
           nIdxTime < m_nNbTimes
        && IsBefore(m_pdTimes[nIdxTime], (*itercbEvents)->m_dTime); 
        nIdxTime++)
      ;
   
    bMeetGridNode = false;
    switch( (*itercbEvents)->m_cbeventType )
    {
      case CBET_StartConversion:
      case CBET_EndConversion:
      case CBET_StartCall:
      case CBET_EndCall:
      case CBET_MonoDateCall:
      case CBET_MonoDateConversion:
      case CBET_Maturity:
      case CBET_ChangeOfRoot:
      case CBET_StartAccretion:
      case CBET_PayCoupon:
      case CBET_CBOptionMaturity:
        bMeetGridNode = true;
        break;
      case CBET_StartOfYear:
      case CBET_PayDividend:
        if( m_cbParams.HasNewShare() )
          bMeetGridNode = true;
        break;
      default:
        bMeetGridNode = false;
    }

    if(bMeetGridNode)
    {  
      dDateEventNow = (*itercbEvents)->m_dTime;

      if( !AreTimesEqual(dDateEventOld, dDateEventNow) )
      {
        // We find the end of the sub grid   
        m_pGrids[m_nNbGrids].m_pdTimes = m_pdTimes.Get() + nIdxBeginNow;         
        m_pGrids[m_nNbGrids].m_nNbTimes        = nIdxTime - nIdxBeginNow + 1;
        m_pGrids[m_nNbGrids].m_nIdxTimeBegin   = nIdxBeginNow;

        nIdxBeginNow                           = nIdxTime;
        dDateEventOld                          = dDateEventNow;
        m_nNbGrids++;
      }
/*
      GridInfo& gridInfo = pGridInfo[m_nNbGrids - 1];
      switch( (*itercbEvents)->m_cbeventType )
      {
        case CBET_EndConversion:
          gridInfo.m_bHasConversion = false;
          m_cbParams.DisableConversion();
          break;
        case CBET_StartConversion:
          m_cbParams.SetConversionIndex( (*itercbEvents)->m_iIndex );
          break;
        case CBET_EndCall:
          m_cbParams.DisableCall();
          break;
        case CBET_StartCall:
          m_cbParams.SetCallIndex( (*itercbEvents)->m_iIndex );
          break;
      } 
*/
    } // bMeetGridNode
  }

  m_pGrids.resize(m_nNbGrids);

  if( m_cbParams.GetMeshParams()->GetUniformSpaceGrid() )
  {
    // TODO:
  }
  else
    AdaptativeGrid();

  //Some initialization.
  m_nNbSMax = m_pGrids[0].m_nNbS;
  for(size_t nIdx = 1; nIdx < m_nNbGrids; nIdx++)
    if(m_pGrids[nIdx].m_nNbS > m_nNbSMax)
      m_nNbSMax = m_pGrids[nIdx].m_nNbS;
  
  m_pdS = Array<double>(m_nNbSMax);
  m_pdLogS = Array<double>(m_nNbSMax); 
  
  m_pdOldS = Array<double>(m_nNbSMax);
  m_pdOldLogS = Array<double>(m_nNbSMax);
}

void CBMeshManager::ComputeRootsAndUpdateGrids()
{

# ifndef NDEBUG
  {
    m_cbParams.m_bStateFunctionShouldBeCalledWithTimeArgument = true;
  }
# endif

  std::vector<GridWithCV>::iterator pGrid;

  double
    dTime;

  size_t
    nI,
    nIdxTime;

  for(pGrid = m_pGrids.begin(); pGrid != m_pGrids.end(); ++pGrid)
  {
    // reference to the grid that we will work on
    GridWithCV& grid = *pGrid;
    
    dTime = grid.m_pdTimes[0];

    bool 
      bHasCall = m_cbParams.HasContinousCallAt(dTime, true),
      bHasConversion = m_cbParams.HasContinousConversionAt(dTime, true);  
   
    // only when there is something special that we need to do a 
    // change of variable, actually not needed since the ctor sets it to false
    // but we do it here anyway 
    grid.m_bChangeOfVar = false;

    bool
      bHasRoots;

    bHasRoots
      = bHasCall && bHasConversion 
          // without notice
      && !m_cbParams.GetCalls()->HasNoticePeriod()
          // without make whole
      && (
          m_cbParams.GetCalls()->GetTriggerRate(m_cbParams.GetIndexCall()) > 0
          || !m_cbParams.GetCalls()->HasMakeWhole()
         );

    if (bHasCall && bHasConversion) // have some roots to be added to mesh
    {
      size_t
        nIdxRealRoot,
        nNbRoots;
      
      Root
        pRoots[10],
        pRootsOld[10];
      
      bool bGridNeedChangeOfVar = true;
      
      /*
         we will suppose that we should always do change of variable 
      */

      // There is a fixed boundary.
      // Verify if we have to allocate the roots table.
      // TODO: miIndiceDividendeAnnee ????
      if( //(m_cbParams.GetIsRealNewShare() && miIndiceDividendeAnnee >= 0) || 
         bGridNeedChangeOfVar )
      {
        grid.m_pdCorrectionS = new double [grid.m_nNbTimes + 2];

        // for now, the array stores fixed boundary at each time point
        double *pdRootsChangeVar = grid.m_pdCorrectionS;
        
        ComputeRoots(grid.m_pdTimes[0], nNbRoots, pRootsOld);
        
        const bool bMinus = false;

        ComputeRoots(grid.m_pdTimes[1], nNbRoots, pRoots, bMinus);

        for(nIdxRealRoot = 0; nIdxRealRoot < nNbRoots; nIdxRealRoot++)
        {
          if(pRootsOld[nIdxRealRoot].m_iRootType == RootType_Real && 
             pRoots[nIdxRealRoot].m_iRootType == RootType_Real)
            break;
        }
        pdRootsChangeVar[0] = pRootsOld[nIdxRealRoot].m_dRoot;
        pdRootsChangeVar[1] = pRoots[nIdxRealRoot].m_dRoot;

        // find the fix boundary at the biginning
        grid.m_dFixedBoundary = pdRootsChangeVar[0];
        grid.m_bChangeOfVar = true;
        
        for(nIdxTime = 2; nIdxTime < grid.m_nNbTimes; nIdxTime++)
        {
          ComputeRoots(grid.m_pdTimes[nIdxTime], nNbRoots, pRoots, bMinus);

          pdRootsChangeVar[nIdxTime] = pRoots[nIdxRealRoot].m_dRoot;
        }
      }
      else // bGridNeedChangeOfVar
      {    
        ComputeRoots(grid.m_pdTimes[0], nNbRoots, pRoots);

        for(nI = 0; nI < nNbRoots; nI++)
          if(pRoots[nI].m_iRootType != RootType_Virtual)
          {
            grid.m_dFixedBoundary = pRoots[nI].m_dRoot;
            break;
          }      
      } // IsGridNeedChangeOfVar()

    } // isCall && isConversion

  } // for iIdxGrid

  
# ifndef NDEBUG
  {
    m_cbParams.m_bStateFunctionShouldBeCalledWithTimeArgument = false;
  }
# endif

}

// maximum number of special space points for a sub-space mesh
const size_t NB_SEPCIAL_SPACE_POINTS = 10;

void CBMeshManager::GetSpecialSpacePointsForGrid
                    (
                      GridWithCV& grid,                           
                      double* pdEventLogSpot,
                      size_t& nNbEventSpot
                    )
{
  ASSERT_MSG(nNbEventSpot <= NB_SEPCIAL_SPACE_POINTS,
             "input number of special space points overflows.");

  // For now, this function is implemented very simply:
  // we only consider fixed boundary as special points
  //
  // idealy, we may need that the special points of grid[i+1]
  // exist also in the space mesh of grid[i] at its end time   

  if(nNbEventSpot > 0)
  {
    double dCorrection = 1;
    if(grid.m_pdCorrectionS)
      dCorrection = grid.m_pdCorrectionS[0] 
                  / grid.m_pdCorrectionS[grid.m_nNbTimes - 1];
    for(size_t nIdx = 0; nIdx < nNbEventSpot; nIdx++)
      pdEventLogSpot[nIdx] = log(pdEventLogSpot[nIdx] * dCorrection);

  }

  if(grid.m_dFixedBoundary > 0)
  {
    pdEventLogSpot[nNbEventSpot++] = log(grid.m_dFixedBoundary);
  }

  // Add conversion trigger(if any) as Spot event for initial time index
  // Note that this still can't guarantee that the trigger will be in the 
  // space mesh at each time step, but it still helps to have more points 
  // around the trigger, especially when there is a jump around the
  // trigger (dividend for example can introduce a jump)

  // note : conversion index is set in ComputeRootsAndUpdateGrids by the call
  // to HasContinousCallAt. This is weired, but... 
  size_t nIdxConversion = m_cbParams.GetIndexConversion();
  if (nIdxConversion != INVALIDINDEX)
  {
    double dTriggerRate = m_cbParams.GetConversions()
                        ->GetTriggerRate(nIdxConversion);

    double dTime = grid.m_pdTimes[0];
    
    finance::TriggerAsPercentageOf 
      trigger = m_cbParams.GetConversions()->GetTriggerAsPercentageOf();

    double dTriggerLevel = dTriggerRate * m_cbParams.GetConversions()
                         ->GetConversionPrice(dTime, trigger);

    if (dTriggerLevel > 0)
      pdEventLogSpot[nNbEventSpot++] = log(dTriggerLevel);
  }

  //call part
  size_t nIdxCall = m_cbParams.GetIndexCall();
  if ( nIdxCall != INVALIDINDEX && nIdxConversion != INVALIDINDEX)
  {
    double dTriggerRate = m_cbParams.GetCalls()->GetTriggerRate(nIdxCall);

    double dTime =grid.m_pdTimes[0];

    finance::TriggerAsPercentageOf
      trigger = m_cbParams.GetCalls()->GetTriggerAsPercentageOf();

    double dTriggerLevel = dTriggerRate * m_cbParams.GetConversions()
                         ->GetConversionPrice(dTime, trigger);

    if (dTriggerLevel > 0 )
      pdEventLogSpot[nNbEventSpot++] = log(dTriggerLevel);

  }

}

void CBMeshManager::MakeSpaceMeshForGrid
                    (
                      GridWithCV& grid, 
                      double dDelta,
                      double* pdEventLogSpot,
                      size_t nNbEventSpot,
                      double dLogSpotInitial
                    )
{
  int iTmp;
  size_t nI;

  grid.m_bSpaceGridByRef = false;
          
  iTmp = GenMeshSize
         ( 
           1, 
           dDelta / m_cbParams.GetMeshParams()->GetSpaceCompression(), 
           dDelta, 
           m_cbParams.GetMeshParams()->GetSpaceStretch(), 
           pdEventLogSpot, nNbEventSpot
         );

  grid.m_pdLogS = new double [iTmp];

  GenMesh
  ( 
    1, 
    dDelta / m_cbParams.GetMeshParams()->GetSpaceCompression(), 
    dDelta, 
    m_cbParams.GetMeshParams()->GetSpaceStretch(), 
    pdEventLogSpot, nNbEventSpot,grid.m_pdLogS, 
    iTmp
  );

  grid.m_nNbS = iTmp;

  iTmp = grid.GetIndexLog(dLogSpotInitial);

  double 
    *pd; 
  int 
    iSize;

  WidenGrid(grid.m_pdLogS, 
            grid.m_nNbS, 
            1.5,
            2. * 85. / (m_cbParams.GetNumParams()->GetNbSpaceSteps() - 1),
            dLogSpotInitial - 10,
            2,
            2,
            dLogSpotInitial + 10,
            &pd,
            iSize,
            iTmp);

  delete [] grid.m_pdLogS;
  grid.m_pdLogS = pd;

  grid.m_nNbS = iSize;

  grid.m_pdS = new double [grid.m_nNbS];   
  for(nI = 0; nI < grid.m_nNbS; nI++)
    grid.m_pdS[nI] = exp(grid.m_pdLogS[nI]);    

  grid.m_nNbCurrentS = iSize;
}


void CBMeshManager::ComputeMeshHelperData
                    (
                      double& dLogSpotMin,
                      double& dLogSpotMax,
                      double& dLogSpotInitial,
                      double& dDelta
                    )
{
  double
    dVariation,
    dInitialSpot = m_cbParams.GetSpotSharePrice();

  dLogSpotInitial = log(dInitialSpot);
  //____________________________________________________________________
  //
  // determination of spotmin and spotmax
  // The computation domain contains, at least, one year diffusion
  // it contains also the spot (InitialSpot/100)
  //____________________________________________________________________
  
  // Get the diffusion size from the used numerical model
  //using the original unshifted volatility 
  double
    dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
        (m_cbParams.GetStoppingTime(), dInitialSpot);
  
  dVariation = m_cbParams.GetMeshParams()->GetGridSpan() * 0.5 
             * m_cbParams.GetDiffusionSize(dSquaredTotalVol);  
  if(dVariation < 1)
    dVariation = 1;

  dLogSpotMin = dLogSpotInitial - dVariation;
  dLogSpotMax = dLogSpotInitial + dVariation;
  
  dDelta = (dLogSpotMax - dLogSpotMin) / 
            (m_cbParams.GetNumParams()->GetNbSpaceSteps() - 1);
}


void CBMeshManager::AdaptativeGrid()
{ 
  size_t
    nI,
    nIdxTime,
    nNbEventSpotOld = 0,
    nNbEventSpot = 0;
  
  double
    dDelta,
    dLogSpotInitial,
    dLogSpotMin,
    dLogSpotMax;

  bool
    bEventsChanged = true;

  double
    pdEventLogSpotOld[NB_SEPCIAL_SPACE_POINTS],
    pdEventLogSpot[NB_SEPCIAL_SPACE_POINTS];
  
  ComputeRootsAndUpdateGrids();

  ComputeMeshHelperData(dLogSpotMin, dLogSpotMax, dLogSpotInitial, dDelta);

  for(size_t nIdxGrid = 0; nIdxGrid < m_nNbGrids; nIdxGrid++)
  {
    /*
      we do the following things here
      1. calculate the special spot points for this grid
      2. build space mesh if we have to, otherwise we just
         copy the mesh of the last grid.
      3. we complete the grid information such as m_nNbCurrentSpots,
         m_pdCorrectionSpot.
      */
    
    GridWithCV&
      grid = m_pGrids[nIdxGrid];

    double pdEventsTmp[NB_SEPCIAL_SPACE_POINTS];
    size_t nNbEventsTmp = 0;

    if(nIdxGrid + 1 < m_nNbGrids &&
       m_pGrids[nIdxGrid + 1].m_dFixedBoundary > 0 )
    {
      pdEventsTmp[nNbEventsTmp++] = m_pGrids[nIdxGrid + 1].m_dFixedBoundary;
    }

    GetSpecialSpacePointsForGrid(grid, pdEventsTmp, nNbEventsTmp);

    // remove the points outside [dLogSpotMin, dLogSMax]
    nNbEventSpot = 0;
    pdEventLogSpot[nNbEventSpot++] = dLogSpotMin;
    pdEventLogSpot[nNbEventSpot++] = dLogSpotInitial;
    pdEventLogSpot[nNbEventSpot++] = dLogSpotMax;

    for(nI = 0; nI < nNbEventsTmp; nI++)
      if(pdEventsTmp[nI] < dLogSpotMax && pdEventsTmp[nI] > dLogSpotMin)
        pdEventLogSpot[nNbEventSpot++] = pdEventsTmp[nI];

    double 
      *pdEnd;
    std::sort(pdEventLogSpot, pdEventLogSpot + nNbEventSpot);
    pdEnd = std::unique(pdEventLogSpot, pdEventLogSpot + nNbEventSpot, 
                        numeric::IsEqual);

    nNbEventSpot = pdEnd - pdEventLogSpot;

    // check if the events are different from the last grid's
    bEventsChanged = true;
    if(nNbEventSpot == nNbEventSpotOld)
    {
      for(nI = 0; nI < nNbEventSpot; nI++)
        if(fabs(pdEventLogSpot[nI] - pdEventLogSpotOld[nI]) > 1.e-12)
          break;

      if(nI == nNbEventSpot)
        bEventsChanged = false;
    }

    if(!bEventsChanged) // just copy from previous space mesh
    {
      grid.m_bSpaceGridByRef   = true;
      grid.m_pdS = m_pGrids[nIdxGrid - 1].m_pdS;
      grid.m_pdLogS = m_pGrids[nIdxGrid - 1].m_pdLogS;
      grid.m_nNbS = m_pGrids[nIdxGrid - 1].m_nNbS;           
    }
    else // make the space mesh 
    { 
      MakeSpaceMeshForGrid(
                           grid, 
                           dDelta,
                           pdEventLogSpot,
                           nNbEventSpot,
                           dLogSpotInitial
                           );

    }

    if(grid.m_dFixedBoundary > 0)
    {
      grid.m_nNbCurrentS = grid.GetIndexLog(log(grid.m_dFixedBoundary));
      
      if ( grid.m_nNbCurrentS == INVALIDINDEX )
        // the fixed boundary is not in the space mesh
        grid.m_nNbCurrentS = grid.m_nNbS;
      else
        grid.m_nNbCurrentS++;
    }

    // update in case of change of variable inside the grid
    if(grid.m_bChangeOfVar)
    {
      double *pdRootsChangeVar = grid.m_pdCorrectionS;
      
      for(nIdxTime = 1; nIdxTime < grid.m_nNbTimes; nIdxTime++)
        pdRootsChangeVar[nIdxTime] /= pdRootsChangeVar[0];
      
      pdRootsChangeVar[0] = 1;       
    }

    for(nI = 0; nI < nNbEventSpot; nI++)
      pdEventLogSpotOld[nI] = pdEventLogSpot[nI];
    nNbEventSpotOld = nNbEventSpot;

    if( m_cbParams.HasNewShare() )
    {
      // TODO: NewShare.
      
    }
  }
}

void CBMeshManager::SetupRates()
{
  BackwardMeshManager::SetupRates();

  if( m_cbParams.GetCBLike().IsCrossCurrency() )
  {
    // Allocate the storage for the final derivative rate array in the 
    //cross currency case
    m_pdDerivativeRates = Array<double>(m_nNbTimes);
   
    Array<double>
      pdDerivativeRatesTmp = Array<double>(m_nNbTimes); 
    
    /*
      $ e^{ - \int_t^T r(s) ds } $ is a trivial solution of the backward PDE,
      so does $ e^{ \int_{ t_{ref} }^t r(s) ds } $ with $t$ the time variable
      since $ e^{ \int_{ t_{ref} }^t r(s) ds } 
            = e^{ - \int_t^T r(s) ds } * e^{ \int_{ t_{ref} }^T r(s) ds $
    */  
    m_cbParams.GetCBLike().GetDerivativeCurve()->GetCompoundFactor
      (m_pdTimes.Get(), pdDerivativeRatesTmp.Get(), m_nNbTimes);
    
    double
      dInvDeltaT, dOldInvDeltaT = 0.,
      dTmp1, dTmp2;
    
    for(size_t nIdxTime = m_nNbTimes - 2; nIdxTime < m_nNbTimes - 1; nIdxTime--)
    {
      dInvDeltaT = 1.0 / ( m_pdTimes[nIdxTime + 1] - m_pdTimes[nIdxTime] );

      if( m_pSchemeTypes[nIdxTime] == SchemeType_ThreeLevel )
      {
        double dC1 = (dInvDeltaT + dOldInvDeltaT);
        double dC2 = dOldInvDeltaT / dInvDeltaT / ( m_pdTimes[nIdxTime + 2] - 
                    m_pdTimes[nIdxTime] );
        
        dTmp1 = pdDerivativeRatesTmp[nIdxTime + 1] / 
                pdDerivativeRatesTmp[nIdxTime];   
        dTmp2 = pdDerivativeRatesTmp[nIdxTime + 2] / 
                pdDerivativeRatesTmp[nIdxTime + 1]; 

        m_pdDerivativeRates[nIdxTime] = dTmp1 * (dC1 - dC2 * dTmp2) - dC1 + dC2;
      }
      else if( m_pSchemeTypes[nIdxTime] == SchemeType_Implicit )
      {
        dTmp1 = pdDerivativeRatesTmp[nIdxTime + 1] / 
                pdDerivativeRatesTmp[nIdxTime];   
        m_pdDerivativeRates[nIdxTime] = (dTmp1 - 1.) * dInvDeltaT;
      }
      else if( m_pSchemeTypes[nIdxTime] == SchemeType_CrankNicolson )
      {
        dTmp1 = pdDerivativeRatesTmp[nIdxTime + 1] / 
                pdDerivativeRatesTmp[nIdxTime];  
        dTmp2 = m_pdDerivativeRates[nIdxTime + 1];

        m_pdDerivativeRates[nIdxTime] = 
          (dTmp1 * (1. - 0.5 * dTmp2 / dInvDeltaT) - 1.) * 2.0 * dInvDeltaT;
      }

      dOldInvDeltaT = dInvDeltaT;
    }

    // Compute the foreign exchange rates
    m_pdFXRates = Array<double>(m_nNbTimes);

    if (m_cbParams.GetCBLike().IsFixedQuanto() )
    {
      double
        dFixedFXRate = m_cbParams.GetCBLike().GetConversions()
                                              ->GetFixedFXRate();
      for(size_t nIdxTime = 0; nIdxTime < m_nNbTimes; nIdxTime++)
        m_pdFXRates[nIdxTime] = dFixedFXRate;
    }
    else
    {    
      Array<double> pdRatesTmp = Array<double>(m_nNbTimes);   
      
      m_cbParams.GetYieldCurve()->GetCompoundFactor(m_pdTimes.Get(), 
        pdRatesTmp.Get(), m_nNbTimes);   
      
      m_cbParams.GetCBLike().GetDerivativeCurve()->GetCompoundFactor
        (m_pdTimes.Get(), pdDerivativeRatesTmp.Get(), m_nNbTimes);
      
      double dTmp;
      
      m_pdFXRates[0] = m_cbParams.GetCBLike().GetSpotFXRate();
      
      dTmp = m_pdFXRates[0] * pdRatesTmp[0] / pdDerivativeRatesTmp[0];
      for(size_t nIdxTime = 1; nIdxTime < m_nNbTimes; nIdxTime++)
        m_pdFXRates[nIdxTime] = dTmp * pdDerivativeRatesTmp[nIdxTime] /
                                pdRatesTmp[nIdxTime];
    }

  }
}


} //namespace pricing

} //namespace ito33

