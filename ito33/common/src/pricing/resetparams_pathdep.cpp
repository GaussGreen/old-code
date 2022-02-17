/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/resetparams_pathdep.cpp
// Purpose:     PathDep stuff related to resettable pricing
// Created:     2005/06/30
// RCS-ID:      $Id: resetparams_pathdep.cpp,v 1.9 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/binarysearch.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/resetflooredby.h"

#include "ito33/finance/dividends.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/resetevent.h"
#include "ito33/pricing/resetparams.h"

namespace ito33
{

namespace pricing
{
  using namespace finance;


bool ResetParams::IsPathDependent() const
{
  if (!HasActiveResetDate() )
    return false;

  // It is path dependent if at least one of the following is true
  // - has a cash or pseudocash dividend
  // - initial conversion ratio rule and more than one reset date
  if ( m_pDividends 
       && m_pDividends->HasCashBetween
          (GetDateFrom(m_dValuationTime), GetDateFrom(m_dStoppingTime))) 
  {
    return true;
  } 

  // if conversion by initial ratio only one
  // reset date is allowed
  if (   m_reset.GetFlooredBy() == ResetFlooredBy_InitialConversionPrice 
      && m_reset.GetResetTimes().size() > 1 )
    return true;

  return false;
}

std::vector<double> ResetParams::ConstructPathDepGrid(Model& model) const
{
  std::vector<double> pdGridY;

  // Basic idea: Use a cbmeshmanager to create a grid. Limit the
  // grid by the min and max conversion prices as determined by the floor
  // and cap values. Since conversion prices are reset around parity
  // (i.e. around the spot price at time t), it makes sense to center
  // the conversion grid around the starting spot price
  
  double dCurrentConvPrice = m_reset.GetCurrentConversionPrice();
  double dInitialConvPrice = m_reset.GetInitialConversionPrice();

  // Get the mesh manager.  Work on clone to be safe 
  AutoPtr<CBLikeParams> cblikeparams( Clone() );
  CBMeshManager cbmeshes(*cblikeparams, model);

  cblikeparams->Init();

  cbmeshes.SetupMe();
  cbmeshes.SetInitialState();

  size_t nGridSize = cbmeshes.GetNbS();
  const double* pdGrid = cbmeshes.GetS();
    
  // find the max and min values for the conversion ratio by assuming
  // the cap/floor are hit each reset date
  bool bIsResetFlooredByPrevailing = false; 
  if ( m_reset.GetFlooredBy() == ResetFlooredBy_PrevailingConversionPrice )
    bIsResetFlooredByPrevailing = true;

  const std::vector<double>& pdResetTimes = m_reset.GetResetTimes();
  const std::vector<double>& pdResetCaps = m_reset.GetCapRates();
  const std::vector<double>& pdResetFloors = m_reset.GetFloorRates();
  const std::vector<double>& pdMultiplier = m_reset.GetMultipliers();

  size_t nNbResetDates = pdResetTimes.size();
  double dMinFloor = dCurrentConvPrice;
  double dMaxCap = dCurrentConvPrice;
  
  if ( !bIsResetFlooredByPrevailing )
    dMinFloor = dInitialConvPrice;

  std::vector<double> pdSpecialPoints;

  pdSpecialPoints.push_back(dCurrentConvPrice);

  for (size_t nIdx = 0; nIdx < nNbResetDates; nIdx++)
  {   
    // Add these boundary values to the mesh. Duplicates are removed below.
    dMaxCap *= pdResetCaps[nIdx];
    pdSpecialPoints.push_back(dMaxCap);

    if ( bIsResetFlooredByPrevailing )
    {
      dMinFloor *= pdResetFloors[nIdx];
      pdSpecialPoints.push_back(dMinFloor);
    }
    else
    {
      // floor is calculated wrt initial ratio. Does not chain
      // at each reset date as for prevailing
      double dTmpFloor = dInitialConvPrice * pdResetFloors[nIdx];
      pdSpecialPoints.push_back(dTmpFloor);
      if ( dTmpFloor < dMinFloor )
        dMinFloor = dTmpFloor;
    }

  } // loop over reset dates

  // sort the vector of special points
  std::sort( pdSpecialPoints.begin(), pdSpecialPoints.end() );

  // Shift the grid by the initial multiplier. This works well for
  // the first reset period. 
  // TODO: shift grids in the reset events for later periods?
  size_t nIdxSpot = BinSearch(pdGrid, nGridSize, m_dSpot - 1.e-8);  
  
  double dMult = m_dSpot * pdMultiplier[0] / pdGrid[nIdxSpot];
  //dMult = dSpot / pdGrid[nIdxSpot - 1];
  //dMult = dCurrentConvPrice / pdGrid[nIdxSpot - 1];

  // make the Y grid of conversion prices. Limit the range due to caps/floors
  pdGridY.clear();
  for (size_t nIdS = 0; nIdS < nGridSize; nIdS++)
  {
    double dConvPrice = pdGrid[nIdS] * dMult;

    if (    numeric::IsEqualOrGreater(dConvPrice, dMinFloor) 
         && numeric::IsEqualOrLess(dConvPrice, dMaxCap) )
      pdGridY.push_back(dConvPrice);
  }

  // create temporary vector
  std::vector<double> pdGridNewY( pdGridY.size() + pdSpecialPoints.size() );

  // merge y grid and special point container
  std::merge(pdGridY.begin(), pdGridY.end(),
             pdSpecialPoints.begin(), pdSpecialPoints.end(),
             pdGridNewY.begin());

  // remove to keep only unique elements
  std::vector<double>::iterator 
    newEnd1 = std::unique(pdGridNewY.begin(), pdGridNewY.end(),
                          numeric::IsEqual);

  pdGridNewY.erase(newEnd1, pdGridNewY.end());

  pdGridY.clear();
  pdGridY = pdGridNewY;

  return pdGridY;
} 

size_t ResetParams::GetPathToSave(const std::vector<double>& grid) const
{
  double dCurrentConvPrice = m_reset.GetCurrentConversionPrice();
  
  const size_t nNbPaths = grid.size();
  size_t nIdPath;

  // find where the initial conversion price ended up
  for (nIdPath = 0; nIdPath < nNbPaths; nIdPath++)
    if ( numeric::IsEqual(grid[nIdPath], dCurrentConvPrice) )
      break;

  ASSERT_MSG(nIdPath < nNbPaths, "Cannot find the path to save");

  return nIdPath;
}

std::vector< AutoPtr<CBLikeParams> >
ResetParams::ConstructPaths(const std::vector<double>& grid) const
{
  const size_t nNbPaths = grid.size();

  // Clone the params and create the paths
  std::vector< AutoPtr<CBLikeParams> > pCBLikeParams(nNbPaths);
  
  for (size_t nIdPath = 0 ; nIdPath < nNbPaths; nIdPath++)
  {
    // Need access to pricing reset, so keep correct type for now
    AutoPtr<ResetParams> pClone( Clone() );

    // Set the new ratio in both conversions and the reset
    double dNewRatio = m_reset.GetNominal() / grid[nIdPath];
    pClone->GetReset().SetCurrentConversionPrice( grid[nIdPath] );
    pClone->GetConversions()->SetRatios(dNewRatio);

    // store as cblikeparams
    pCBLikeParams[nIdPath] = AutoPtr<CBLikeParams>( pClone.release() );

  } //end loop constructing paths

  return pCBLikeParams;
}

// Note: for a cross currency, we use the conversion price converted 
// into the currency of derivative by using the fixed FX rate as Y grid.
// this requires no change to ResetEvent and a few minor changes here
std::list< shared_ptr<PathDepEvent> >
ResetParams::ConstructPathDepEvents
(const std::vector< AutoPtr<CBLikeParams> >& /* paths */) const
{
  // Create the list of events
  std::list< shared_ptr<PathDepEvent> > pathDepEvents;
    
  double dCurrentConvPrice = m_reset.GetCurrentConversionPrice();
  double dInitialConvPrice = m_reset.GetInitialConversionPrice();

  const std::vector<double>& pdCapRates = m_reset.GetCapRates();
  const std::vector<double>& pdFloorRates = m_reset.GetFloorRates();
  const std::vector<double>& pdMultipliers = m_reset.GetMultipliers();
  const std::vector<double>& pdResetTimes = m_reset.GetResetTimes();

  bool bIsResetFlooredByPrevious = false; 
      
  if ( m_reset.GetFlooredBy() == ResetFlooredBy_PrevailingConversionPrice )
    bIsResetFlooredByPrevious = true;

  // need to track the last reset event so paths can be turned off
  bool bIsLastReset = true;

  for (size_t nIdx = 0; nIdx < pdResetTimes.size(); nIdx++)
  {
    double dEventTime = pdResetTimes[nIdx];

    // do not include reset dates on or before valuation time
    // or on or after maturity
    if ( numeric::IsEqualOrBefore(dEventTime, m_dValuationTime) || 
         numeric::IsEqualOrAfter(dEventTime, m_dStoppingTime) )
      continue;

    shared_ptr<PathDepEvent> pEvent; 

    // the coefficient that the spot will be multiplied by
    double 
      dMultiplier = pdMultipliers[nIdx] * GetFXRate(dEventTime);

    pEvent = make_ptr( new ResetEvent(dEventTime,
                                      pdCapRates[nIdx],
                                      pdFloorRates[nIdx],
                                      dMultiplier,
                                      dInitialConvPrice,
                                      dCurrentConvPrice,
                                      bIsResetFlooredByPrevious,
                                      bIsLastReset) );
 
    pathDepEvents.push_back(pEvent);

    // dates are ordered forward in time, so when pricing backward the 
    // first event is the last reset event
    bIsLastReset = false;
  }

  return pathDepEvents;
}

void ResetParams::InitPaths(PathDepStructure& /* pathDepStruct */)
{
}


} // namespace pricing

} // namespace ito33
