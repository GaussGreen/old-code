/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/varianceswapparams.cpp
// Purpose:     variance swap params class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapparams.cpp,v 1.33 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/predicatedouble.h"
#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/mesh/optionspacemesh.h"
#include "ito33/numeric/gammadistribution.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/payoffconstant.h"
#include "ito33/finance/returntype.h"
#include "ito33/finance/swappayofftype.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapevent.h"
#include "ito33/pricing/varianceswapevent_special.h"
#include "ito33/pricing/constructeventtimes.h"
#include "ito33/pricing/conditionalswapevent.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::VarianceSwapParams);

namespace pricing
{

typedef shared_ptr<VarianceSwapEvent> VarSwapEventPtr;

AutoPtr<VarianceSwapParams> VarianceSwapParams::Clone()
{
  // Copy the underlying contract. The new params class will manage the
  // memory
  AutoPtr<VarianceSwap> clonedVarianceSwap( new VarianceSwap(m_varianceSwap) );

  // Construct and setup the cloned params
  AutoPtr<VarianceSwapParams> 
    clonedParams(new VarianceSwapParams(clonedVarianceSwap) );

  clonedParams->SetNumParams(      m_pNumParams );
  clonedParams->SetMeshParams(     m_pMeshParams );
  clonedParams->SetYieldCurve(     GetYieldCurve() );
  clonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  clonedParams->SetForeignCurve(   GetForeignCurve() );
  clonedParams->SetDividends(      GetDividends() );
  clonedParams->SetValuationTime(  GetValuationTime() );
  clonedParams->SetSpotSharePrice( GetSpotSharePrice() );  

  return clonedParams;

} // Clone()


void VarianceSwapParams::Init()
{
  Params::Init();

  ConstructDividendEvents();

  ConstructPayoff();
}


void VarianceSwapParams::ConstructPayoff()
{
  double dValue = 0.0;

  // Check if valuing the fixed leg of conditional swap
  if ( m_dConditionalPayoff >= 0.0)
    dValue = ComputeConditionalPayoff();
  else
  {
    dValue = m_varianceSwap.GetPayoffValue( m_dAvgSqrReturns );

    if ( m_varianceSwap.GetSwapPayoffType() == finance::SwapPayoff_Call )
      dValue = std::max( dValue, 0.0 );
    else if ( m_varianceSwap.GetSwapPayoffType() == finance::SwapPayoff_Put )
      dValue = std::max( -dValue, 0.0 );
  }

   m_pPayoff = make_ptr( new finance::PayoffConstant(dValue) );

} // ConstructPayoff()


bool VarianceSwapParams::HasSimilarityReduction()
{
  // Check for cash or pseudocash dividends
  if (   GetDividends() 
      && GetDividends()->HasCashBetween(GetDateFrom(GetValuationTime()),
                                        GetDateFrom(GetStoppingTime()) ) ) 
    return false;

  // Gamma swaps do not admit a similarity reduction due to the S_i / S_0 term
  if ( m_varianceSwap.GetSwapPayoffType() == finance::SwapPayoff_Gamma )
    return false;

  if ( m_varianceSwap.GetSwapPayoffType() != finance::SwapPayoff_Standard )
    return false;

  // Cannot use a similarity reduction if the variance swap has corridors 
  // (the corridor will either be missed entirely or always used if there 
  // is only one previous path).  This also covers Conditional payoffs.
  if (   m_varianceSwap.GetUpCorridorBarrier() > 0.0 
      || m_varianceSwap.GetDownCorridorBarrier() > 0.0 )
    return false;

  return true;
}


std::vector<double> 
VarianceSwapParams::ConstructMasterGrid(const Model& model)
{
  // the grid to be built and returned
  std::vector<double> pdMasterGrid( GenerateSpaceMesh(model) );
  
  // Recenter around the previous share price
  double dPreviousSharePrice = m_varianceSwap.GetPreviousSharePrice();
  double dScale = dPreviousSharePrice / GetSpotSharePrice();

  size_t nNbPoints = pdMasterGrid.size();
  size_t nIdxMax = 0;
  for (size_t nIdx = 0; nIdx < nNbPoints; nIdx++)
  {
    pdMasterGrid[nIdx] = dScale * exp(pdMasterGrid[nIdx]);

    // Limit the grid to 50 times the previous share price. This
    // still represents a huge return.
    if ( pdMasterGrid[nIdx] < dPreviousSharePrice * 50)
      nIdxMax++;
  }

  std::vector<double> pdReturnGrid(nIdxMax);
  for (size_t nIdx = 0; nIdx < nIdxMax; nIdx++)
    pdReturnGrid[nIdx] = pdMasterGrid[nIdx];  

  return pdReturnGrid;
}


std::vector<double> 
VarianceSwapParams::ConstructAvgSqrReturnGrid()
{
  // In theory, we need a Z grid that extends to the maximum possible
  // squared return.  Otherwise, an event will compute a new Z value
  // greater than the max grid value, and the algorithm must extrapolate.
  // In practice, there are several things to consider:
  // 1) for variance swaps, everything is linear in Z, so linear
  //    interpolation and extrapolation is exact
  // 2) the maximum possible squared return on a large S grid will be
  //    prohibitively large.  The grid will be capped.
  // 3) the largest log return is much smaller than the largest 
  //    actual return.  Different caps should be used.

  double dAvgSqrReturn = m_varianceSwap.GetCurrentAvgSqrReturn();
  std::vector<double> pdGrid;

  // indicate if we can use only two points in the Z-grid direction
  if ( IsSpecialVarianceSwap() )
  {
    if ( dAvgSqrReturn > 0.0 )
    {
      pdGrid.push_back( dAvgSqrReturn );
      pdGrid.push_back( 10. * dAvgSqrReturn );
    }
    else
    {
      // need this since this is the path that will be saved
      pdGrid.push_back(0.0); 
      pdGrid.push_back(1.e-3);
    }

    return pdGrid;
  }

  // Determine the max Z value.  Assume the largest reasonable jump
  // in price is a factor of 50.  Adjust for caps.
  double dR = 50.0;
  if (m_varianceSwap.GetReturnType() == finance::Return_Log)
    dR = log(50.0);
  double dZMax = dR * dR;

  double dCap = m_varianceSwap.GetCapMultiplier();
  if (dCap > 0.0)
  {
    // Find when the cap on the volatility strike is the same as the
    // pricing Z state variable
    double dVolStrike = m_varianceSwap.GetVolatilityStrike();
    size_t nFreq = m_varianceSwap.GetAnnualReturnFrequency();

    double dCapMax = dCap * dCap * dVolStrike * dVolStrike / double(nFreq);

    // The cap is on the final realized volatility (or variance), not
    // on the individual samples.  The realized volatility can thus
    // be above the cap during the life of the contract.  The price
    // will still be close to zero above the cap though, so only
    // extend to 10 times past this level, and never past the return
    // level computed above (ie. in case the cap is huge).  The magic
    // number 10 was determined by experiment. In theory, it should
    // depend on the number of samples remaining, since with more
    // samples, the more likely it is that the final realized volatility
    // will be under the cap.
    dZMax = std::min(10.0*dCapMax, dZMax);
  }

  // Generate the grid in three stages.  We want to use the chi^2 
  // distribution, want to extend to dZMax, and need to include the
  // current average return.  However, the chi^2 functions of degree 
  // one have problems beyond about 30, so generate a mesh to 
  // dZSmallMax and then extend to dZMax.  After, scale to force the
  // current average.
  double dZSmallMax = std::min(5e-4, dZMax);
  double dChi2CdfMax = numeric::Chi2Cdf(1, dZSmallMax);

  // Don't need as many points as the 'S' mesh, but do enforce a minimum
  size_t nNbPoints = GetNumParams()->GetNbSpaceSteps() / 4;
  nNbPoints = nNbPoints < 34 ? 34 : nNbPoints;

  pdGrid.push_back(0.);

  // Generate the mesh 
  size_t nIdZ;
  for (nIdZ = 1; nIdZ < nNbPoints; nIdZ++)
  {    
    double dFraction = double(nIdZ)/(nNbPoints-1);
    double dZ = numeric::InvChi2Cdf( 1, dFraction * dChi2CdfMax );

    pdGrid.push_back( dZ );
  }

  // Extend to dZMax. The last point in the generated grid above
  // is not always reliable.
  nNbPoints = pdGrid.size();
  double dFact = pdGrid[nNbPoints - 2] / pdGrid[nNbPoints - 3];
  dFact = std::max(dFact, 1.01);
  dFact = std::min(dFact, 2.0);

  double dZ = pdGrid[nNbPoints-2] * dFact;
  pdGrid[nNbPoints - 1] = dZ;
  while (dZ < dZMax)
  {
    dFact *= 1.02;
    if (dFact > 2.0)
      dFact = 2.0;

    dZ *= dFact;
    pdGrid.push_back(dZ);    
  }

  // Scale grid to make sure the current return is in the grid
  nNbPoints = pdGrid.size();
  if (dAvgSqrReturn > DOUBLETOLERANCE)
  {
    if (dAvgSqrReturn < pdGrid[1])
    {
      pdGrid[1] = dAvgSqrReturn;
    }
    else if (dAvgSqrReturn >= pdGrid[nNbPoints - 1])
    {
      pdGrid[nNbPoints-1] = dAvgSqrReturn;
    }
    else
    {
      // find scale factor. Scale up or down, whichever gives less change.
      double dScale = 1.0;
      size_t nIdx = 1;
      while (pdGrid[nIdx] < dAvgSqrReturn)
        nIdx++;

      if ( (dAvgSqrReturn - pdGrid[nIdx-1]) > (pdGrid[nIdx] - dAvgSqrReturn) )
        dScale = dAvgSqrReturn / pdGrid[nIdx];        
      else
        dScale = dAvgSqrReturn / pdGrid[nIdx-1];
        
      // actually do the scaling
      for (nIdZ = 1; nIdZ < nNbPoints; nIdZ++)
        pdGrid[nIdZ] *= dScale;
    }

  } // if current return is not zero

  // Add the max point into the grid.  Usually not needed, since the max
  // point is somewhat arbitrary, but can be useful if the scaling for 
  // the current average decreased the last point
  double dLast = pdGrid[pdGrid.size()-1];
  if ( dLast < dZMax && !numeric::IsEqual(dLast, dAvgSqrReturn) )
    pdGrid[pdGrid.size()-1] = dZMax;

  return pdGrid;
}


std::vector<double> 
VarianceSwapParams::ConstructPreviousSpotGrid
                    (const Model& model, bool bIsSimilarityReduction)
{
  // the grid to be built and returned
  std::vector<double> pdPreviousSpotGrid;

  // A similarity reduction can be used
  if ( bIsSimilarityReduction )
  {
    double dPreviousSharePrice = m_varianceSwap.GetPreviousSharePrice();
    pdPreviousSpotGrid.push_back(dPreviousSharePrice);
    
    // Sometimes to avoid extrapolation, an extra path is solved.  Then
    // interpolation between the two paths can be done instead of a 
    // similarity extrapolation.  This is not needed for variance swaps,
    // but this note is left as a reminder.
    // double dPreviousSpot = dSpot / 10.0;    
    // pdPreviousSpotGrid.push_back(dPreviousSpot);

    return pdPreviousSpotGrid;
  }
  
  // If no similarity reduction, make the 'P' grid the same as the master
  // 'S' grid
  std::vector<double> pdGrid = ConstructMasterGrid(model);

  // Only use roughly half the points
  double dPreviousSharePrice = m_varianceSwap.GetPreviousSharePrice();
  size_t nIdxPrevious = 0;
  size_t nNbPrevious = pdGrid.size();
  for (size_t nIdx = 0; nIdx < nNbPrevious; nIdx++)
  {
    if ( numeric::IsEqual(dPreviousSharePrice, pdGrid[nIdx]) )
    {
      nIdxPrevious = nIdx;
      break;
    }
  }
  ASSERT_MSG(nIdxPrevious != nNbPrevious, "path to save not found.");

  size_t nStart = nIdxPrevious % 2;

  size_t nNbNew = nNbPrevious / 2;
  if (nNbPrevious % 2 == 1 && nStart == 0)
    nNbNew += 1;

  std::vector<double> pdReturnGrid(nNbNew);
  for (size_t nIdx = nStart, nIdxNew = 0; 
       nIdx < nNbPrevious; 
       nIdx += 2, nIdxNew++)
  {
    pdReturnGrid[nIdxNew] = pdGrid[nIdx];
  }

  return pdReturnGrid;
}


std::vector<double> VarianceSwapParams::ConstructConditionalGrid()
{
  // Construct grid from current conditional count to the number
  // of samples.  Since events always look to higher counts, paths
  // from 0 to current count are not needed.
  size_t nNbSamples = m_varianceSwap.GetNbSamplingReturns();
  int iCurrentCount = m_varianceSwap.GetCurrentConditionalCount();
  size_t nSize = nNbSamples - iCurrentCount;

  std::vector<double> pdYGrid(nSize);
  
  for (size_t nIdx = 0; nIdx < nSize; nIdx++)
    pdYGrid[nIdx] = nIdx + iCurrentCount;

  return pdYGrid;
}


std::vector< AutoPtr<pricing::VarianceSwapParams> > 
VarianceSwapParams::ConstructParams(
  const std::vector<double>& pdAvgSqrReturnGrid, 
  const std::vector<double>& pdPreviousSpotGrid)
{
  size_t nNbReturns = pdAvgSqrReturnGrid.size();
  size_t nNbPrevious = pdPreviousSpotGrid.size();

  size_t nNbPaths = nNbReturns * nNbPrevious;

  // Clone the params and create the paths by looping over the return
  // grid (outer loop) and previous spot grid (inner loop)
  std::vector< AutoPtr<VarianceSwapParams> > pParams(nNbPaths);

  size_t nIdxPath = 0;
  for (size_t nIdxReturn = 0; nIdxReturn < nNbReturns; nIdxReturn++)
  {
    for (size_t nIdxPrevious = 0; nIdxPrevious < nNbPrevious; nIdxPrevious++)
    {
      AutoPtr<VarianceSwapParams> pVarSwapTmp = Clone();
  
      pVarSwapTmp->SetAvgSqrReturn( pdAvgSqrReturnGrid[nIdxReturn] );
      pVarSwapTmp->SetPreviousSpot( pdPreviousSpotGrid[nIdxPrevious] );

      pParams[nIdxPath] = pVarSwapTmp;

      nIdxPath++;

    } // loop over previous grid

  } // loop over return grid

  return pParams;
}


std::vector< AutoPtr<pricing::VarianceSwapParams> > 
VarianceSwapParams::ConstructConditionalParams(
  const std::vector<double>& pdYGrid)
{
  // Clone the params and create the paths by looping over the Y 
  // (conditional count) grid
  size_t nNbPaths = pdYGrid.size();
  std::vector< AutoPtr<VarianceSwapParams> > pParams(nNbPaths);

  size_t nIdxPath;
  for (nIdxPath = 0; nIdxPath < nNbPaths; nIdxPath++)
  {
    AutoPtr<VarianceSwapParams> pVarSwapTmp = Clone();
  
    pVarSwapTmp->SetConditionalPayoff( pdYGrid[nIdxPath] );

    pParams[nIdxPath] = pVarSwapTmp;

  } // loop over paths

  return pParams;
}


size_t VarianceSwapParams::GetPathToSave(
  const std::vector<double>& pdAvgSqrReturnGrid, 
  const std::vector<double>& pdPreviousSpotGrid)
{
  // Return is the outer variable, previous is the inner.
  // At the valuation date, we want the current average of the squared
  // returns (either zero or set by user if valuation date is past the
  // contract start date), and the current spot share price
  double dAvgSqrReturn = m_varianceSwap.GetCurrentAvgSqrReturn();

  double dPSharePrice = m_varianceSwap.GetPreviousSharePrice();

  size_t nNbPrevious = pdPreviousSpotGrid.size();
  size_t nNbReturn = pdAvgSqrReturnGrid.size();

  // find in (inner) previous grid
  size_t nIdxPrevious = nNbPrevious;
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbPrevious; nIdx++)
  {
    if ( numeric::IsEqual(dPSharePrice, pdPreviousSpotGrid[nIdx]) )
    {
      nIdxPrevious = nIdx;
      break;
    }
  }
  ASSERT_MSG(nIdxPrevious != nNbPrevious, "path to save not found.");

  // find in (outer) return grid
  size_t nIdxReturn = nNbReturn;
  for (nIdx = 0; nIdx < nNbReturn; nIdx++)
  {
    if ( numeric::IsEqual(dAvgSqrReturn, pdAvgSqrReturnGrid[nIdx]) )
    {
      nIdxReturn = nIdx;
      break;
    }
  }
  ASSERT_MSG(nIdxReturn != nNbReturn, "path to save not found.");

  return (nIdxReturn * nNbPrevious) + nIdxPrevious;
}


size_t VarianceSwapParams::GetConditionalPathToSave(
  const std::vector<double>& pdYGrid)
{

  // Get and then find the current count
  double dCurrentCount = (double) m_varianceSwap.GetCurrentConditionalCount();
  
  size_t nNbPaths = pdYGrid.size();
  size_t nIdx = 0;
  size_t nPathToSave = 0;
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    if ( numeric::IsEqual(pdYGrid[nIdx], dCurrentCount) )
    {
      nPathToSave = nIdx;
      break;
    }
  }
  ASSERT_MSG(nIdx != nNbPaths, "conditional path to save not found.");

  return nPathToSave;
}


std::list< shared_ptr<pricing::PathDepEvent> > 
VarianceSwapParams::ConstructEvents(bool bIsSimilarityReduction)
{
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;

  double dValuationTime     = GetValuationTime();
  double dStartSamplingTime = m_varianceSwap.GetStartTimeOfSamplingPeriod();
  double dStopSamplingTime  = m_varianceSwap.GetMaturityTime();
  size_t nNbSampling        = m_varianceSwap.GetNbSamplingReturns();
  size_t nNbSamplesUsed     = m_varianceSwap.GetNbSamplesUsed();

  std::vector<double> pdEventTimes;
  std::vector<size_t> pnSamplingNumbers;

  ConstructSamplingEventTimes(dValuationTime, dStartSamplingTime, 
    dStopSamplingTime, nNbSampling, nNbSamplesUsed, 
    pdEventTimes, pnSamplingNumbers);

  // Get the return type
  finance::ReturnType returnType = m_varianceSwap.GetReturnType();

  // If the valuation date is before the sampling start time add a fake event
  // with sampling number 0.  This tells the event to copy the solution along
  // the diagonal S=P for Z = CurrentAvg
  shared_ptr<VarianceSwapEvent> pEvent;
  bool bIsFirstEvent = true;
  
  bool bIsSpecial = IsSpecialVarianceSwap();

  if ( dValuationTime < dStartSamplingTime )
  {
    if ( bIsSpecial )
      pEvent = VarSwapEventPtr
               ( new VarianceSwapEventSpecial
                     ( dStartSamplingTime, 0,
                       returnType, bIsFirstEvent,
                       bIsSimilarityReduction,
                       m_varianceSwap ) );
    else
      pEvent = VarSwapEventPtr( new VarianceSwapEvent
                                    ( dStartSamplingTime, 0,
                                      returnType, bIsFirstEvent,
                                      bIsSimilarityReduction,
                                      m_varianceSwap ) );
   
    pathDepEvents.push_back(pEvent);

    bIsFirstEvent = false;
  }
  
  // Get the first real event directly  
  if ( bIsSpecial ) 
    pEvent = VarSwapEventPtr( new VarianceSwapEventSpecial
                                  ( pdEventTimes[0], pnSamplingNumbers[0],
                                    returnType, bIsFirstEvent,
                                    bIsSimilarityReduction,
                                    m_varianceSwap ) );
  else
    pEvent = VarSwapEventPtr( new VarianceSwapEvent
                                  ( pdEventTimes[0], pnSamplingNumbers[0],
                                    returnType, bIsFirstEvent , 
                                    bIsSimilarityReduction,
                                    m_varianceSwap ) );
  pathDepEvents.push_back(pEvent);

  bIsFirstEvent = false;

  for (size_t nIdx = 1; nIdx < pdEventTimes.size(); nIdx++)
  {
    shared_ptr<VarianceSwapEvent> pEvent;

    if ( bIsSpecial )
      pEvent = VarSwapEventPtr( new VarianceSwapEventSpecial
                                    ( pdEventTimes[nIdx], 
                                      pnSamplingNumbers[nIdx],
                                      returnType, bIsFirstEvent,
                                      bIsSimilarityReduction,
                                      m_varianceSwap ) );
    else  
      pEvent = VarSwapEventPtr( new VarianceSwapEvent
                                    ( pdEventTimes[nIdx], 
                                      pnSamplingNumbers[nIdx],
                                      returnType, bIsFirstEvent, 
                                      bIsSimilarityReduction,
                                      m_varianceSwap ) );

    pathDepEvents.push_back(pEvent);
  }

  return pathDepEvents;
}


std::list< shared_ptr<pricing::PathDepEvent> > 
VarianceSwapParams::ConstructConditionalEvents()
{
  // Use the same event times as for the floating leg of the swap
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;

  double dValuationTime     = GetValuationTime();
  double dStartSamplingTime = m_varianceSwap.GetStartTimeOfSamplingPeriod();
  double dStopSamplingTime  = m_varianceSwap.GetMaturityTime();
  size_t nNbSampling        = m_varianceSwap.GetNbSamplingReturns();
  size_t nNbSamplesUsed     = m_varianceSwap.GetNbSamplesUsed();

  std::vector<double> pdEventTimes;
  std::vector<size_t> pnSamplingNumbers;

  ConstructSamplingEventTimes(dValuationTime, dStartSamplingTime, 
    dStopSamplingTime, nNbSampling, nNbSamplesUsed, 
    pdEventTimes, pnSamplingNumbers);
  
  if ( pdEventTimes.empty() )
    return pathDepEvents;

  // Add the first event manually
  shared_ptr<ConditionalSwapEvent>
    pEvent( new ConditionalSwapEvent(pdEventTimes[0], true, m_varianceSwap) );

  pathDepEvents.push_back(pEvent);

  // The remaining events
  for (size_t nIdx = 1; nIdx < pdEventTimes.size(); nIdx++)
  {
    pEvent = make_ptr( new ConditionalSwapEvent
                           ( pdEventTimes[nIdx], false, m_varianceSwap ) );

    pathDepEvents.push_back(pEvent);
  }

  return pathDepEvents;

}


double VarianceSwapParams::ComputeConditionalPayoff()
{
  // See Philippe's paper for the payoff
  size_t nNbSamples = m_varianceSwap.GetNbSamplingReturns();
  double dStrike = m_varianceSwap.GetVolatilityStrike();
  double dValue = m_dConditionalPayoff / nNbSamples * dStrike;
  if ( m_varianceSwap.GetSwapType() == finance::Swap_Variance )
    dValue *= dStrike; 

  return dValue;
}


bool VarianceSwapParams::IsSpecialVarianceSwap() const
{
  // A standard variance swap allows for some pricing simplifications and/or
  // analytic formulas. For example, the solution is linear in terms of the 
  // average of the squared returns, so the PDE approach only needs two paths 
  // in the Z direction.
  if (    m_varianceSwap.IsVanilla()
       && m_varianceSwap.GetSwapPayoffType() == finance::SwapPayoff_Standard )
    return true;

  return false;
}

bool VarianceSwapParams::IsAnalytical() const
{
  finance::SwapPayoffType payoffType = m_varianceSwap.GetSwapPayoffType();

  // Some special variance swaps can be solved analytically
  if (    m_varianceSwap.IsVanilla()
       && (    payoffType == finance::SwapPayoff_Standard 
            || payoffType == finance::SwapPayoff_Gamma) )
    return true;

  return false;
}

bool VarianceSwapParams::IsForwardStarting() const
{
  double dValuationTime     = GetValuationTime();
  double dStartSamplingTime = m_varianceSwap.GetStartTimeOfSamplingPeriod();

  if ( numeric::IsEqualOrAfter(dValuationTime, dStartSamplingTime) )
    return false;
  
  return true;
}

} // namespace pricing

} // namespace ito33
