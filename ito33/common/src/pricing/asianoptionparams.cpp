/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/asianoptionparams.cpp
// Author:      ITO 33 Canada
// Created:     April 6, 2005
// RCS-ID:      $Id: asianoptionparams.cpp,v 1.8 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/binarysearch.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/payoffoption.h"
#include "ito33/finance/payoffconstant.h"
#include "ito33/finance/dividends.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/asianoption.h"
#include "ito33/pricing/averagingevent.h"
#include "ito33/pricing/optionliketype.h"
#include "ito33/pricing/asianoptionparams.h"
#include "ito33/pricing/asianoptionmeshmanager.h"
#include "ito33/pricing/constructeventtimes.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::pricing;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::AsianOptionParams);
}

AutoPtr<AsianOptionParams> AsianOptionParams::Clone()
{
  // Copy the underlying contract. The new params class will manage the
  // memory
  AutoPtr<pricing::AsianOption> clonedAsianOption( new AsianOption(m_asianOption) );

  // Construct and setup the cloned cb params
  AutoPtr<pricing::AsianOptionParams> 
    clonedParams(new AsianOptionParams(clonedAsianOption) );

  clonedParams->SetNumParams(      m_pNumParams );
  clonedParams->SetMeshParams(     m_pMeshParams );
  clonedParams->SetYieldCurve(     GetYieldCurve() );
  clonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  clonedParams->SetForeignCurve(   GetForeignCurve() );
  clonedParams->SetDividends(      GetDividends() );
  clonedParams->SetValuationTime(  GetValuationTime() );
  clonedParams->SetSpotSharePrice( GetSpotSharePrice() );  

  return clonedParams;
} //AsianOptionParams::Clone()

void AsianOptionParams::ConstructPayoff()
{
  /*
    The correct payoff has been set at the pricing level
    so the strike should never be negative
    if fixed strike option then strike == strike :)
    if floating strike option then strike == averaging.
  */

  /*
    we support the m_dAverage = 0 only because that in 
    AsianOption::ConstructAveragingGrid() we calls Init() of the default 
    params which doesn't really play a role in engine run but whose 
    m_dAverage is initialized to 0.
  */
  ASSERT_MSG( m_dAverage >= 0,
              "Negative or zero Average for Asian option."
            );

  switch ( m_asianOption.GetOptionType() )
  {
  case pricing::AsianOption_FloatingStrikeCall:
    m_pPayoff = shared_ptr<Payoff>( new PayoffCall(m_dAverage) );
    break;

  case pricing::AsianOption_FixedStrikeCall:  
    ASSERT_MSG( m_dStrike >= 0,
                "Negative strike for Asian option.");
    m_pPayoff = shared_ptr<Payoff>
                ( new PayoffConstant
                      (m_dAverage < m_dStrike ? 0. : m_dAverage - m_dStrike) );
    break;
     
  case pricing::AsianOption_FloatingStrikePut:  
    m_pPayoff = shared_ptr<Payoff>( new PayoffPut(m_dAverage) );
    break;
  
  case pricing::AsianOption_FixedStrikePut:
    ASSERT_MSG( m_dStrike >= 0,
                "Negative strike for Asian option.");
    m_pPayoff = shared_ptr<Payoff>
                ( new PayoffConstant
                      (m_dAverage > m_dStrike ? 0. : m_dStrike - m_dAverage) );
    break;
     
  default: 
    FAIL("No such Asian option type.");   
  
  }//switch

}//ConstructPayoff()


void AsianOptionParams::ConstructPathDepGrid(std::vector<double>& pdGridY,
       Model& model, bool bHasSimilarityReduction)
{

  double dCurrentAverage  = GetAsianOption().GetCurrentAverage();

  if ( bHasSimilarityReduction )
  {
    //pdGridY.push_back(0.);
    pdGridY.push_back( dCurrentAverage );
    return;
  }

  //construct the grid using the asian option mesh manager
  AutoPtr<pricing::AsianOptionParams> params = Clone();
  
  // For floating the strike is not set
  // so we need to adjust it otherwise
  // the average mesh will be wrong
  if ( !m_asianOption.HasFixedStrike() )
  {
    // Current average has been set
    if ( m_asianOption.GetCurrentAverage() > 0 )
      params->SetStrike( m_asianOption.GetCurrentAverage() );
    else
      params->SetStrike( GetSpotSharePrice() );
  }

  pricing::AsianOptionMeshManager optionMeshes(*params, model);

  params->Init();

  optionMeshes.SetupMe();
  optionMeshes.SetInitialState();

  size_t nIdS;
  size_t nGridSize        = optionMeshes.GetNbS();
  const double* pdGrid    = optionMeshes.GetS();

  for ( nIdS = 0; nIdS < nGridSize; nIdS++)
  {
    //add current average if not in the grid
    if ( (pdGrid[nIdS] > dCurrentAverage) 
      && (pdGrid[nIdS - 1] < dCurrentAverage) && 
        !numeric::IsEqual(pdGrid[nIdS], dCurrentAverage) 
        && !numeric::IsEqual( pdGrid[nIdS-1], dCurrentAverage) )
      pdGridY.push_back( dCurrentAverage ); 

    pdGridY.push_back( pdGrid[nIdS] );
  }
}

bool AsianOptionParams::IsPathDependent() const
{
  // It is path dependent if at least one of the following is true
  // - has a cash or pseudocash dividend
  // - initial conversion ratio rule and more than one reset date
  // - stock dependent vol and hazard rate

  if (   GetDividends() && GetDividends()->HasCashBetween
         ( GetDateFrom( GetValuationTime() ),
           GetDateFrom( GetStoppingTime() ) ) ) 
  {
    return true;
  } 

  // Only for floating strike
  if ( GetAsianOption().GetOptionType() 
                   == pricing::AsianOption_FixedStrikeCall
      || 
      GetAsianOption().GetOptionType() 
                   == pricing::AsianOption_FixedStrikePut )
  {
    return true;
  }

  return false;
}

std::list< shared_ptr<pricing::PathDepEvent> > 
AsianOptionParams::ConstructEvents(bool bHasSimilarityReduction)
{
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;

  const pricing::AsianOption &asianOpt = GetAsianOption();

  double dStartSamplingTime = asianOpt.GetAverageStartTime();
  double dStopSamplingTime  = asianOpt.GetAverageEndTime();
  double dValuationTime     = GetValuationTime();
  size_t nNbSampling        = asianOpt.GetNbSamplingAverages();
  size_t nNbSamplesUsed     = asianOpt.GetNbSamplesUsed();

  std::vector<double> pdEventTimes;
  std::vector<size_t> pnSamplingNumbers;

  ConstructSamplingEventTimes(dValuationTime, dStartSamplingTime, dStopSamplingTime,
    nNbSampling, nNbSamplesUsed, pdEventTimes, pnSamplingNumbers);

  // First event
  bool bIsFirstEvent = true;
  shared_ptr<pricing::AveragingEvent> pEvent
      ( 
        new pricing::AveragingEvent( pdEventTimes[0], 
                                      pnSamplingNumbers[0], 
                                      bIsFirstEvent, 
                                      bHasSimilarityReduction )
        );

  pathDepEvents.push_back(pEvent);

  bIsFirstEvent = false;
  for (size_t nIdx = 1; nIdx < pdEventTimes.size() ; nIdx++)
  {
     shared_ptr<pricing::AveragingEvent> pEvent
            ( 
              new pricing::AveragingEvent( pdEventTimes[nIdx], 
                                           pnSamplingNumbers[nIdx], 
                                           bIsFirstEvent, 
                                           bHasSimilarityReduction )
              );

     pathDepEvents.push_back(pEvent);
  }

  return pathDepEvents;
}


std::vector< AutoPtr<pricing::OptionParams> > 
AsianOptionParams::ConstructParam(const std::vector<double>& pdGridY)
{

  size_t nNbPaths = pdGridY.size();
  size_t nIdPath;

  // Clone the params and create the paths
  std::vector< AutoPtr<pricing::OptionParams> > pAsianOptionParams(nNbPaths);

  for ( nIdPath = 0 ; nIdPath < nNbPaths ; nIdPath++)
  {
    AutoPtr<pricing::AsianOptionParams> pAsianTmp = Clone();
  
    // Fixed call or put
    if ( m_asianOption.HasFixedStrike() )
      pAsianTmp->SetStrike( pAsianTmp->GetAsianOption().GetStrike() );
    else // floating
       pAsianTmp->SetStrike( pdGridY[nIdPath] );

    pAsianTmp->SetAverage( pdGridY[nIdPath] );
  
    pricing::AsianOptionParams *asianTmp = pAsianTmp.release();
    AutoPtr<pricing::OptionParams> pp(asianTmp);
    pAsianOptionParams[nIdPath] = pp;

  } //end loop constructing paths

  return pAsianOptionParams;
}

size_t
AsianOptionParams::GetPathToSave(const std::vector<double>& pdGridY) const
{
  ASSERT_MSG( !pdGridY.empty(), "Asian params: error in averaging grid.");

  // similarity reduction is used
  if ( pdGridY.size() == 1 )
    return 0;

  // The path to save is the one containing the current average
  double dCurrentAverage = GetAsianOption().GetCurrentAverage();

  size_t nIdx = BinSearch(&pdGridY[0],pdGridY.size(),dCurrentAverage);

  if ( !numeric::IsEqual( pdGridY[nIdx], dCurrentAverage) )
   nIdx--;

  return nIdx;
}
