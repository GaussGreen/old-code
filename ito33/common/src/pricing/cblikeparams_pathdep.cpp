/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cblikeparams_pathdep.cpp
// Purpose:     PathDep stuff related to convertible like pricing
// Created:     2005/06/29
// RCS-ID:      $Id: cblikeparams_pathdep.cpp,v 1.3 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"

#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/anytimeforeverevent.h"
#include "ito33/pricing/quarterlyevent.h"
#include "ito33/pricing/ndayscallevent.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/conversionprovisions.h"
#include "ito33/pricing/callprovisions.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/finance/bondlike/cocotype.h"

namespace ito33
{

namespace pricing
{

typedef shared_ptr<PathDepEvent> PathDepEventPtr;

  using namespace finance;
  using namespace numeric;


std::vector< AutoPtr<CBLikeParams> >
CBLikeParams::ConstructPaths(const std::vector<double>& grid) const
{
  // Determine what types of path dependency are present 
  bool bHasCoCo = HasPathDepCoCo();
  bool bHasCoCall = HasPathDepCall();
 
  const size_t nNbPaths = grid.size();

  // Construct the params for each path
  std::vector< AutoPtr<CBLikeParams> > cblikeParams(nNbPaths);

  size_t nIdPath;
  for ( nIdPath = 0 ; nIdPath < nNbPaths; nIdPath++ )
  {
    // Clone the params, since params cannot be shared
    cblikeParams[nIdPath] = AutoPtr<CBLikeParams>( Clone() );

    // Manually set the analysis times
    cblikeParams[nIdPath]->SetAnalysisTime(m_dAnalysisTime);
  }

  // Update the params depending on the contingent clause.
  if ( bHasCoCo && !bHasCoCall )
  {
    // Contingent conversion only

    // Path 0 only has continuously observed conversion windows. 
    // Paths 1 and 2 have all periods, but no triggers, since they
    // represent "convert as of check date" problems
    cblikeParams[0]->GetConversions()->DeactivateDiscretePeriods();
    cblikeParams[1]->GetConversions()->RemoveTriggers();
    cblikeParams[2]->GetConversions()->RemoveTriggers();
  }
  else if ( !bHasCoCo && bHasCoCall ) // Contingent calls only
  {
    // Only the top pde (above trigger for N days) can call
    //need also to deactive the notice period for path 1..N-2.
    //to avoid extra work
    for ( nIdPath = 0 ; nIdPath < nNbPaths - 1; nIdPath++) 
    {
      cblikeParams[nIdPath]->GetCalls()->DeactivateSoftCall();

      if ( nIdPath >= 1 )
        cblikeParams[nIdPath]->GetCalls()->DeactivateNoticePeriod();
    }
  }
  else
  {
    // Both contingent conversion and calls

    size_t nCallPaths = GetCalls()->GetTriggerPeriod() + 1;
    for ( nIdPath = 0 ; nIdPath < nNbPaths; nIdPath++)
    {
      // Handle the CoCo part
      // The first nCallPath paths only have continuously observed conversion 
      // windows. The remaining paths (quarterly and forever) have all 
      // periods,  but no triggers, since they represent "convert as of check 
      // date"  problems
      if ( nIdPath < nCallPaths)
        cblikeParams[nIdPath]->GetConversions()->DeactivateDiscretePeriods();
      else
        cblikeParams[nIdPath]->GetConversions()->RemoveTriggers();

      // Handle the contingent call window part
      // Only the top path in each call group can call 
      if ( (nIdPath + 1) % nCallPaths != 0 )
        cblikeParams[nIdPath]->GetCalls()->DeactivateSoftCall();

      if ( (nIdPath+1)% nCallPaths != 0 && nIdPath%nCallPaths != 0)
        cblikeParams[nIdPath]->GetCalls()->DeactivateNoticePeriod();

    } //end loop over path
  }

  return cblikeParams;
}

std::list< shared_ptr<PathDepEvent> >
CBLikeParams::ConstructPathDepEvents
(const std::vector< AutoPtr<CBLikeParams> >& paths) const
{
  ASSERT_MSG(IsPathDependent(),
             "No path dependent clause found in PricePathDep");
 
  std::list< shared_ptr<PathDepEvent> > pathDepEventList;
  const size_t nNbPaths = paths.size();

  // Determine what types of path dependency are present 
  bool bHasCoCo = HasPathDepCoCo();
  bool bHasCoCall = HasPathDepCall();

  CallProvisions* pCall = GetCalls();
  ConversionProvisions* pConversion = GetConversions();

  // Need to make sure call/conversion windows are within

  // Conversion provisions are needed by the events to determine the trigger
  // level (not rate, which is constant) at a given time. The top-most
  // path should have all conversion periods, and will be initialized,
  // so can be used for this purpose. The provisions in m_cbparams
  // won't be initialized, so should not be used.
  ConversionProvisions*
    pInitializedConversion = paths[nNbPaths - 1]->GetConversions();

  if ( bHasCoCall )
  {
    // Setup the contingent call events

    // Get the call periods, and related data
    size_t nNbCalls        = pCall->GetNbCalls();
    size_t nCallDays       = pCall->GetTriggerPeriod();
    size_t nTriggerHistory = pCall->GetTriggerHistory(); 

    // Scan the periods for a contingent call clause
    size_t nIdx;
    for (nIdx = 0; nIdx < nNbCalls; nIdx++)
    {      
      // Extract the period information
      double dStartTime      = pCall->GetStartTime(nIdx);
      double dEndTime        = pCall->GetEndTime(nIdx);
      double dTriggerRate    = pCall->GetTriggerRate(nIdx);
      
      // Check if the current call period is contingent
      if ( dTriggerRate > 0 ) 
      {
        // make sure it is within the pricing range. Adjust if
        // partially in the range
        if ( (   IsEqualOrAfter(dStartTime, m_dStoppingTime)
              || IsEqualOrBefore(dEndTime, m_dValuationTime) )
          && !AreTimesEqual(dStartTime, m_dStoppingTime) )
          continue;

        if ( IsAfter(dEndTime, m_dStoppingTime) )
          dEndTime = m_dStoppingTime;

        if ( IsBefore(dStartTime, m_dValuationTime) )
          dStartTime = m_dValuationTime;

        shared_ptr<PathDepEvent> event;

        // Note: the call dates have already been shifted
        // when cb calls was constructed
   
        // Create an event at each day within the period
        double dEventTime = dStartTime - ONEDAY;

        while ( fabs( dEventTime - dEndTime ) > TIMETOLERANCE )
        {    
          dEventTime = dEventTime + ONEDAY;

          shared_ptr<PathDepEvent> pEvent;

          pEvent = PathDepEventPtr( new NDaysCallEvent
                                        (dStartTime, dEndTime, dEventTime, 
                                         dTriggerRate, pInitializedConversion,
                                         nCallDays, nTriggerHistory) );
        
          pathDepEventList.push_back(pEvent);       
         
        } // end while loop creating events in current window

      } // end check triggerrate

    } // end for loop over call windows

  } // if contingent calls

  if ( bHasCoCo )
  {
    // Setup the contingent conversion events

    // Get the conversion periods.
    size_t nNbProvisions = pConversion->GetNbConversions();

    // if contingent calls are present, the events need to know
    size_t nCallDays = 0;
    if ( bHasCoCall)
      nCallDays = pCall->GetTriggerPeriod();

    // Scan the periods for a coco clause
    for (size_t nIdx = 0; nIdx < nNbProvisions; nIdx++)
    {      
      // Extract the period information
      double dStartTime   = pConversion->GetStartTime(nIdx);
      double dEndTime     = pConversion->GetEndTime(nIdx);
      CoCoType cocoType = pConversion->GetPeriodType(nIdx);
      double dTriggerRate = pConversion->GetTriggerRate(nIdx);

      // make sure it is within the pricing range. Adjust if
      // partially in the range
      if ( (   IsEqualOrAfter(dStartTime, m_dStoppingTime)
            || IsEqualOrBefore(dEndTime, m_dValuationTime) )
        && !AreTimesEqual(dStartTime, m_dStoppingTime) )
        continue;

      if ( IsAfter(dEndTime, m_dStoppingTime) )
        dEndTime = m_dStoppingTime;

      // If the quarterly observation is before the observation
      // time, then do not include an event
      if ( IsBefore(dStartTime, m_dValuationTime) )
      {
        if ( cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate
          || cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter)
          continue;
        else        
          dStartTime = m_dValuationTime;
      }

      shared_ptr<PathDepEvent> event;

      // Use the above info to construct the appropriate event
      if ( cocoType == CoCoType_CheckAnyTimeAndConvertAsOfCheckDate )
      {
        event = PathDepEventPtr( new AnyTimeForeverEvent
                                     (dStartTime, dEndTime, dStartTime, 
                                      dTriggerRate, pInitializedConversion,
                                      nCallDays) );

        pathDepEventList.push_back(event);

        event = PathDepEventPtr( new AnyTimeForeverEvent
                                     (dStartTime, dEndTime, dEndTime, 
                                      dTriggerRate, pInitializedConversion,
                                      nCallDays) );

        pathDepEventList.push_back(event);
      }
      else if ( cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate )
      {
        // only one event, since observation is discrete
        event = PathDepEventPtr( new AnyTimeForeverEvent
                                     (dStartTime, dStartTime, dStartTime, 
                                      dTriggerRate, pInitializedConversion,
                                      nCallDays) );

        pathDepEventList.push_back(event);
      }
      else if (cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter)
      {
        event = PathDepEventPtr( new QuarterlyEvent
                                     (dStartTime, false, 
                                      dTriggerRate, pInitializedConversion,
                                      nCallDays) );

        pathDepEventList.push_back(event);

        event = PathDepEventPtr( new QuarterlyEvent
                                     (dEndTime, true, 
                                      dTriggerRate, pInitializedConversion, 
                                      nCallDays) );

        pathDepEventList.push_back(event);
      }

    } //end for loop

  } // if coco

  // Sorting is required if contingent conversion and call events are mixed
  if ( bHasCoCo && bHasCoCall )
    pathDepEventList.sort(std::greater< shared_ptr<PathDepEvent> >());

  return pathDepEventList;
}

std::vector<double> 
CBLikeParams::ConstructPathDepGrid(Model& /* model */) const
{
  bool bHasCoCo = HasPathDepCoCo();
  bool bHasCoCall = HasPathDepCall();
  
  // Calculate the number of paths
  size_t nNbPaths = 1;
  size_t nTriggerPeriod = 0;

  if ( bHasCoCo )
    nNbPaths = 3;

  if ( bHasCoCall )
  {
    nTriggerPeriod = GetCalls()->GetTriggerPeriod();

    nNbPaths *= nTriggerPeriod + 1;
  }

  std::vector<double> pdGridY(nNbPaths);
  for ( size_t nIdPath = 0 ; nIdPath < nNbPaths; nIdPath++ )
  {
    // The state variable for the contingent clauses is not needed/used, 
    // so just initialize with a counter
    pdGridY[nIdPath] = nIdPath;
  }

  return pdGridY;
}

size_t CBLikeParams::GetPathToSave(const std::vector<double>& /* grid */) const
{
  return 0;
}

void CBLikeParams::InitPaths(PathDepStructure& pathDepStruct)
{
  // Only the base paths need to be active at the beginning
  // for contingent call problems
  if ( HasPathDepCall() )
  {
    size_t nCallPaths = GetCalls()->GetTriggerPeriod() + 1;
    for (size_t nIdPath = 0; nIdPath < m_pathDepEvents.size(); nIdPath++)
      if ( (nIdPath % nCallPaths) != 0 )
        pathDepStruct.m_pbIsActive[nIdPath] = false;
  }
}


} // namespace pricing

} // namespace ito33
