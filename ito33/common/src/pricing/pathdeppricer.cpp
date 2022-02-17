/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pathdep/pathdeppricer.h
// Purpose:     path dependent pricer class
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: pathdeppricer.cpp,v 1.14 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <list>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"
#include "ito33/pricing/steppertridiag.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/params.h"
#include "ito33/pricing/pathdeppricer.h"

namespace ito33
{

using numeric::AreTimesEqual;

namespace pricing
{


void PathDepPricer::Price(
     pricing::PathDepStructure& pathDepStruct, 
     const std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
     bool bSaveOutputAfterEvents)
{ 

  // Reset in case this function is called more than once  
  m_continuousEvents.clear();

  // Read events in reverse
  std::list< shared_ptr<PathDepEvent> >::const_reverse_iterator iterEvents;
  iterEvents = pathDepEvents.rbegin();

  // Valuation (stop) time should be the same in all paths
  double dStopTime = pathDepStruct.m_path[0].params->GetValuationTime();

  // Only save data after the path dep events are done, and only if
  // instructed to save output after the last event
  m_bSaveOutput = false;

  // Loop over the events. Be careful with the final time, which is not
  // currently recorded as a path dependent event. Multiple events at the
  // same time are supported, but are assumed to be already properly sorted.
  double dEventTime = -1.0;
  while (dEventTime != dStopTime)
  {
    if ( iterEvents != pathDepEvents.rend() )
    {
      dEventTime = (*iterEvents)->GetTime();
    }
    else
    {
      dEventTime = dStopTime;

      // Only save output if instructed
      if ( bSaveOutputAfterEvents )
      {
        m_bSaveOutput = true;
        pathDepStruct.InitPathToSave();
      }

      // Turn off paths except the path to save if more pricing is needed
      TurnOffPaths(pathDepStruct);
    }

    // Advance to the event time.  Multiple events at the same time are all
    // applied below.
    AdvanceToTime( dEventTime, pathDepStruct );
      
    // Apply the current events and advance
    while ( iterEvents != pathDepEvents.rend() 
          && AreTimesEqual(dEventTime, (*iterEvents)->GetTime()) )
    {
      ApplyEvent(dEventTime, *iterEvents, pathDepStruct);
      iterEvents++;
    } 

  } // while timestepping

  // If output was never saved, must be path dep from valuation to maturity.
  // Only save the final values.
  if ( m_bSaveOutput == false)
  {
    m_bSaveOutput = true;
    pathDepStruct.InitPathToSave();
  }

  pathDepStruct.Finalize();
}


void PathDepPricer::AdvanceToTime(double dStopTime, PathDepStructure& pathDepStruct)
{
  // Price until the stop time, which is typically an event time.
  // If there are no continuous events, the solution of all paths
  // can simply be advanced to the stop time. Otherwise, need to advance
  // all paths step by step while applying the continuous events
  size_t nNbPaths = pathDepStruct.m_path.size();
  if ( m_continuousEvents.size() == 0 )
  {
    // Should be more efficient to solve each path for as long as
    // possible, instead of proceeding in lock-step
    for (size_t nIdx = 0; nIdx < nNbPaths; nIdx++)
    {
      while 
      (
      !AreTimesEqual(pathDepStruct.m_path[nIdx].meshes->GetTime(), dStopTime)
      && pathDepStruct.m_path[nIdx].meshes->TryGoAhead() 
      )
      {
        AdvanceOneStep(nIdx, pathDepStruct);

        if ( m_bSaveOutput == true )
          pathDepStruct.UpdateOutput(nIdx);
      }
    } // for loop over paths
  }
  else
  {
    // Go step by step
    double dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();
    while ( !AreTimesEqual(dCurrentTime, dStopTime) )
    {
      AdvanceOneStepForAll(pathDepStruct);
      ApplyContinuousEvents(pathDepStruct);
      dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();

      ASSERT_MSG( m_bSaveOutput == false, 
                  "Cannot save output while continuous events are active");

    } // while timestepping step by step

  } // if continuous events or not
}


void PathDepPricer::AdvanceOneStepForAll(PathDepStructure& pathDepStruct)
{

  // Make one step in each path. Assume that the time meshes are
  // synchronized
  size_t nNbPaths = pathDepStruct.m_path.size();
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)  
  {
    if ( pathDepStruct.m_path[nIdx].meshes->TryGoAhead() )
    {
      AdvanceOneStep(nIdx, pathDepStruct);
    }
  }

  // Some debugging checks
  #ifndef NDEBUG
    double dTimeCheck = pathDepStruct.m_path[0].meshes->GetTime();
    for (nIdx = 1; nIdx < nNbPaths; nIdx++)
    {
      double dTime = pathDepStruct.m_path[nIdx].meshes->GetTime(); 
      ASSERT_MSG( AreTimesEqual(dTimeCheck, dTime),
        "Time mismatch in PathDepPricer::AdvanceOneStep");
    }
  #endif
}


void PathDepPricer::AdvanceOneStep(size_t nIdx, PathDepStructure& pathDepStruct)
{
  // This code is essentially the non-path-dependent engine.

  // This step is done outside the if loop to be sure
  // that everything stays in sync. In particular events
  // are still advanced appropriately.
  pathDepStruct.m_path[nIdx].instdata->UpdateBeforeStep();
  
  // Only do the solve if the path is active. Still need the other lines
  // to advance the time, update params, etc
  if (pathDepStruct.m_pbIsActive[nIdx] == true) 
  {
    pathDepStruct.m_path[nIdx].stepper->Run();

    if ( pathDepStruct.m_path[nIdx].params->HasEventsNow()
         && m_bSaveOutput == true )
      pathDepStruct.UpdateOutput(nIdx);
  }

  pathDepStruct.m_path[nIdx].instdata->DoEvents();
}


void PathDepPricer::ApplyContinuousEvents(PathDepStructure& pathDepStruct)
{
  std::list< shared_ptr<pricing::PathDepEvent> >::iterator iterContinuousEvents;
  iterContinuousEvents = m_continuousEvents.begin();

  while ( iterContinuousEvents != m_continuousEvents.end() )
  {
    (*iterContinuousEvents)->Apply(pathDepStruct);
    ++iterContinuousEvents;
  }
}


void PathDepPricer::RemoveContinuousEvents(double dEventTime)
{
  // Note only one item is removed!
  std::list< shared_ptr<pricing::PathDepEvent> >::iterator iterContinuousEvents;
  iterContinuousEvents = m_continuousEvents.begin();
  while ( iterContinuousEvents != m_continuousEvents.end() )
  {
    if ( AreTimesEqual( (*iterContinuousEvents)->GetStartTime(), dEventTime) )
    {
      m_continuousEvents.erase(iterContinuousEvents);
      break;
    }
    ++iterContinuousEvents;
  }
}


void PathDepPricer::ApplyEvent(double dEventTime,
                               shared_ptr<PathDepEvent> pEvent,
                               PathDepStructure& pathDepStruct)
{
  if ( pEvent->IsContinuous() == false ||
       AreTimesEqual( pEvent->GetEndTime(), pEvent->GetStartTime() ) )
  {
    pEvent->Apply(pathDepStruct);
  }
  else
  {
    // continuous event. Check if at start or end or window
    if ( AreTimesEqual( pEvent->GetEndTime(), dEventTime ) )
    {
      m_continuousEvents.push_back( pEvent );
      pEvent->Apply(pathDepStruct);
    }
    else if ( AreTimesEqual( pEvent->GetStartTime(), dEventTime ) )
    {
      RemoveContinuousEvents(dEventTime);   
    }
  }
}


void PathDepPricer::TurnOffPaths(PathDepStructure& pathDepStruct)
{
  // Turn off all the paths
  size_t nPath = pathDepStruct.m_path.size();    
  for ( size_t nIdPath = 0; nIdPath < nPath; nIdPath++ )  
   pathDepStruct.m_pbIsActive[nIdPath] = false;

  // Enable the path to save
  size_t nPathToSave = pathDepStruct.GetPathToSave();
  pathDepStruct.m_pbIsActive[nPathToSave] = true;
}

} // namespace pricing

} // namespace ito33
