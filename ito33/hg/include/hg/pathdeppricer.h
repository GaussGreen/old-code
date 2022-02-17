/////////////////////////////////////////////////////////////////////////////
// Name:        hg/pathdeppricer.h
// Purpose:     path dependent pricer class for HG
// Created:     2006/03/02
// RCS-ID:      $Id: pathdeppricer.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/pathdeppricer.h
    @brief Path dependent pricer class for HG.

    Implementation of the HG pricer class for path dependent contracts.
*/

#ifndef _HG_PATHDEPPRICER_H_ 
#define _HG_PATHDEPPRICER_H_

#include "ito33/sharedptr.h"
#include "ito33/list.h"

#include "ito33/pricing/pathdeppricer.h"

#include "ito33/numeric/predicatetime.h"

#include "hg/pathdepstructure.h"


namespace ito33
{

namespace pricing
{
  class PathDepEvent;
}

namespace hg
{
 
/// Path dependent pricer class for HG
template<class InstDataType, 
         class MeshManagerType, 
         class ParamsType,
         class NumOutputType, 
         class StepperType>
class PathDepPricer : pricing::PathDepPricer
{
public:

  /// Constructor.  Nothing to do.
  PathDepPricer(PathDepStructure<InstDataType, 
                                 MeshManagerType, 
                                 ParamsType, 
                                 NumOutputType, 
                                 StepperType>& pathDepStruct) 
    : pricing::PathDepPricer(),
      m_pathDepStruct(pathDepStruct)
  { }
  
  /** 
      Prices the contract.
      
      Output is saved in the pathDepStruct.

      @param pathDepEvents the list of path dependent events needed for pricing
      @param bSaveOutputAfterEvents whether or not to save data in the structure
                                    after the last event has been applied
   */
  void 
  Price(const std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
        bool bSaveOutputAfterEvents = true)
  {
    pricing::PathDepPricer::Price(m_pathDepStruct, pathDepEvents, 
                                  bSaveOutputAfterEvents);
  }
 
protected:

  /**
      Applies the event.

      Virtual function in base class.  This HG version applies the 
      event to each regime.

      @param dEventTime the event time (used for continuous events)
      @param pEvent the event to apply
      @param pathDepStruct contains all information about the paths
   */
  void ApplyEvent(double dEventTime,
                  shared_ptr<pricing::PathDepEvent> pEvent,
                  pricing::PathDepStructure& pathDepStruct)
  {
    if ( pEvent->IsContinuous() == false ||
      numeric::AreTimesEqual( pEvent->GetEndTime(), pEvent->GetStartTime() ) )
    {
      size_t nNbRegimes = m_pathDepStruct.GetNbRegimes();
      for (size_t nIdxRegime = 0; nIdxRegime < nNbRegimes; nIdxRegime++)
      {
        m_pathDepStruct.SetCurrentRegime(nIdxRegime);
        pEvent->Apply(pathDepStruct);
      }
    }
    else
    {
      // continuous event. Check if at start or end or window
      if ( numeric::AreTimesEqual( pEvent->GetEndTime(), dEventTime ) )
      {
        m_continuousEvents.push_back( pEvent );

        size_t nNbRegimes = m_pathDepStruct.GetNbRegimes();
        for (size_t nIdxRegime = 0; nIdxRegime < nNbRegimes; nIdxRegime++)
        {
          m_pathDepStruct.SetCurrentRegime(nIdxRegime);
          pEvent->Apply(pathDepStruct);
        }
      }
      else if ( numeric::AreTimesEqual( pEvent->GetStartTime(), dEventTime ) )
      {
        RemoveContinuousEvents(dEventTime);   
      }
    } // if continuous event

  }

  /**
      Advance all paths to the specified time.

      Virtual in base class.  This version checks for the end of subgrids.

      @param dStopTime is the time to stop timestepping at
      @param pathDepStruct contains all information about the paths
  */
  void AdvanceToTime(double dStopTime, 
                     pricing::PathDepStructure& pathDepStruct)
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
        !numeric::AreTimesEqual(pathDepStruct.m_path[nIdx].meshes->GetTime(), 
                                dStopTime)
        && pathDepStruct.m_path[nIdx].meshes->TryGoAhead() 
        )
        {
          // Check if at end of subgrid (typically only relevant for cb-like
          // problems)
          if( pathDepStruct.m_path[nIdx].meshes->IsEndOfGrid())
          {
            pathDepStruct.m_path[nIdx].instdata->UpdateAtEndOfGrid();
            
            if (m_bSaveOutput)
              pathDepStruct.UpdateOutputEndOfGrid(nIdx);            
          }
          else
          {
            AdvanceOneStep(nIdx, pathDepStruct);

            if ( m_bSaveOutput )
              pathDepStruct.UpdateOutput(nIdx);
          }
        }
      } // for loop over paths
    }
    else
    {
      // Go step by step
      double dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();
      while ( !numeric::AreTimesEqual(dCurrentTime, dStopTime) )
      {
        AdvanceOneStepForAll(pathDepStruct);
        ApplyContinuousEvents(pathDepStruct);
        dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();

        ASSERT_MSG( m_bSaveOutput == false, 
                    "Cannot save output while continuous events are active");

      } // while timestepping step by step

    } // if continuous events or not

  }

private:

  /// The path dep structure (HG type needed needed when applying events)
  PathDepStructure<InstDataType, 
                   MeshManagerType, 
                   ParamsType, 
                   NumOutputType, 
                   StepperType>& m_pathDepStruct;

  NO_COPY_CLASS(PathDepPricer);

}; // class PathDepPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_PATHDEPPRICER_H_
