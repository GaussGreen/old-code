/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/pathdeppricer.h
// Purpose:     path dependent pricer class
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: pathdeppricer.h,v 1.12 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/pathdeppricer.h
    @brief Path dependent pricer class.

    Implementation of the pricer class for path dependent contracts.
*/

#ifndef _ITO33_PRICING_PATHDEPPRICER_H_ 
#define _ITO33_PRICING_PATHDEPPRICER_H_

#include "ito33/sharedptr.h"
#include "ito33/list.h"

namespace ito33
{

namespace pricing
{
  class PathDepStructure;
  class PathDepEvent;

/// Path dependent pricer class
class PathDepPricer
{
public:

  /// Constructor.  Nothing to do.
  PathDepPricer() { }
  
  /** 
      Prices the contract.
      
      Output is saved in the pathDepStruct.

      @param pathDepStruct the structure holding all pricing data for each path
      @param pathDepEvents the events to apply
      @param bSaveOutputAfterEvents whether or not to save data in the structure
                                    after the last event has been applied
   */
  void Price(PathDepStructure& pathDepStruct, 
             const std::list< shared_ptr<PathDepEvent> >& pathDepEvents,
             bool bSaveOutputAfterEvents = true);


protected:

  /**
      Advances the solution for all paths by one timestep.

      @param pathDepStruct contains all information about the paths
   */
  void AdvanceOneStepForAll(PathDepStructure& pathDepStruct);

  /**
      Advances all paths to the specified time.

      @param dStopTime is the time to stop timestepping at
      @param pathDepStruct contains all information about the paths
   */
  virtual void AdvanceToTime(double dStopTime, PathDepStructure& pathDepStruct);

  /**
      Advances the solution for the specified path by one timestep.

      @param nIdx is the path to solve for
      @param pathDepStruct contains all information about the paths
   */
  void AdvanceOneStep(size_t nIdx, PathDepStructure& pathDepStruct);

  /**
      Applies the continuous events at the current time.

      @param pathDepStruct contains all information about the paths
   */
  void ApplyContinuousEvents(PathDepStructure& pathDepStruct);

  /**
      Removes the continuous events at the specified time

      @param dEventTime is the time at which to remove events
   */
  void RemoveContinuousEvents(double dEventTime);

  /**
      Applies the event.

      @param dEventTime the event time (used for continuous events)
      @param pEvent the event to apply
      @param pathDepStruct contains all information about the paths
   */
  virtual void ApplyEvent(double dEventTime,
                          shared_ptr<PathDepEvent> pEvent,
                          PathDepStructure& pathDepStruct);

  /**
      Turn off all paths except the path to save.

      If past the last event, no point in solving paths that will not
      contribute anything to the path to save.

      @param pathDepStruct contains all information about the paths
   */
  void TurnOffPaths(PathDepStructure& pathDepStruct);


  /// List of continuous events
  std::list< shared_ptr<PathDepEvent> > m_continuousEvents;

  /// Whether or not output should be saved
  bool m_bSaveOutput;

private:

  NO_COPY_CLASS(PathDepPricer);

}; // class PathDepPricer


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_PATHDEPPRICER_H_
