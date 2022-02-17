/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/pathdepevent.h
// Purpose:     a generic path dependent event class
// Author:      Yann and David 
// Created:     8/08/2004
// RCS-ID:      $Id: pathdepevent.h,v 1.16 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/pathdepevent.h
    @brief Declaration of the base path dependent event class.
 */

#ifndef _ITO33_PRICING_PATHDEPEVENT_H_
#define _ITO33_PRICING_PATHDEPEVENT_H_

#include "ito33/beforestd.h"
#include <cmath>
#include <functional>
#include "ito33/afterstd.h"

#include "ito33/constants.h"
#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/backwardmeshmanager_withspacemesh.h"

namespace ito33
{

namespace pricing
{
 
  class PathDepStructure;

/** 
    Base path dependent event class.
    
    All path dependent events must be derived from this base class.      
 */
class PathDepEvent
{
public:
 
  /**
      ctor

      @param dStartTime start time of the event (continuous events)
      @param dEndTime end time of the event
      @param dTime time of the event
   */
  PathDepEvent(double dStartTime, double dEndTime, double dTime)
             : m_dStartTime(dStartTime), m_dEndTime(dEndTime), m_dTime(dTime)
  { }

  /// dtor
  virtual ~PathDepEvent() { }

  /**
      Function to apply the event to a path dependent structure

      @param pathDepStruct the path dependent structure containing prices, grids
   */
  void Apply(PathDepStructure& pathDepStruct) 
  {
    double dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();
    
    if ( fabs(m_dStartTime - m_dEndTime) < TIMETOLERANCE)
    {
      // If discrete event, then apply only
      ApplyAtTime(pathDepStruct);
    }
    else if ( fabs(dCurrentTime - m_dStartTime ) < TIMETOLERANCE)
    {
      // If continuous event and at start time, must apply first, then
      // do whatever is needed for event clean-up (assumes pricing
      // is backwards in time)
      ApplyAtTime(pathDepStruct);
      ApplyAtStartTime(pathDepStruct);
    }
    else if ( fabs( dCurrentTime - m_dEndTime ) < TIMETOLERANCE )
    {
      // If continuous event and at end time, must initialize the event
      // first, then apply (assumes pricing is backwards in time)
      ApplyAtEndTime(pathDepStruct);
      ApplyAtTime(pathDepStruct);
    }
    else
    {
      // If continuous event and not at start or end times, just need
      // to apply the event
      ApplyAtTime(pathDepStruct);
    } 
  }
 
  /**
      Function which returns the time at which the event takes place.

      @return double
   */
  double GetTime() const 
  { 
    return m_dTime; 
  }

  /**
      Function which returns the start time of the event.

      @return double
   */
  double GetStartTime() const
  {
    return m_dStartTime;
  }

  /**
      Function which returns the end time of the event.

      @return double
   */
  double GetEndTime() const
  {
    return m_dEndTime;
  }

  /**
      Indicates whether or not the event is done continuously
      i.e at each timestep or discretely.
   */
  virtual bool IsContinuous() const
  {
    //discrete events
     if ( fabs( m_dStartTime - m_dEndTime ) < TIMETOLERANCE )
       return false;

     return true;
  }

  /**
      less comparison operator.
    
      @param event1 the first event
      @param event2 the second event to compare
   */
  friend bool operator < (const PathDepEvent& event1, const PathDepEvent& event2)
  {
    return event1.m_dTime < event2.m_dTime;
  }

  /**
      great comparison operator.
    
      @param event1 the first event
      @param event2 the second event to compare
   */
  friend bool operator > (const PathDepEvent& event1, const PathDepEvent& event2)
  {
    return event1.m_dTime > event2.m_dTime;
  }

protected:

  /// time of the event
  double m_dTime;

  /// start of window
  double m_dStartTime;

  ///end of window
  double m_dEndTime;

  /**
      Function to apply event at the start time.

      @param pathDepStruct the path dependent structure containing prices, grids
   */
  virtual void ApplyAtStartTime(PathDepStructure& pathDepStruct) const = 0;


  /**
      Function to apply event End Time.

      @param pathDepStruct the path dependent structure containing prices, grids
   */
  virtual void ApplyAtEndTime(PathDepStructure& pathDepStruct) const =  0;

  /**
      Function to apply event at time m_dTime.

      @param pathDepStruct the path dependent structure containing prices, grids
   */
  virtual void ApplyAtTime(PathDepStructure& pathDepStruct) const = 0;

};

} // namespace pricing

} // ito33 namespace

namespace std
{

/*
  A specialization of greater for shared pointers of PathDepEvent,
  used for ascending sort of path dependent events. 
*/
template<>
inline bool greater< ito33::shared_ptr<ito33::pricing::PathDepEvent> >::operator()
     (const ito33::shared_ptr<ito33::pricing::PathDepEvent>& pEvent1, 
      const ito33::shared_ptr<ito33::pricing::PathDepEvent>& pEvent2) const
{
  return *pEvent1 < *pEvent2;
}

} // namespace std

#endif  // _ITO33_PRICING_PATHDEPEVENT_H_
