/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/event.h
// Purpose:     a generic event class
// Author:      David Pooley
// Created:     2003/08/14
// RCS-ID:      $Id: event.h,v 1.16 2004/11/12 17:02:19 zhang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/event.h
    @brief Declaration of the base event class.

    event defines the events that can happen during pricing (dividend,
    coupon, discrete barrier, etc.) 
*/

#ifndef _ITO33_PRICING_EVENT_H_
#define _ITO33_PRICING_EVENT_H_

#include "ito33/common.h"
#include "ito33/constants.h"

#include <cmath>

namespace ito33
{

namespace pricing
{

/** Event type. 

    Used for sorting the events, assuming a forward time direction. Each
    type represents an event category. For example, ET_Dividend can be
    a fixed dividend or a proportional dividend.
 */
enum EventType
{
  ET_Reset,
  ET_Dividend,
  ET_Payment,
  
  ET_Max
};


/** Base event class
    
    All events must be derived from  this base class.  The only difference
    between events is the implementation of the ApplyToPrice and 
    ApplyToGreek functions.
*/
class Event
{
public:
 
  /**
    ctor

    @param dTime time of the event
    */
  Event(double dTime) : m_dTime(dTime), m_IsAppliedAfterConstraints(false) { }

  /// dtor
  virtual ~Event() { }

  /**
    Function to apply the event to an array of price values 

    @param pdS the spot prices
    @param pdValues array of values as a function against spot. It is
        normally the price value, but it can also be the quantity with
        the same feature as the price (ex. the straight bond)
    @param nNbS the size of the array
    */
  virtual void ApplyToPrice(const double *pdS,
                            double* pdValues,
                            size_t nNbS) const = 0;

  /**
    Function to apply the event to Greek values pdValues on grid pdS

    @param pdS the spot prices
    @param pdValues array of Greek values
    @param nNbS the size of the array
    */
  virtual void ApplyToGreek(const double *pdS,
                            double* pdValues,
                            size_t nNbS) const = 0;

  /// Returns the time at which the event takes place
  double GetTime() const { return m_dTime; }

  /// Returns the type (category) of the event
  EventType GetType() const { return m_eventType; }

  /**
   Indicate if the event should be applied after constraints. By default,
   this is false.  Events such as resets must set this to true.

   @return false if the event is applied before constraints, true otherwise
 */
  bool IsAppliedAfterConstraints() const
  {
    return m_IsAppliedAfterConstraints;
  }

  /**
    less comparison operator
    
    @param event1 the first event
    @param event2 the second event to compare
    */
  friend bool operator < (const Event& event1, const Event& event2)
  {
    // If events occur at the same time, check event type for the order
    if ( fabs(event1.m_dTime - event2.m_dTime) < TIMETOLERANCE )
      return  event1.m_eventType < event2.m_eventType;
    else
      return event1.m_dTime < event2.m_dTime;
  }

  /**
    great comparison operator
    
    @param event1 the first event
    @param event2 the second event to compare
    */
  friend bool operator > (const Event& event1, const Event& event2)
  {
    // If events occur at the same time, check event type for the order
    if ( fabs(event1.m_dTime - event2.m_dTime) < TIMETOLERANCE )
      return  event1.m_eventType > event2.m_eventType;
    else
      return event1.m_dTime > event2.m_dTime;
  }

protected:

  /// time of the event
  double m_dTime;

  /// type of the event
  EventType m_eventType;

  /// bIsAfterConstraints
  bool m_IsAppliedAfterConstraints;

};


} // namespace pricing

} // ito33 namespace

#endif

