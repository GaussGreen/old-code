/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/baseeventmanager.h
// Purpose:     Base template class to manage events
// Author:      Nabil
// Created:     2004/02/18
// RCS-ID:      $Id: baseeventmanager.h,v 1.25 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/baseeventmanager.h
    @brief Template class to manage events

    Template class to manage events (such as : dividends, coupons, etc) 
    (and also, for the convertible bond: puts provisions, 
     monodates calls provisions and monodates conversions provisions)
    for all contract types

    @todo the iterator inside of the eventmanager may be moved outside so that
          an event manager can be used for backward and forward problems at
          the same time. And the current time can be maintained then by Params.
          A temporary solution is to add more checks on m_dTime.

    @todo Adding same event at same time should probably be forbided.
 */

#ifndef _ITO33_PRICING_BASEEVENTMANAGER_H_
#define _ITO33_PRICING_BASEEVENTMANAGER_H_

#include "ito33/list.h"
#include "ito33/debug.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/problemtype.h"

namespace ito33
{

namespace pricing
{


/**
    Template class to manage events.

    Template class to manage events (such as : dividends, barriers, etc) 
    (and also, for the convertible bond: puts provisions, 
     monodates calls provisions and monodates conversions provisions)
    for all contract types.
    
    When this template class is instantiated for a given event type, the
    comparator for the shared pointer of this event type needs to be
    implemented in the header file, see eventmanager.h for an example.

    @param T is an event class.
 */
template <class T>
class BaseEventManager
{
public:
  
  /// Container for the list of events
  typedef std::list< shared_ptr<T> > Events;

  /** 
      Constructor, sets the problem type to backward.
   */
  BaseEventManager();

  /**
      @param problemType if forward, then traverse the list of events forward
             in time. If backward, traverse the list backwards
 
      @todo need probably to reinitialize the internal state variable
   */
  void SetProblemType(ProblemType problemType)
  {
    m_problemType = problemType;
  }

  /// Default dtor is ok

  /**
      Adds an event to the event list.

      Note that this function will not check if there has two dividend events
      at a same date, the verification is done by the class Dividends.

      the same is true for other kinds of events
    
      @param pEvent shared pointer to the event to be added
   */
  void AddEvent(const shared_ptr<T>& pEvent);

  /**
      Sets the initial state.

      @param dTime A time used to help determine the initial state.
   */
  void SetInitialState(double dTime);

  /** 
      Sets the current time of the PDE solve

      @param dTime The current time
   */
  void SetCurrentTime(double dTime);

  /** 
      Checks if any events occur at the current time
 
      Returns true if at least one event occurs at the current time,
      false otherwise.  Assumes that the list of events is sorted,
      and that external calls are ordered.
   */
  bool HasEventsNow() const;

  /** 
      Gets eventually an event occuring at the current time

      Returns an event occuring at the current time (if any), or NULL if 
      no more events occur at the current time.  Assumes that the list of 
      events is sorted, and that external calls are ordered.
   */
  const T* GetEvent() const;

  /** 
      Gets the event list. 

      Even if we can use the GetEvent() function to traverse the events, it's 
      a bit strange when event manager is used not for applying the events:
     
      First, it can only traverse in the way defined by the problem type
      (forward or backward)

      Second, SetInitialState needs to be called before and after (may not 
      in this function).
      
      @return the reference to the event list
   */
  const Events& GetAll() const 
  { 
    Sort();

    return m_events; 
  }

  /**
      Initialization. Clean up so that an event manager can be reused.
   */
  void Init() 
  { 
    m_events.clear();

    // By default, set current time to a non valid value
    m_dTime = -1.;

    // By default, no sorting has occured
    m_bIsNotSorted = true;

    // By default m_iterEventPointer point to an invalid position
    m_iterEventPointer = m_events.end();
  }


protected:

  /// Sort the list of underlying events
  void Sort() const
  {
    // If the list is not sorted, sort the list
    if ( m_bIsNotSorted )
    {
      const_cast<BaseEventManager *>(this)->DoSort();
    }
  }

  /// Sort the list of underlying events
  void DoSort();

  /// current time
  double m_dTime;

  /// The currently pointed to event
  mutable typename Events::const_iterator m_iterEventPointer;

  /// Does the iterator go forward, or backward?
  ProblemType m_problemType;

  /// The container for the underlying events
  Events m_events;

  /// Has the list been sorted?
  bool m_bIsNotSorted;

}; // classs BaseEventManager


template <class T> 
BaseEventManager<T>::BaseEventManager()
  : m_problemType(ProblemType_Backward)
{
  Init();
}

template <class T> 
void BaseEventManager<T>::AddEvent(const shared_ptr<T>& pEvent)
{
  ASSERT_MSG(m_bIsNotSorted,
             "can't add event when the list is already sorted");

  m_events.push_back(pEvent);
}

template <class T> 
void BaseEventManager<T>::SetCurrentTime(double dTime)
{
  m_dTime = dTime;
}

template <class T> 
bool BaseEventManager<T>::HasEventsNow() const
{
  // make sure the list is sorted
  Sort();

  // Check if all events have been processed
  if ( m_iterEventPointer == m_events.end() )
    return false;

  return numeric::AreTimesEqual( (*m_iterEventPointer)->GetTime(), m_dTime );
}

template <class T> 
const T* BaseEventManager<T>::GetEvent() const
{
  // Make sure the list has been sorted
  Sort();

  // Check if we already iterated through all the events.
  // Valid check even if iterating backwards (see below).
  if ( m_iterEventPointer == m_events.end() )
    return 0;

  // Check if the currently pointed to event occurs at the current time
  if ( numeric::AreTimesEqual( (*m_iterEventPointer)->GetTime(), m_dTime ) )
  {
    // Return the current event, and update the iterator according
    // to the direction the events are traversed
    shared_ptr<T> pEvent = *m_iterEventPointer;

    if ( m_problemType == ProblemType_Forward )
      ++m_iterEventPointer;
    else
    {
      if ( m_iterEventPointer == m_events.begin() )
        m_iterEventPointer = m_events.end();
      else
        --m_iterEventPointer;
    }

    return pEvent.get();
  }

  return 0;
}

template <class T> 
void BaseEventManager<T>::SetInitialState(double dTime)
{
  Sort();

  // Set-up the iterator
  if ( m_problemType == ProblemType_Forward )
    m_iterEventPointer = m_events.begin();
  else  if ( !m_events.empty() )     
    m_iterEventPointer = --(m_events.end());

  SetCurrentTime(dTime);
}

template <class T> 
void BaseEventManager<T>::DoSort() 
{
  m_bIsNotSorted = false;

  // handle special case of no events
  if ( m_events.size() )
    m_events.sort(std::greater< shared_ptr<T> >());
}


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_BASEEVENTMANAGER_H_
