/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/eventmanager.h
// Purpose:     Manage events (dividends, coupons, etc)
// Author:      David Pooley
// Created:     2003/10/03
// RCS-ID:      $Id: eventmanager.h,v 1.22 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/eventmanager.h
    @brief Class to manage events (dividends, coupons, etc) that change locally
           the prices
 */

#ifndef _ITO33_PRICING_EVENTMANAGER_H_
#define _ITO33_PRICING_EVENTMANAGER_H_

#include "ito33/beforestd.h"
#include <functional>
#include "ito33/afterstd.h"

#include "ito33/debug.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/baseeventmanager.h"

namespace std
{

/*
   A specialization of greater for shared pointers of Event,
   used for ascending sort of events. 

   Note that this specialization can't be put into source file since the 
   compiler may use only this header file for specialization, then it will give 
   an compile error.

   Even if a raw pointer is used, the compiler may still get wrong since it may use
   the default operator for pointer comparison.
*/
template<>
inline bool greater< ito33::shared_ptr<ito33::pricing::Event> >::operator()
     (const ito33::shared_ptr<ito33::pricing::Event>& pEvent1, 
      const ito33::shared_ptr<ito33::pricing::Event>& pEvent2) const
{
  return *pEvent1 < *pEvent2;
}

} // namespace std

namespace ito33
{

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace pricing
{

enum ProblemType;


/// Class to manage events (dividends, coupons, etc) for all contract types 
class EventManager : public BaseEventManager< Event > 
{
public:
  
  /** 
      EventManager constructor, calls the ctor of the base class
   */
  EventManager() : BaseEventManager<Event>()
  {
  }

  // Default dtor is ok

  /** 
      Get a list of event times as special times for time mesh.

      @param specialTimes special times to be filled
   */
  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

  /** 
      Gets eventually an event occuring at the current time.

      Returns an event occuring at the current time (if any), or NULL if 
      no more events occur at the current time.  Assumes that the list of 
      events is sorted, and that external calls are ordered.
   */
  const Event* GetEvent() const
  {
    const Event* pEvent = BaseEventManager< Event >::GetEvent();

    if (  pEvent == 0 || !pEvent->IsAppliedAfterConstraints()  )
      return pEvent;
    else
    {
      ASSERT_MSG(m_problemType == ProblemType_Backward,
                 "Only backward event can be applied after constraints.");

      // let iterEventPointer point to pEvent
      if ( m_iterEventPointer != m_events.end() )
        ++m_iterEventPointer;
      else
        m_iterEventPointer = m_events.begin();
      return 0;
    }
  }


  /** 
      Gets eventually an event occuring at the current time but should
      be applied after constraints (ResetEventOneD etc.)

      Returns an event occuring at the current time (if any), or NULL if 
      no more events occur at the current time.  Assumes that the list of 
      events is sorted, and that external calls are ordered.
   */
  const Event* GetEventAppliedAfterConstraints() const
  {
    const Event* pEvent = BaseEventManager< Event >::GetEvent();

    ASSERT_MSG( !(pEvent != 0 && !pEvent->IsAppliedAfterConstraints()),
                "A normal basic event is being treated after constraints.\n"\
                "Tips: Check whether all EventAppliedAfterConstraints\n"\
                "has lower order than other basic events.");

    if ( pEvent != 0 && pEvent->IsAppliedAfterConstraints() )
      return pEvent;
    else
      return 0;
  }

}; // class EventManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_EVENTMANAGER_H_
