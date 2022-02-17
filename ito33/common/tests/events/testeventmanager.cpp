#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/problemtype.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/dividendevents.h"
#include "ito33/pricing/paymentevent.h"

#include "ito33/tests/testeventmanager.h"

using namespace std;
using namespace ito33;
using namespace ito33::pricing;


// No events
void EventManagerTest::TestEmpty(EventManager& eventManager)
{
  const Event* event = eventManager.GetEvent();
  CPPUNIT_ASSERT(event == NULL);

  // all tests should show no events
  numeric::mesh::SpecialTimes specialTimes;
  eventManager.GetSpecialTimes(specialTimes);
  
  CPPUNIT_ASSERT(specialTimes.empty());

  bool bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == false);

  eventManager.SetCurrentTime(0.0);
  bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == false);
}


// No events, backward event manager
void EventManagerTest::TestEmptyBackward()
{
  EventManager eventManager;

  eventManager.SetProblemType(ProblemType_Backward);

  TestEmpty(eventManager);
}

// No events, forward event manager
void EventManagerTest::TestEmptyForward()
{
  EventManager eventManager;

  eventManager.SetProblemType(ProblemType_Forward);

  TestEmpty(eventManager);
}

// Single event (at time 1.0)
void EventManagerTest::TestSingle(EventManager& eventManager)
{
  const Event* event;
  bool bEventsNow;

  // If time is not 1.0, no event should be present
  eventManager.SetCurrentTime(0.0);
  event = eventManager.GetEvent();
  CPPUNIT_ASSERT(event == NULL);

  bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == false);

  eventManager.SetCurrentTime(2.0);
  event = eventManager.GetEvent();
  CPPUNIT_ASSERT(event == NULL);

  bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == false);

  // Check if it can find the event at time 1.0
  eventManager.SetCurrentTime(1.0);

  bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == true);

  event = eventManager.GetEvent();
  CPPUNIT_ASSERT(event != NULL);

  // Once the event has been obtained, no more events should be found
  bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == false);

  // The event should have been processed, so cannot be retrieved
  event = eventManager.GetEvent();
  CPPUNIT_ASSERT(event == NULL);

  // reset the time, but the event should have already been processed
  eventManager.SetCurrentTime(1.0);

  bEventsNow = eventManager.HasEventsNow();
  CPPUNIT_ASSERT(bEventsNow == false);

  event = eventManager.GetEvent();
  CPPUNIT_ASSERT(event == NULL);

  // Check for list of times
  numeric::mesh::SpecialTimes specialTimes;
  eventManager.GetSpecialTimes(specialTimes);
  
  CPPUNIT_ASSERT(specialTimes.size() == 1);
}


// Single event, time 1.0, backward event manager
void EventManagerTest::TestSingleBackward()
{
  EventManager eventManager;
  
  eventManager.SetProblemType(ProblemType_Backward);

  shared_ptr<Event> pevent( new YieldDividendEvent(1.0, 0.05) );
  eventManager.AddEvent(pevent);
  eventManager.SetInitialState(2.0);

  TestSingle(eventManager);
}

// Single event, time 1.0, forward event manager
void EventManagerTest::TestSingleForward()
{
  EventManager eventManager;
  
  eventManager.SetProblemType(ProblemType_Forward);
  
  shared_ptr<Event> pevent( new YieldDividendEvent(1.0, 0.05) );
  eventManager.AddEvent(pevent);
  eventManager.SetInitialState(0.0);

  TestSingle(eventManager);
}

// test multiple events, forward time direction
void EventManagerTest::TestMultipleForward()
{
  EventManager eventManager;

  eventManager.SetProblemType(ProblemType_Forward);
  {
    shared_ptr<Event> pevent;
    pevent = shared_ptr<Event>( new YieldDividendEvent(3.0, 0.03) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.01) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
  }
  
  eventManager.SetInitialState(0.0);


  // Nothing should happen if the current time does not match the event time
  eventManager.SetCurrentTime(0.5);

  const Event* pevent;
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);

  // Start iterating through the events using GetEvent
  eventManager.SetCurrentTime(1.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 1.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);

  // should be 2 events at time 2
  eventManager.SetCurrentTime(2.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 2.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 2.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);

  // only 1 event at time 3
  eventManager.SetCurrentTime(3.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 3.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);
}


void EventManagerTest::TestMultipleBackward()
{
  EventManager eventManager;

  eventManager.SetProblemType(ProblemType_Backward);
  {
    shared_ptr<Event> pevent;
    pevent = shared_ptr<Event>( new YieldDividendEvent(3.0, 0.03) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.01) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
  }
  
  eventManager.SetInitialState(4.0);
    
  // Nothing should happen if the current time does not match the event time
  eventManager.SetCurrentTime(3.5);

  const Event* pevent;

  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);

  // Start iterating through the events using GetEvent
  eventManager.SetCurrentTime(3.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 3.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);

  // Should be 2 events at time 2.0
  eventManager.SetCurrentTime(2.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 2.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 2.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);

  // Only one event at time 1.0
  eventManager.SetCurrentTime(1.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent != NULL);
  CPPUNIT_ASSERT( fabs(pevent->GetTime() - 1.0) < TIMETOLERANCE);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT(pevent == NULL);
}

// Primarily test the GetUniqueEventTimes function
void EventManagerTest::TestEventTimes()
{
  shared_ptr<Event> pevent;
  EventManager eventManager;
  
  eventManager.SetProblemType(ProblemType_Backward);

  // Add a bunch of events, but only 4 unique times
  pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.05) );
  eventManager.AddEvent(pevent);
  pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.06) );
  eventManager.AddEvent(pevent);
  pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.05) );
  eventManager.AddEvent(pevent);
  pevent = shared_ptr<Event>( new YieldDividendEvent(3.0, 0.05) );
  eventManager.AddEvent(pevent);
  pevent = shared_ptr<Event>( new YieldDividendEvent(3.0, 0.05) );
  eventManager.AddEvent(pevent);
  pevent = shared_ptr<Event>( new YieldDividendEvent(3.0, 0.05) );
  eventManager.AddEvent(pevent);
  pevent = shared_ptr<Event>( new YieldDividendEvent(4.0, 0.05) );
  eventManager.AddEvent(pevent);

  double pdTimes[] = { 2.0, 2.0, 1.0, 3.0, 3.0, 3.0, 4.0};

  numeric::mesh::SpecialTimes specialTimes;
  eventManager.GetSpecialTimes(specialTimes);

  CPPUNIT_ASSERT(specialTimes.size() == 7);
  
  // make sure the unique time list is sorted
  numeric::mesh::SpecialTimes::iterator iter;
  
  size_t nIdx = 0;
  for (iter = specialTimes.begin(); iter != specialTimes.end(); ++iter)
    CPPUNIT_ASSERT( numeric::AreTimesEqual(iter->GetTime(), pdTimes[nIdx++]) );
}

//// Primarily test the ApplyToPrice function
//void EventManagerTest::TestApplyToPrice()
//{
//  shared_ptr<Event> pevent;
//
//  EventManager eventManager;
//
//  eventManager.SetProblemType(ProblemType_Backward);
//
//  pevent = shared_ptr<Event>( new YieldDividendEvent(3.0, 0.05) );
//  eventManager.AddEvent(pevent);
//  pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.08) );
//  eventManager.AddEvent(pevent);
//  pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.05) );
//  eventManager.AddEvent(pevent);
//  pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.08) );
//  eventManager.AddEvent(pevent);
//
//  // Set the initial state of the event manager
//  eventManager.SetInitialState(5.0);
//
//  // Nothing should happen if the current time does not match the event time
//  eventManager.SetCurrentTime(4.0);
//
//  double pdS[4] = {1.0, 2.0, 3.0, 4.0};
//  double pdValues[4] = {10.0, 20.0, 30.0, 40.0};
//
//  eventManager.ApplyToPrice(pdS, pdValues, 4);
//
//  size_t nIdx;
//  for (nIdx = 0; nIdx < 4; nIdx++)
//    CPPUNIT_ASSERT( fabs(pdValues[nIdx] - (nIdx*10.0+10.0) ) < 1.e-14 );
//
//  // If the current time matches the event times, then the event should be 
//  // processed
//  eventManager.SetCurrentTime(3.0);
//
//  eventManager.ApplyToPrice(pdS, pdValues, 4);
//
//  for (nIdx = 0; nIdx < 4; nIdx++)
//    CPPUNIT_ASSERT( fabs(pdValues[nIdx] - (nIdx*10.0+10.0) ) > 1.e-14 );
//
//  // If the current time matches the event time, then the event should be 
//  // processed (different event than before)
//  eventManager.SetCurrentTime(2.0);
//  double pdValues2[4] = {10.0, 20.0, 30.0, 40.0};
//
//  eventManager.ApplyToPrice(pdS, pdValues2, 4);
//
//  for (nIdx = 0; nIdx < 4; nIdx++)
//  {
//    CPPUNIT_ASSERT( fabs(pdValues2[nIdx] - (nIdx*10.0+10.0) ) > 1.e-14 );
//    CPPUNIT_ASSERT( fabs(pdValues2[nIdx] - pdValues[nIdx] ) > 1.e-14 );
//  }
//
//  // make sure that two events at the same time are handled
//  eventManager.SetCurrentTime(1.0);
//  double pdValues3[4] = {10.0, 20.0, 30.0, 40.0};
//
//  eventManager.ApplyToPrice(pdS, pdValues3, 4);
//
//  for (nIdx = 0; nIdx < 4; nIdx++)
//  {
//    CPPUNIT_ASSERT( fabs(pdValues3[nIdx] - (nIdx*10.0+10.0) ) > 1.e-14 );
//    CPPUNIT_ASSERT( fabs(pdValues3[nIdx] - pdValues2[nIdx] ) > 1.e-14 );
//    CPPUNIT_ASSERT( fabs(pdValues3[nIdx] - pdValues[nIdx] ) > 1.e-14 );
//  }
//}

// Test the sorting of events, forward time direction
void EventManagerTest::TestForwardSort()
{
  EventManager eventManager;

  eventManager.SetProblemType(ProblemType_Forward);

  {
    // Add several events in random order, including two different
    // events (different types) at the same time
    shared_ptr<Event> pevent;
    pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.01) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new PaymentEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(4.0, 0.03) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new PaymentEvent(3.0, 0.03) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new PaymentEvent(4.0, 0.04) );
    eventManager.AddEvent(pevent);
  }

  
  const Event* pevent;
  // go through the events using the GetEvent function (and
  // setting times). In the forward direction, dividends come
  // before payments.
  eventManager.SetInitialState(0.1);
  eventManager.SetCurrentTime(1.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Dividend );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );

  eventManager.SetCurrentTime(2.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Dividend );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Payment );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );

  eventManager.SetCurrentTime(3.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Payment );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );
 
  eventManager.SetCurrentTime(4.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Dividend );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Payment );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );

}

// Test the sorting of events, backward time direction
void EventManagerTest::TestBackwardSort()
{
  EventManager eventManager;

  eventManager.SetProblemType(ProblemType_Backward);

  {
    // Add several events in random order, including two different
    // events (different types) at the same time. In the backward
    // direction, payments come before dividends.
    shared_ptr<Event> pevent;
    pevent = shared_ptr<Event>( new YieldDividendEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(1.0, 0.01) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new PaymentEvent(2.0, 0.02) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new YieldDividendEvent(4.0, 0.03) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new PaymentEvent(3.0, 0.03) );
    eventManager.AddEvent(pevent);
    pevent = shared_ptr<Event>( new PaymentEvent(4.0, 0.04) );
    eventManager.AddEvent(pevent);
  }

  const Event* pevent;
  // go through the events using the GetEvent function (and
  // setting times)
  eventManager.SetInitialState(4.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Payment );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Dividend );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );

  eventManager.SetCurrentTime(3.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Payment );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );
 
  eventManager.SetCurrentTime(2.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Payment );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Dividend );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );

  eventManager.SetCurrentTime(1.0);
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent != NULL );
  CPPUNIT_ASSERT( pevent->GetType() == ET_Dividend );
  pevent = eventManager.GetEvent();
  CPPUNIT_ASSERT( pevent == NULL );

}
