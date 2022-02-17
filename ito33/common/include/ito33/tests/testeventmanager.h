#ifndef _ITO33_TEST_EVENT_MANAGER_H_
#define _ITO33_TEST_EVENT_MANAGER_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"

#include "ito33/pricing/eventmanager.h"

// normally there should be no "using" in header but this is just a test...
using namespace ito33;
using namespace ito33::pricing;


// Test class
class EventManagerTest : public CppUnit::TestCase { 

public: 
  EventManagerTest( ) {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( EventManagerTest );
    CPPUNIT_TEST ( TestEmptyForward );
    CPPUNIT_TEST ( TestEmptyBackward );
    CPPUNIT_TEST ( TestSingleForward );
    CPPUNIT_TEST ( TestSingleBackward );
    CPPUNIT_TEST ( TestMultipleBackward );
    CPPUNIT_TEST ( TestMultipleForward );
    CPPUNIT_TEST ( TestEventTimes );
   // CPPUNIT_TEST ( TestApplyToPrice );
    CPPUNIT_TEST ( TestForwardSort );
    CPPUNIT_TEST ( TestBackwardSort );
  CPPUNIT_TEST_SUITE_END();

  // helper function for some of the tests
  void TestEmpty(EventManager& eventManager);
  void TestSingle(EventManager& eventManager);

  // The test functions
  void TestEmptyForward();
  void TestEmptyBackward();
  void TestSingleForward();
  void TestSingleBackward();
  void TestMultipleForward();
  void TestMultipleBackward();
  void TestEventTimes();
  void TestApplyToPrice();
  void TestForwardSort();
  void TestBackwardSort();

  NO_COPY_CLASS( EventManagerTest );
};
#endif
