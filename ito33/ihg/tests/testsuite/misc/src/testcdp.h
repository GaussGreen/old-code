/////////////////////////////////////////////////////////////////////////////
// Name:        testcdp.h
// Purpose:     acceptance tests for cumulative default probabilities
// Created:     13/01/2006
// RCS-ID:      $Id: testcdp.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"


// ----------------------------------------------------------------------------
// Cumulative default probability acceptance tests
// ----------------------------------------------------------------------------

namespace ito33
{

namespace finance
{ 
  class SessionData;
}

namespace ihg
{
class HazardRateTimeOnly;

namespace CDPTEST
{

class CDPTest : public CppUnit::TestCase
{

public:
  
  CDPTest() 
  {
    // Use one session for most (all?) tests
    m_pSessionData = MakeSessionData();

    // Use the same dates for most (all?) tests
    m_pDates = MakeDates();
  }

private:

  CPPUNIT_TEST_SUITE( CDPTest );
    CPPUNIT_TEST( ZeroHR );
    CPPUNIT_TEST( ConstHR );
    CPPUNIT_TEST( TimeOnlyHR1 );
    CPPUNIT_TEST( TimeOnlyHR2 );
    CPPUNIT_TEST( IncreasingHR );
  CPPUNIT_TEST_SUITE_END();

  // zero hazard rate
  void ZeroHR();

  // constant hazard rate
  void ConstHR();

  // time only hazard rate (1 hr date)
  void TimeOnlyHR1();

  // time only hazard rate (several hr dates)
  void TimeOnlyHR2();

  // increasing in value hazard rate
  void IncreasingHR();
  
  // _____________ Helper functions ________________
  
  // Create a session
  shared_ptr<finance::SessionData> MakeSessionData();

  // Create the dates used for testing
  std::vector<Date> MakeDates();

  // Analytic formula for time only hazard rates
  double CDPForTimeOnlyHR(shared_ptr<HazardRateTimeOnly> pHR,
                          double dStartTime,
                          double dEndTime);

  // _____________ Helper variables ________________  

  // The session used testing
  shared_ptr<finance::SessionData> m_pSessionData;

  // The dates used for testing
  std::vector<Date> m_pDates;


  NO_COPY_CLASS(CDPTest);
};

} // namespae CDPTest

} // namespace ihg

} // namespace ito33
