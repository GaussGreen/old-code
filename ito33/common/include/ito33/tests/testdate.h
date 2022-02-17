/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testdate.h
// Purpose:     header file for date test
// Author:      Nabil 
// Created:     02.05.2006
// RCS-ID:      $Id: testdate.h,v 1.1 2006/05/11 10:36:49 nabil Exp $
// Copyright:   (c) 2004-2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <string>
#include <map>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/date.h"
#include "ito33/exception.h"
#include "ito33/vector.h"

// ----------------------------------------------------------------------------
// test class
// ----------------------------------------------------------------------------

class DateTestCase : public CppUnit::TestCase
{
public:
  DateTestCase() { Init(); }

private:
  CPPUNIT_TEST_SUITE( DateTestCase );
    
    CPPUNIT_TEST( TestDayCountConvention );
    CPPUNIT_TEST( TestExcel );
    CPPUNIT_TEST( TestDaysDiffWithDayCount );
    CPPUNIT_TEST( TestYearsDiffWithDayCount );
    
#ifdef _WIN32
    CPPUNIT_TEST( TestOleDate );
#endif // _WIN32
  CPPUNIT_TEST_SUITE_END();

  void TestDayCountConvention();

  void TestExcel();

  void TestDaysDiffWithDayCount();

  void TestYearsDiffWithDayCount();

#ifdef _WIN32
  void TestOleDate();
#endif // _WIN32
  
  void Init();

  //std::vector<TestDaysDiffWithDCC> 
  //  ReadTestDaysDiffWithDCC(std::string& sDaysDiffWithDCCFile);

  std::map<std::string, ito33::Date::DayCountConvention> m_dccValueNames;

  NO_COPY_CLASS(DateTestCase);
}; //DateTestCase
