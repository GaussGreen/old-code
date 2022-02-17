/////////////////////////////////////////////////////////////////////////////
// Name:        tests/date/testdate.cpp
// Purpose:     test file for date class
// Author:      Nabil
// Created:     02.05.2006
// RCS-ID:      $Id: testdate.cpp,v 1.1 2006/05/11 10:38:02 nabil Exp $
// Copyright:   (c) 2004-2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"
#include "ito33/vector.h"

#include "ito33/tests/testdate.h"

using namespace ito33;
using namespace std;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( DateTestCase, "DateTestCase" );

#define DCC_30360                   "30360"
#define DCC_30E360                  "30E360"
#define DCC_30U360                  "30U360"
#define DCC_ActAct                  "actact"
//#define DCC_ActAct_ISDA             "actact_ISDA"
#define DCC_Act360                  "act360"
#define DCC_Act365                  "act365"
#define DCC_Act365L                 "act365L"
// Day count conventions with NO EOM (Non End Of Month)
#define DCC_30360_NO_EOM            "30360_NO_EOM"
#define DCC_30E360_NO_EOM           "30E360_NO_EOM"
#define DCC_30U360_NO_EOM           "30U360_NO_EOM"
#define DCC_ActAct_NO_EOM           "actact_NO_EOM"
//#define DCC_ActAct_ISDA_NO_EOM      "actact_ISDA_NO_EOM"
#define DCC_Act360_NO_EOM           "act360_NO_EOM"
#define DCC_Act365_NO_EOM           "act365_NO_EOM"
#define DCC_Act365L_NO_EOM          "act365L_NO_EOM"
#define DCC_Max                     "max"

void DateTestCase::Init()
{
  m_dccValueNames.insert(std::make_pair(DCC_30360, 
    Date::DayCountConvention_30360));
  m_dccValueNames.insert(std::make_pair(DCC_30E360, 
    Date::DayCountConvention_30E360));
  m_dccValueNames.insert(std::make_pair(DCC_30U360, 
    Date::DayCountConvention_30U360));
  m_dccValueNames.insert(std::make_pair(DCC_ActAct, 
    Date::DayCountConvention_ActAct));
  /*m_dccValueNames.insert(std::make_pair(DCC_ActAct_ISDA, 
    Date::DayCountConvention_ActAct_ISDA));*/
  m_dccValueNames.insert(std::make_pair(DCC_Act360, 
    Date::DayCountConvention_Act360));
  m_dccValueNames.insert(std::make_pair(DCC_Act365, 
    Date::DayCountConvention_Act365));
  m_dccValueNames.insert(std::make_pair(DCC_Act365L, 
    Date::DayCountConvention_Act365L));
  
  // Day count conventions with NO EOM (Non End Of Month)
  m_dccValueNames.insert(std::make_pair(DCC_30360_NO_EOM, 
    Date::DayCountConvention_30360_NO_EOM));
  m_dccValueNames.insert(std::make_pair(DCC_30E360_NO_EOM, 
    Date::DayCountConvention_30E360_NO_EOM));
  m_dccValueNames.insert(std::make_pair(DCC_30U360_NO_EOM, 
    Date::DayCountConvention_30U360_NO_EOM));
  m_dccValueNames.insert(std::make_pair(DCC_ActAct_NO_EOM, 
    Date::DayCountConvention_ActAct_NO_EOM));
  /*m_dccValueNames.insert(std::make_pair(DCC_ActAct_ISDA_NO_EOM, 
    Date::DayCountConvention_ActAct_ISDA_NO_EOM));*/
  m_dccValueNames.insert(std::make_pair(DCC_Act360_NO_EOM, 
    Date::DayCountConvention_Act360_NO_EOM));
  m_dccValueNames.insert(std::make_pair(DCC_Act365_NO_EOM, 
    Date::DayCountConvention_Act365_NO_EOM));
  m_dccValueNames.insert(std::make_pair(DCC_Act365L_NO_EOM, 
    Date::DayCountConvention_Act365L_NO_EOM));
  m_dccValueNames.insert(std::make_pair(DCC_Max, 
    Date::DayCountConvention_Max));
}

void DateTestCase::TestDayCountConvention()
{
  Date date1("2006/06/23");
  Date date2("2020/06/23");
  Date::DayCountConvention dcc = Date::DayCountConvention_30360;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(Date::YearsDiff(date2, date1, dcc), -14, 1.e-6);
}

void DateTestCase::TestExcel()
{
    static const unsigned long serials[] =
    {
        // note that day 60 does _not_ work here because it corresponds to a
        // non existing date of Feb 29, 2000 which is internally mapped to Feb
        // 28, 2003
        0, 1, 30, 59, /* 60, */ 61, 30000,
    };

    for ( size_t n = 0; n < SIZEOF(serials); ++n ) {
        Date dt(serials[n]);
        CPPUNIT_ASSERT( dt.GetExcel() == serials[n] );
    }

    // the degenerate case
    Date dt;
    dt.SetExcel(60);
    CPPUNIT_ASSERT( dt.GetExcel() == 59 );
}

#ifdef _WIN32

void DateTestCase::TestOleDate()
{
    // only put integer values here as the Date class doesn't store the time
    // part
    static const struct OleDateTestData
    {
        double dateOle;
        Date::Tm tm;
    } testdata[] =
    {
        { -6666.0,   { 1881, Date::Sep, 29 } },
        { -666.0,    { 1898, Date::Mar,  4 } },
        { -3.0,      { 1899, Date::Dec, 27 } },
        { -2.0,      { 1899, Date::Dec, 28 } },
        { -1.0,      { 1899, Date::Dec, 29 } },
        {  0.0,      { 1899, Date::Dec, 30 } },
        {  0.0,      { 1899, Date::Dec, 30 } },
        {  1.0,      { 1899, Date::Dec, 31 } },
        {  2.0,      { 1900, Date::Jan,  1 } },
        {  3.0,      { 1900, Date::Jan,  2 } },
        {  666.0,    { 1901, Date::Oct, 27 } },
        {  6666.0,   { 1918, Date::Apr,  1 } },
        {  66666.0,  { 2082, Date::Jul,  9 } },

        {  27909.0,  { 1976, Date::May, 29 } },
        {  36585.0,  { 2000, Date::Feb, 29 } },
        {  38046.0,  { 2004, Date::Feb, 29 } },
        {  38176.0,  { 2004, Date::Jul,  8 } },
    };

    Date dt;
    for ( size_t n = 0; n < SIZEOF(testdata); ++n ) {
        const OleDateTestData& data = testdata[n];
        dt.SetTm(data.tm);
        CPPUNIT_ASSERT( dt.GetOleDate() == data.dateOle );

        dt.SetOleDate(data.dateOle);
        CPPUNIT_ASSERT( dt.GetOleDate() == data.dateOle );
    }
}

#endif // _WIN32

void DateTestCase::TestDaysDiffWithDayCount()
{
  struct TestDaysDiffWithDCC
  {
    Date::Tm m_tm1;
    
    Date::Tm m_tm2;

    Date::DayCountConvention m_dcc;

    long m_iExpectedDaysDiff;
  } 
  pTests[] =
  {
    { {2003, Date::Oct, 31}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_30360, 60 },
    { {2003, Date::Oct, 31}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_30E360, 60 },
    { {2003, Date::Oct, 31}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_30U360, 60 },
    { {2003, Date::Oct, 31}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_ActAct, 61 },
    //
    { {2003, Date::Oct, 29}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_30360, 62 },
    { {2003, Date::Oct, 29}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_30E360, 61 },
    { {2003, Date::Oct, 29}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_30U360, 62 },
    { {2003, Date::Oct, 29}, {2003, Date::Dec, 31}, 
      Date::DayCountConvention_ActAct, 63 },
    //
    { {2003, Date::Feb, 28}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_30360, 361 },
    { {2003, Date::Feb, 28}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_30E360, 361 },
    { {2003, Date::Feb, 28}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_30U360, 360 },
    { {2003, Date::Feb, 28}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_ActAct, 366 },
    //
    { {2003, Date::Feb, 27}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_30360, 362},
    { {2003, Date::Feb, 27}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_30E360, 362},
    { {2003, Date::Feb, 27}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_30U360, 362},
    { {2003, Date::Feb, 27}, {2004, Date::Feb, 29}, 
      Date::DayCountConvention_ActAct, 367},
  };

  size_t nNbTests = SIZEOF(pTests);
  
  //Date date;

  long iNbDays;

  Date
    date1,
    date2;

  for(size_t nIdx = 0; nIdx < nNbTests; ++nIdx)
  {
    date1.SetTm(pTests[nIdx].m_tm1);
    date2.SetTm(pTests[nIdx].m_tm2);

    iNbDays = Date::DaysDiffWithDayCount(date1, date2, pTests[nIdx].m_dcc);

    CPPUNIT_ASSERT( pTests[nIdx].m_iExpectedDaysDiff == iNbDays );
  }
}

void DateTestCase::TestYearsDiffWithDayCount()
{
  Date 
    date1("2005/05/31"),
    date2("2006/05/31");

  vector<Date::DayCountConvention> dcc;
  dcc.push_back(Date::DayCountConvention_30360);
  dcc.push_back(Date::DayCountConvention_30E360);
  dcc.push_back(Date::DayCountConvention_30U360);
  dcc.push_back(Date::DayCountConvention_ActAct);
  // NO EOM
  dcc.push_back(Date::DayCountConvention_30360_NO_EOM);
  dcc.push_back(Date::DayCountConvention_30E360_NO_EOM);
  dcc.push_back(Date::DayCountConvention_30U360_NO_EOM);
  dcc.push_back(Date::DayCountConvention_ActAct_NO_EOM);

  double dYears;

  for (vector<Date::DayCountConvention>::const_iterator pdcc = dcc.begin();
       pdcc != dcc.end();
       ++pdcc)
  {
    dYears = Date::YearsDiff(date1, date2, *pdcc);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dYears, 1, 1.e-6);
  }
}
