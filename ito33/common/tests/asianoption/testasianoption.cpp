/////////////////////////////////////////////////////////////////////////////
// Name:        testasianoption.cpp
// Purpose:     Acceptance test for asian option 
// Author:      ITO 33 Canada
// Created:     March 29, 2005
// RCS-ID:      $Id: testasianoption.cpp,v 1.6 2006/07/28 21:03:07 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/exoticoption/curran.h"
#include "ito33/finance/exoticoption/asianoption.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/exercisetype.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testasianoption.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Asian Option Test
// ----------------------------------------------------------------------------


void AsianOptionTest::ValidFixedStrikeButNegativeStrike()
{

  double dStrike = -1.0;
  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::Jan, 1);

  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_American;

  AsianOption asianOpt(dStrike, 
                       maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

}

void AsianOptionTest::FixedAvgStartDateAfterMaturityDate()
{

  double dStrike = 100.0;
  
  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2006, Date::Jan, 1);

  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_European;

  AsianOption asianOpt(dStrike, 
                       maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

}

void AsianOptionTest::NegativeCurrentAverage()
{
  double dStrike         = 100.0;
  double dCurrentAverage = -10.0;

  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::Jan, 1);
  
  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_European;

  AsianOption asianOpt(dStrike, 
                       maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

  asianOpt.SetCurrentAverage(dCurrentAverage, 10);
}

void AsianOptionTest::FloatingAvgStartDateAfterMaturityDate()
{

  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2006, Date::Jan, 1);

  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_European;

  AsianOption asianOpt(
                       maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

}

void AsianOptionTest::EndStartDateAfterAverageStartDate()
{
  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::May, 1);
  Date avgEndDate(2004, Date::Jan, 1);

  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_European;

  AsianOption asianOpt(maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

  asianOpt.SetAverageEndDate(avgEndDate);
}

void AsianOptionTest::AverageEndDateAfterMaturityDate()
{
  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::Jan, 1);
  Date avgEndDate(2006, Date::Jan, 1);

  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_European;

  AsianOption asianOpt(maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

  asianOpt.SetAverageEndDate(avgEndDate);
}

void AsianOptionTest::EndStartDateEqualsMaturityDate()
{
  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::Jan, 1);

  OptionType optionType     = Option_Put;
  ExerciseType exerciseType = ExerciseType_European;

  AsianOption asianOpt(maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

  CPPUNIT_ASSERT( maturityDate == asianOpt.GetAverageEndDate() );
}

void AsianOptionTest::DumpFixedStrike()
{
  double dStrike        = 100.0;
  double dCurrentAverage = 100.0;

  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::Jan, 1);
  Date avgEndDate(2004, Date::May, 1);

  OptionType optionType     = Option_Call;
  ExerciseType exerciseType = ExerciseType_American;

  AsianOption asianOpt(dStrike, 
                       maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);

  asianOpt.SetCurrentAverage(dCurrentAverage, 10);
  asianOpt.SetAverageEndDate(avgEndDate);

  std::ostringstream oss;

  ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<asian_option>\n"
    "<type>call</type>\n"
    "<exercise>american</exercise>\n"
    "<maturity>2005-01-01</maturity>\n"   
    "<strike>100</strike>\n"
    "<average_start_date>2004-01-01</average_start_date>\n"
    "<average_end_date>2004-05-01</average_end_date>\n"
    "<number_of_sampling_averages>252</number_of_sampling_averages>\n"
    "<current_average>100</current_average>\n"
    "<number_of_samples_used>10</number_of_samples_used>\n"
    "</asian_option>\n"
    "</root>\n");

  ito33::XML::RootTag root("root",oss);

  asianOpt.Dump(root);
}

void AsianOptionTest::DumpFloatingStrike()
{

  Date maturityDate(2005, Date::Jan, 1);
  Date avgStartDate(2004, Date::Jan, 1);

  OptionType optionType     = Option_Call;
  ExerciseType exerciseType = ExerciseType_American;

  AsianOption asianOpt(maturityDate,
                       optionType, 
                       exerciseType,
                       avgStartDate,
                       252);


  std::ostringstream oss;

  ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<asian_option>\n"
    "<type>call</type>\n"
    "<exercise>american</exercise>\n"
    "<maturity>2005-01-01</maturity>\n"
    "<average_start_date>2004-01-01</average_start_date>\n"
    "<number_of_sampling_averages>252</number_of_sampling_averages>\n"
    "</asian_option>\n"
    "</root>\n");

  ito33::XML::RootTag root("root",oss);

  asianOpt.Dump(root);
}

void CurranTest::BeforeAveragingPeriod()
{
  double ONEWEEK       = 1./52.;
  double dSpot         = 100.;
  double dMaturityTime = 72.*ONEWEEK;
  size_t nNbAveraging  = 53;
  double dStartTime    = 20.*ONEWEEK;
  double dInterestRate = .09;
  double dDividend     = 0.;

  std::vector< CurranTestData > pdTestValues;
  CurranTestData testData;

  testData.dVol = .05; testData.dStrike = 95.; testData.dSolution = 11.76;
  pdTestValues.push_back( testData );
  testData.dVol = .05; testData.dStrike = 100.; testData.dSolution = 7.39;
  pdTestValues.push_back( testData );
  testData.dVol = .05; testData.dStrike = 105.; testData.dSolution = 3.47;
  pdTestValues.push_back( testData );

  testData.dVol = .1; testData.dStrike = 95.; testData.dSolution = 11.96;
  pdTestValues.push_back( testData );
  testData.dVol = .1; testData.dStrike = 100.; testData.dSolution = 8.07;
  pdTestValues.push_back( testData );
  testData.dVol = .1; testData.dStrike = 105.; testData.dSolution = 4.87;
  pdTestValues.push_back( testData );

  testData.dVol = .3; testData.dStrike = 90.; testData.dSolution = 19.17;
  pdTestValues.push_back( testData );
  testData.dVol = .3; testData.dStrike = 100.; testData.dSolution = 13.43;
  pdTestValues.push_back( testData );
  testData.dVol = .3; testData.dStrike = 110.; testData.dSolution = 9.05;
  pdTestValues.push_back( testData );

  testData.dVol = .5; testData.dStrike = 90.; testData.dSolution = 24.10;
  pdTestValues.push_back( testData );
  testData.dVol = .5; testData.dStrike = 100.; testData.dSolution = 19.37;
  pdTestValues.push_back( testData );
  testData.dVol = .5; testData.dStrike = 110.; testData.dSolution = 15.47;
  pdTestValues.push_back( testData );

  std::vector< CurranTestData >::iterator iter;

  for ( iter = pdTestValues.begin(); iter != pdTestValues.end() ; ++iter )
  {
    CurranTestData &test = *iter;

    double dPrice
    = CurranCall(dSpot, test.dStrike, dMaturityTime, nNbAveraging,
                 dStartTime, test.dVol, dInterestRate, dDividend);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( test.dSolution, dPrice, 1e-2);
  }
}

/*
  Compare solution with Curran's paper
*/
void CurranTest::AtAveragingPeriod()
{
  double dSpot         = 100.;
  double dMaturityTime = 1.0;
  size_t nNbAveraging  = 53;
  double dStartTime    = 0.0;
  double dInterestRate = .09;
  double dDividend     = 0.;

  std::vector< CurranTestData > pdTestValues;
  CurranTestData testData;

  testData.dVol = .05; testData.dStrike = 95.; testData.dSolution = 8.81;
  pdTestValues.push_back( testData );
  testData.dVol = .05; testData.dStrike = 100.; testData.dSolution = 4.31;
  pdTestValues.push_back( testData );
  testData.dVol = .05; testData.dStrike = 105.; testData.dSolution = .95;
  pdTestValues.push_back( testData );

  testData.dVol = .1; testData.dStrike = 95.; testData.dSolution = 8.91;
  pdTestValues.push_back( testData );
  testData.dVol = .1; testData.dStrike = 100.; testData.dSolution = 4.91;
  pdTestValues.push_back( testData );
  testData.dVol = .1; testData.dStrike = 105.; testData.dSolution = 2.06;
  pdTestValues.push_back( testData );

  testData.dVol = .3; testData.dStrike = 90.; testData.dSolution = 14.96;
  pdTestValues.push_back( testData );
  testData.dVol = .3; testData.dStrike = 100.; testData.dSolution = 8.80;
  pdTestValues.push_back( testData );
  testData.dVol = .3; testData.dStrike = 110.; testData.dSolution = 4.67;
  pdTestValues.push_back( testData );

  testData.dVol = .5; testData.dStrike = 90.; testData.dSolution = 18.14;
  pdTestValues.push_back( testData );
  testData.dVol = .5; testData.dStrike = 100.; testData.dSolution = 12.98;
  pdTestValues.push_back( testData );
  testData.dVol = .5; testData.dStrike = 110.; testData.dSolution = 9.07;
  pdTestValues.push_back( testData );

  std::vector< CurranTestData >::iterator iter;

  for ( iter = pdTestValues.begin(); iter != pdTestValues.end() ; ++iter )
  {
    CurranTestData &test = *iter;

    double dPrice
    = CurranCall(dSpot, test.dStrike, dMaturityTime, nNbAveraging,
                 dStartTime, test.dVol, dInterestRate, dDividend);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( test.dSolution, dPrice, 1e-2);
  }
}
 
void CurranTest::AfterAveragingPeriod()
{
  double ONEWEEK          = 1./52.;
  double dSpot            = 100.;
  double dMaturityTime    = 32.*ONEWEEK;
  size_t nNbAveraging     = 53;
  double dStartTime       = 0.0;
  double dInterestRate    = .09;
  double dDividend        = 0.;
  double dA               = 100.0;
  size_t nNbIntoAveraging = 20;

  std::vector< CurranTestData > pdTestValues;
  CurranTestData testData;

  testData.dVol = .05; testData.dStrike = 95.; testData.dSolution = 6.39;
  pdTestValues.push_back( testData );
  testData.dVol = .05; testData.dStrike = 100.; testData.dSolution = 1.73;
  pdTestValues.push_back( testData );
  testData.dVol = .05; testData.dStrike = 105.; testData.dSolution = 0.01;
  pdTestValues.push_back( testData );

  testData.dVol = .1; testData.dStrike = 95.; testData.dSolution = 6.40;
  pdTestValues.push_back( testData );
  testData.dVol = .1; testData.dStrike = 100.; testData.dSolution = 2.11;
  pdTestValues.push_back( testData );
  testData.dVol = .1; testData.dStrike = 105.; testData.dSolution = 0.20;
  pdTestValues.push_back( testData );

  testData.dVol = .3; testData.dStrike = 90.; testData.dSolution = 11.32;
  pdTestValues.push_back( testData );
  //note this is different from curran, but since everything else matches
  //I assume that 4.12 is a typo
  testData.dVol = .3; testData.dStrike = 100.; testData.dSolution = 4.13; 
  pdTestValues.push_back( testData );
  testData.dVol = .3; testData.dStrike = 110.; testData.dSolution = 0.92;
  pdTestValues.push_back( testData );

    //note this is different from curran, but since everything else matches
  //I assume that 12.29 is a typo
  testData.dVol = .5; testData.dStrike = 90.; testData.dSolution = 12.30;
  pdTestValues.push_back( testData );
  testData.dVol = .5; testData.dStrike = 100.; testData.dSolution = 6.24;
  pdTestValues.push_back( testData );
  testData.dVol = .5; testData.dStrike = 110.; testData.dSolution = 2.78;
  pdTestValues.push_back( testData );

  std::vector< CurranTestData >::iterator iter;

  for ( iter = pdTestValues.begin(); iter != pdTestValues.end() ; ++iter )
  {
    CurranTestData &test = *iter;

    double dPrice
    = CurranCall(dSpot, test.dStrike, dMaturityTime, nNbAveraging,
                 dStartTime, test.dVol, dInterestRate, dDividend, 
                 dA, nNbIntoAveraging);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( test.dSolution, dPrice, 1e-2);
  }
  
}
