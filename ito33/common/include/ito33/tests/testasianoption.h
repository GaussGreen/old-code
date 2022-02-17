/////////////////////////////////////////////////////////////////////////////
// Name:        testasinaoption.h
// Purpose:     acceptance test for asian option
// Author:      ITO 33 Canada
// Created:     March 29, 2005
// RCS-ID:      $Id: testasianoption.h,v 1.3 2006/06/03 16:24:05 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class AsianOptionTest : public CppUnit::TestCase
{

public:
  AsianOptionTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( AsianOptionTest );
    
    //Fixed Strike Asian option

    CPPUNIT_TEST_EXCEPTION(ValidFixedStrikeButNegativeStrike, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION(FixedAvgStartDateAfterMaturityDate,ito33::Exception);
  

    CPPUNIT_TEST( DumpFixedStrike );

    //Floating Strike Asian option

    //Request the strike for a floating
    CPPUNIT_TEST_EXCEPTION(FloatingAvgStartDateAfterMaturityDate, ito33::Exception);

    CPPUNIT_TEST( DumpFloatingStrike );

    //common test
    CPPUNIT_TEST_EXCEPTION(NegativeCurrentAverage,ito33::Exception);
    CPPUNIT_TEST_EXCEPTION(EndStartDateAfterAverageStartDate, ito33::Exception);

    CPPUNIT_TEST ( EndStartDateEqualsMaturityDate  );
  CPPUNIT_TEST_SUITE_END();

  void ValidFixedStrikeButNegativeStrike();
  void FixedAvgStartDateAfterMaturityDate();
  void DumpFixedStrike();

  void GetFloatingStrike();
  void FloatingAvgStartDateAfterMaturityDate();
  void DumpFloatingStrike();

  void NegativeCurrentAverage();
  void EndStartDateAfterAverageStartDate();
  void AverageEndDateAfterMaturityDate();
  void EndStartDateEqualsMaturityDate();


  NO_COPY_CLASS(AsianOptionTest);
};



/*
  Test to reproduce the numbers in Michael curran's paper
  "Beyond Average Intelligence"
*/
class CurranTest: public CppUnit::TestCase
{
public:
  CurranTest() {}

private:
  CPPUNIT_TEST_SUITE( CurranTest );
    CPPUNIT_TEST( BeforeAveragingPeriod );
    CPPUNIT_TEST( AtAveragingPeriod );
    CPPUNIT_TEST( AfterAveragingPeriod );
  CPPUNIT_TEST_SUITE_END( );

  void BeforeAveragingPeriod();
  void AtAveragingPeriod();
  void AfterAveragingPeriod();

  struct CurranTestData
  {
    double dVol;
    double dStrike;
    double dSolution;
  };

  NO_COPY_CLASS( CurranTest );

};
