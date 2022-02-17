
/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/cppunit_asianoption_test.h
// Purpose:     Base class for testing IHG projects
// Author:      ITO 33 Canada
// Created:     2005/06/08
// RCS-ID:      $Id: cppunit_asianoption_test.h,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_ASIANOPTION_TEST_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_ASIANOPTION_TEST_H_


#include "ito33/exception.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/optiontype.h"

#include "asianoptioninterface.h"

namespace ito33 
{
namespace ihg
{
namespace test
{

class CppUnitAsianOptionTest: public CppUnit::TestCase
{
 
private:
  shared_ptr<AsianOptionInterface> m_pAsianOptionInterface;

public:
  CppUnitAsianOptionTest() {}

  void setup();


 private:
  CPPUNIT_TEST_SUITE( CppUnitAsianOptionTest );
  
   CPPUNIT_TEST( PriceAlwaysPositive );
   CPPUNIT_TEST_EXCEPTION( PriceAfterMaturityGeneratesAnError , ito33::Exception);
   CPPUNIT_TEST( AmericanGreaterEuropean );
   CPPUNIT_TEST( GammaAlwaysPositive );
   CPPUNIT_TEST( PriceIncreasesWhenVolatilityIncreases );   
   CPPUNIT_TEST( PriceDecreasesWhenVolatilityDecreases ); 
   CPPUNIT_TEST( PutDecreasesWhenStrikeDecreases ); 

  CPPUNIT_TEST_SUITE_END();

  void PriceAfterMaturityGeneratesAnError();
  void PriceAlwaysPositive();
  void AmericanGreaterEuropean();
  void GammaAlwaysPositive();
  void PriceIncreasesWhenVolatilityIncreases();
  void PriceDecreasesWhenVolatilityDecreases();
  void PutDecreasesWhenStrikeDecreases();
  void PriceIncreasesWhenHazardRateIncreases();
  void PriceIncreasesWhenMaturityIncreases();


  NO_COPY_CLASS(CppUnitAsianOptionTest);

};

class Curran
{
public: 
  double m_dVol;
  double m_dStrike;
  double m_dSolution;

  Curran(double dVol, double dStrike, double dSolution): 
       m_dVol(dVol), 
       m_dStrike(dStrike), 
       m_dSolution(dSolution) 
       {}


};

class CppUnitCurranAsianOptionTest: public CppUnit::TestCase
{
 

public:
  CppUnitCurranAsianOptionTest() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitCurranAsianOptionTest );
   CPPUNIT_TEST( CheckCurranBeforeAveragingPeriod );
   CPPUNIT_TEST( CheckCurranIntoAveragingPeriod );
   CPPUNIT_TEST( CheckCurranAtAveragingPeriod );
  CPPUNIT_TEST_SUITE_END();

  void CheckCurranBeforeAveragingPeriod();
  void CheckCurranIntoAveragingPeriod();
  void CheckCurranAtAveragingPeriod();

  NO_COPY_CLASS(CppUnitCurranAsianOptionTest);

};


class CppUnitAsianOptionSpecificTest: public CppUnit::TestCase
{
 
public:
  CppUnitAsianOptionSpecificTest() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitAsianOptionSpecificTest );
    CPPUNIT_TEST( ChangeCurrentAverageFixedCall );
    CPPUNIT_TEST( ChangeCurrentAverageFixedPut );
    CPPUNIT_TEST( ChangeCurrentAverageFloatingCall );
    CPPUNIT_TEST( ChangeCurrentAverageFloatingPut );
  CPPUNIT_TEST_SUITE_END();

  void ChangeCurrentAverageFixedCall();
  void ChangeCurrentAverageFixedPut();
  void ChangeCurrentAverageFloatingCall();
  void ChangeCurrentAverageFloatingPut();



  shared_ptr<finance::SessionData> InitSessionData(Date valuationDate);
  void ChangeCurrentAverage(finance::OptionType optionType, bool bIsFloating);

  NO_COPY_CLASS(CppUnitAsianOptionSpecificTest);

};


} //end namespace test
} //end namespace ihg
}//end namespace ito33

#endif
