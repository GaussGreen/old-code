
/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/cppunit_option_test.h
// Purpose:     Base class for testing IHG projects
// Author:      ITO33 Canada
// Created:     2005/06/08
// RCS-ID:      $Id: cppunit_option_test.h,v 1.5 2006/08/22 15:26:31 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/cppunit_option_test.h
**/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_OPTION_TEST_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_OPTION_TEST_H_


#include "ito33/exception.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "optioninterface.h"

namespace ito33 
{
namespace ihg
{
namespace test
{

class CppUnitOptionTest: public CppUnit::TestCase
{
 
private:
  shared_ptr<OptionInterface> m_pOptionInterface;

public:
  CppUnitOptionTest() {}

  void setup();


 private:
  CPPUNIT_TEST_SUITE( CppUnitOptionTest );

   //internal coherence tests
   CPPUNIT_TEST( PriceAlwaysPositive );
   CPPUNIT_TEST( GammaAlwaysPositive );
   CPPUNIT_TEST( CallIncreasesWhenSpotIncreases ); 
   CPPUNIT_TEST( CallDecreasesWhenSpotDecreases );
   CPPUNIT_TEST( PutDecreasesWhenSpotIncreases );  
   CPPUNIT_TEST( PutIncreasesWhenSpotDecreases );
   CPPUNIT_TEST( PriceIncreasesWhenVolatilityIncreases );   
   CPPUNIT_TEST( PriceDecreasesWhenVolatilityDecreases ); 
   CPPUNIT_TEST( AmericanGreaterEuropean );
   CPPUNIT_TEST( PriceIncreasesWhenHazardRateIncreases );   
   CPPUNIT_TEST( PriceIncreasesWhenMaturityIncreases);
   CPPUNIT_TEST( PriceDecreasesWhenRiskFreeRateDecreases );
   CPPUNIT_TEST( CallDecreasesWhenForeignRateIncreases );
   CPPUNIT_TEST( PutIncreasesWhenForeignRateIncreases );
   CPPUNIT_TEST( DeltaInTheMoneyEqualToOne );
   CPPUNIT_TEST( DeltaOutOfTheMoneyEqualToZero );
   CPPUNIT_TEST( GammaOutInOfTheMoneyEqualToZero ); 
   CPPUNIT_TEST( ThetaNegativeRiskFreeRateZero );
   CPPUNIT_TEST( PutDecreasesWhenStrikeDecreases ); 
   CPPUNIT_TEST( CallGoesToSWhenStrikeDecreases); 

   CPPUNIT_TEST_EXCEPTION( PriceAfterMaturityGeneratesAnError , ito33::Exception);

  CPPUNIT_TEST_SUITE_END();

  void PriceAfterMaturityGeneratesAnError();

  //comparison with black scholes
  void SpotPrice();
  void YieldRate();
  void ForeignRate();
  void Volatility();
  void MaturityDate();
  void Strike();

  //coherence testing
  void PriceAlwaysPositive();
  void GammaAlwaysPositive();
  void PutCallParity();
  void CallIncreasesWhenSpotIncreases();
  void CallDecreasesWhenSpotDecreases();
  void PutDecreasesWhenSpotIncreases();
  void PutIncreasesWhenSpotDecreases();
  void PriceIncreasesWhenVolatilityIncreases();   
  void PriceDecreasesWhenVolatilityDecreases(); 
  void PriceIncreasesWhenHazardRateIncreases();
  void PriceIncreasesWhenMaturityIncreases();
  void PriceDecreasesWhenRiskFreeRateDecreases();
  void CallDecreasesWhenForeignRateIncreases();
  void PutIncreasesWhenForeignRateIncreases();
  void DeltaInTheMoneyEqualToOne();
  void DeltaOutOfTheMoneyEqualToZero();
  void AmericanGreaterEuropean();
  void GammaOutInOfTheMoneyEqualToZero();
  void ThetaNegativeRiskFreeRateZero();
  void PutDecreasesWhenStrikeDecreases();
  void CallGoesToSWhenStrikeDecreases();
  

  NO_COPY_CLASS(CppUnitOptionTest);

};


class CppUnitOptionLegacyTests: public CppUnit::TestCase
{
 
public:
  CppUnitOptionLegacyTests() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitOptionLegacyTests );    
    CPPUNIT_TEST( Test3 );
    CPPUNIT_TEST( Test4 );
    CPPUNIT_TEST( Test5 );
    CPPUNIT_TEST( Test6 );
    CPPUNIT_TEST( Test7 );
    CPPUNIT_TEST( Test8 );
    CPPUNIT_TEST( Test9 );
    CPPUNIT_TEST( Test10 );
    CPPUNIT_TEST( Test11 );
    CPPUNIT_TEST( Test12 );
    CPPUNIT_TEST( Test13 );
    CPPUNIT_TEST( Test14 );
    CPPUNIT_TEST( Test15 );
    CPPUNIT_TEST( Test16 );
    CPPUNIT_TEST( Test17 );
    CPPUNIT_TEST( Test18 );
    CPPUNIT_TEST( Test19 );
    CPPUNIT_TEST( Test20 );
    CPPUNIT_TEST( Test21 );
    CPPUNIT_TEST( Test22 );
  CPPUNIT_TEST_SUITE_END();

  void Test3();
  void Test4();
  void Test5();
  void Test6();
  void Test7();
  void Test8();
  void Test9();
  void Test10();
  void Test11();
  void Test12();
  void Test13();
  void Test14();
  void Test15();
  void Test16();
  void Test17();
  void Test18();
  void Test19();
  void Test20();
  void Test21();
  void Test22();

  NO_COPY_CLASS(CppUnitOptionLegacyTests);

};
} //end namespace test
} //end namespace ihg
}//end namespace ito33

#endif
