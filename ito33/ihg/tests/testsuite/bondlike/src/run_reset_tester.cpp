#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "reset_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern ResetTester resetTester;
extern CBUsingResetPricerTester resetTester100;
extern Reset1DModelTester reset1DTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !resetTester. ##name().Test(testTags));   \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitResetTest : public CppUnit::TestCase
{

public:
  CppUnitResetTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitResetTest );
    CPPUNIT_TEST( TestGetBasicData );
    CPPUNIT_TEST( TestAddDividends );
    CPPUNIT_TEST( TestAddNewShare );
    CPPUNIT_TEST( TestAddCrossCurrency );
    CPPUNIT_TEST( TestAddNewShareAndCrossCurrency );
    CPPUNIT_TEST( PriceDecreaseCapWhenNumberOfResetIncreases );
    CPPUNIT_TEST( PriceIncreaseFloorWhenNumberOfResetIncreases );
    CPPUNIT_TEST( PriceStaySameWhenFloorAndCapAreOne );
    CPPUNIT_TEST( PriceStaySameWithResetDateBeforeValuationDate );
  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );
  IMPLEMENT_TEST( AddDividends );
  IMPLEMENT_TEST( AddNewShare );
  IMPLEMENT_TEST( AddCrossCurrency );
  IMPLEMENT_TEST( AddNewShareAndCrossCurrency );

  void PriceDecreaseCapWhenNumberOfResetIncreases()
  {
   resetTester.PriceDecreaseCapWhenNumberOfResetIncreases();
  }

  void PriceIncreaseFloorWhenNumberOfResetIncreases()
  {
    resetTester.PriceIncreaseFloorWhenNumberOfResetIncreases();
  }

  void PriceStaySameWhenFloorAndCapAreOne()
  {
    resetTester.PriceStaySameWhenFloorAndCapAreOne();
  }

  void PriceStaySameWithResetDateBeforeValuationDate()
  {
    resetTester.PriceStaySameWithResetDateBeforeValuationDate();
  }
};


class CppUnitReset100Test : public CppUnit::TestCase
{

public:
  CppUnitReset100Test()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitReset100Test );

    CPPUNIT_TEST( tt );


  CPPUNIT_TEST_SUITE_END();

private:

  void tt()
  {
    CPPUNIT_ASSERT( !resetTester100.GetBasicData().Test(testTags));
    CPPUNIT_ASSERT( resetTester100.Compare() );
  }
};

class CppUnitReset1DTest : public CppUnit::TestCase
{

public:
  CppUnitReset1DTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitReset1DTest );

    CPPUNIT_TEST( tt );


  CPPUNIT_TEST_SUITE_END();

private:

  void tt()
  {
    CPPUNIT_ASSERT( !reset1DTester.GetBasicData().Test(testTags));
    CPPUNIT_ASSERT( reset1DTester.Compare() );
  }

};

int RunResetTester(bool bIsAcceptance)
{
  const size_t nNbResetCBTests = 3;
  std::vector<size_t> pnIdxResetCBTest;
  const size_t nNbReset1DTests = 6;
  std::vector<size_t> pnIdxReset1DTest;
  const size_t nNbResetGeneralTests = 3;
  std::vector<size_t> pnIdxResetGeneralTest;

  if(bIsAcceptance)
  {
    pnIdxResetCBTest.resize(0);
    pnIdxResetCBTest.push_back(1);
    pnIdxReset1DTest.resize(0);
    pnIdxReset1DTest.push_back(1);
    pnIdxResetGeneralTest.resize(0);
    pnIdxResetGeneralTest.push_back(1);
  }
  else
  {
    size_t n;
    pnIdxResetCBTest.resize(nNbResetCBTests);
    for(n = 0; n < nNbResetCBTests; n++)
      pnIdxResetCBTest[n] = n;

    pnIdxReset1DTest.resize(nNbReset1DTests);
    for(n = 0; n < nNbReset1DTests; n++)
      pnIdxReset1DTest[n] = n + 1;
    
    pnIdxResetGeneralTest.resize(nNbResetGeneralTests);
    for(n = 0; n < nNbResetGeneralTests; n++)
      pnIdxResetGeneralTest[n] = n + 1;
  }
  
  try
  {
    size_t n;
    std::cout << "\nRunning tests for Reset...\n";
    g_outputsCB.Read("xmlfiles/output_Reset.xml");

    ///----------------- use reset pricer to price a cb ----------
    for(n = 0; n < pnIdxResetCBTest.size(); n++)
    {
      resetTester100.Setup(
        String::Printf("xmlfiles/reset1%02d.xml",
                       pnIdxResetCBTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitReset100Test::suite() );
      runner.run("");
    }

    ///----------------- test 1D reset pricer ----------
    for(n = 0; n < pnIdxReset1DTest.size(); n++)
    {
      reset1DTester.Setup(
        String::Printf("xmlfiles/reset1D%02d.xml",
                       pnIdxReset1DTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitReset1DTest::suite() );
      runner.run("");
    }

    ///------------------- general reset test ---------------------  
    for(n = 0; n < pnIdxResetGeneralTest.size(); n++)
    {
      resetTester.Setup(
        String::Printf("xmlfiles/reset%02d.xml",
                       pnIdxResetGeneralTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitResetTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_Reset_new.xml",
                            XML_IHGTEST_RESET_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_Reset.xml",
                              XML_IHGTEST_RESET_OUTPUT,
                              ITO33_XML_ATTR_ROOT_VERSION,
                              ITO33_IHG_VERSION_DOT_STRING);
      }
    }
    
    g_nNbFailedSelfRegressionTests += g_outputsCB.GetNbBasicOutputChange();

    return 0; 

   }
  catch ( ito33::Exception& e )
  {
    printf("ITO33 exception:\n%s\n", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    printf("std exception: %s\n", e.what());
  }
  catch ( ... )
  {
    puts("unknown exception!");
  }

  return 0;
}

