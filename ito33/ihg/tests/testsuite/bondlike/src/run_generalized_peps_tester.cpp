#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "generalized_peps_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern GeneralizedPEPSTester generalized_pepsTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !generalized_pepsTester. ##name().Test(testTags));         \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitGeneralizedPEPSTest : public CppUnit::TestCase
{

public:
  CppUnitGeneralizedPEPSTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitGeneralizedPEPSTest );

    CPPUNIT_TEST( TestGetInitialData );
    CPPUNIT_TEST( TestGetDataWithOptionalConversion );
    CPPUNIT_TEST( TestAddCallExchange);
    CPPUNIT_TEST( TestAddCallFixedCash);
    CPPUNIT_TEST( TestAddCallFixedCashWithNotice);
    CPPUNIT_TEST( TestAddCallFixedCashWithTriggerPeriod);
    CPPUNIT_TEST( TestAddDividends );
    CPPUNIT_TEST( TestAddNewShare );
    CPPUNIT_TEST( TestAddCrossCurrency );
    CPPUNIT_TEST( TestAddNewShareAndCrossCurrency ); 

  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetInitialData );
  IMPLEMENT_TEST( GetDataWithOptionalConversion );
  IMPLEMENT_TEST( AddCallExchange );
  IMPLEMENT_TEST( AddCallFixedCash );
  IMPLEMENT_TEST( AddCallFixedCashWithNotice );
  IMPLEMENT_TEST( AddCallFixedCashWithTriggerPeriod );
  IMPLEMENT_TEST( AddCrossCurrency );
  IMPLEMENT_TEST( AddDividends );
  IMPLEMENT_TEST( AddNewShare );
  IMPLEMENT_TEST( AddNewShareAndCrossCurrency );
};


class CppUnitGeneralizedPEPSAveragingTest : public CppUnit::TestCase
{

public:
  CppUnitGeneralizedPEPSAveragingTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitGeneralizedPEPSAveragingTest );

    CPPUNIT_TEST( TestGetInitialData );

  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetInitialData );
};


int RunGeneralizedPEPSTester(bool bIsAcceptance)
{
  
  const size_t nNbGeneralizedPEPSTests = 3;
  std::vector<size_t> pnIdxGeneralizedPEPSTest;

  if(bIsAcceptance)
  {
    pnIdxGeneralizedPEPSTest.resize(0);
    pnIdxGeneralizedPEPSTest.push_back(1);
  }
  else
  {
    pnIdxGeneralizedPEPSTest.resize(nNbGeneralizedPEPSTests);
    for(size_t n = 0; n < nNbGeneralizedPEPSTests; n++)
      pnIdxGeneralizedPEPSTest[n] = n + 1;
  }

  try
  {
    size_t n;
    std::cout << "\nRuning tests for GeneralizedPEPS...\n";
    g_outputsCB.Read("xmlfiles/output_GeneralizedPEPS.xml");

    
    ///------------------- general generalized_peps test ---------------------  
    for(n = 0; n < pnIdxGeneralizedPEPSTest.size(); n++)
    {
     
      std::cout << "Testing: " <<
        String::Printf("xmlfiles/generalizedpeps%02d.xml", pnIdxGeneralizedPEPSTest[n]).c_str()      
        << std::endl;

      generalized_pepsTester.Setup(
        String::Printf("xmlfiles/generalizedpeps%02d.xml", pnIdxGeneralizedPEPSTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      if ( n == 0 ) 
        runner.addTest( CppUnitGeneralizedPEPSTest::suite() );
      else
        runner.addTest( CppUnitGeneralizedPEPSAveragingTest::suite() );

      runner.run("");

    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_GeneralizedPEPS_new.xml",
                            XML_IHGTEST_RESET_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_GeneralizedPEPS.xml",
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

