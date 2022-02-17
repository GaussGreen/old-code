#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "conversion_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern ConversionTester conversionTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !conversionTester. ##name().Test(testTags));   \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitConversionTest : public CppUnit::TestCase
{

public:
  CppUnitConversionTest()
   {
   }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitConversionTest );
    CPPUNIT_TEST( TestAddAnytimeCheckDate );
    CPPUNIT_TEST( TestAddAnytimeForever );
    CPPUNIT_TEST( TestAddQuarterlyForever );
    CPPUNIT_TEST( TestAddQuarterlyNextQuarter );

    CPPUNIT_TEST( TestForcedConversion );
    CPPUNIT_TEST( TestTriggerContinuity );
    
    CPPUNIT_TEST( TestComparisons );
    CPPUNIT_TEST( TestKeepAccrued );

    CPPUNIT_TEST( TestAnytimeForever );

  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( AddAnytimeCheckDate );
  IMPLEMENT_TEST( AddAnytimeForever );
  IMPLEMENT_TEST( AddQuarterlyForever );
  IMPLEMENT_TEST( AddQuarterlyNextQuarter );

  void TestForcedConversion()
  {
    conversionTester.TestForcedConversion();
  }

  void TestTriggerContinuity()
  {
    conversionTester.TestTriggerContinuity();
  }

  void TestComparisons()
  {
    conversionTester.TestComparisons();
  }

  void TestKeepAccrued()
  {
    conversionTester.TestKeepAccrued();
  }

  void TestAnytimeForever()
  {
    conversionTester.TestAnytimeForever();
  }

};


int RunConversionTester(bool bIsAcceptance)
{
  const size_t nNbConversionTests = 3;
  std::vector<size_t> pnIdxConversionTest;

  if(bIsAcceptance)
  {
    pnIdxConversionTest.resize(0);
    pnIdxConversionTest.push_back(1);
  }
  else
  {
    pnIdxConversionTest.resize(nNbConversionTests);
    for(size_t n = 0; n < nNbConversionTests; n++)
      pnIdxConversionTest[n] = n + 1;
  }

  try
  {  
    std::cout << "\nRunning tests for conversion provision...\n";

    g_outputsCB.Read("xmlfiles/output_conversion.xml");
  
    size_t n;
    for(n = 0; n < pnIdxConversionTest.size(); n++)
    {
      conversionTester.Setup(
        String::Printf
          ("xmlfiles/conversion%02d.xml", pnIdxConversionTest[n]).c_str());
      
      std::cout << "Testing: " << 
        String::Printf("xmlfiles/conversion%02d.xml", pnIdxConversionTest[n]).c_str()
        << std::endl;

      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitConversionTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_conversion_new.xml",
                            XML_IHGTEST_CONVERSION_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_conversion.xml",
                              XML_IHGTEST_CONVERSION_OUTPUT,
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

