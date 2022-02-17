#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "peps_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern PEPSTester pepsTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !pepsTester. ##name().Test(testTags));         \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitPEPSTest : public CppUnit::TestCase
{

public:
  CppUnitPEPSTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitPEPSTest );

    CPPUNIT_TEST( TestGetInitialData );
    CPPUNIT_TEST( TestGetDataWithOptionalConversion );
    CPPUNIT_TEST( TestAddCallFixedShare);
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
  IMPLEMENT_TEST( AddCallFixedShare );
  IMPLEMENT_TEST( AddCallFixedCash );
  IMPLEMENT_TEST( AddCallFixedCashWithNotice );
  IMPLEMENT_TEST( AddCallFixedCashWithTriggerPeriod );
  IMPLEMENT_TEST( AddDividends );
  IMPLEMENT_TEST( AddNewShare );
  IMPLEMENT_TEST( AddNewShareAndCrossCurrency );
  IMPLEMENT_TEST( AddCrossCurrency );
};



int RunPEPSTester(bool bIsAcceptance)
{
  
  const size_t nNbPEPSTests = 3;
  std::vector<size_t> pnIdxPEPSTest;

  if(bIsAcceptance)
  {
    pnIdxPEPSTest.resize(0);
    pnIdxPEPSTest.push_back(1);
  }
  else
  {
    pnIdxPEPSTest.resize(nNbPEPSTests);
    for(size_t n = 0; n < nNbPEPSTests; n++)
      pnIdxPEPSTest[n] = n + 1;
  }

  try
  {
    size_t n;
    std::cout << "\nRunning tests for PEPS...\n";
    g_outputsCB.Read("xmlfiles/output_PEPS.xml");

    
    ///------------------- general peps test ---------------------  
    for(n = 0; n < pnIdxPEPSTest.size(); n++)
    {
      pepsTester.Setup(
        String::Printf("xmlfiles/peps%02d.xml", pnIdxPEPSTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitPEPSTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_PEPS_new.xml",
                            XML_IHGTEST_RESET_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_PEPS.xml",
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

