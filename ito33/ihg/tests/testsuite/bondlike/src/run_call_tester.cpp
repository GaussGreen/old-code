#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"

// local files
#include "call_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern CallTester callTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !callTester. ##name().Test(testTags));         \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitCallTest : public CppUnit::TestCase
{

public:
  CppUnitCallTest() {}


  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitCallTest );
   CPPUNIT_TEST( TestGetBasicData );
   CPPUNIT_TEST( PriceIncreaseAsCallPeriodIncreases );
   CPPUNIT_TEST( PriceDecreaseAsTriggerHistoryIncreases );
  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );

  void PriceDecreaseAsTriggerHistoryIncreases()
  {
    callTester.PriceDecreaseAsTriggerHistoryIncreases();
  }

  void PriceIncreaseAsCallPeriodIncreases()
  {
    callTester.PriceIncreaseAsCallPeriodIncreases();
  }

};


int RunCallTester(bool bIsAcceptance)
{
  const size_t nNbCallTests = 2;
  std::vector<size_t> pnIdxCallTest;

  if(bIsAcceptance)
  {
    pnIdxCallTest.resize(0);
    pnIdxCallTest.push_back(1);
  }
  else
  {
    pnIdxCallTest.resize(nNbCallTests);

    for(size_t n = 0; n < nNbCallTests; n++)
      pnIdxCallTest[n] = n + 1;
  }

  try
  {  
    std::cout << "\nRuning tests for call provision...\n";

    g_outputsCB.Read("xmlfiles/output_call.xml");
  
    size_t n;
    for(n = 0; n < pnIdxCallTest.size(); n++)
    {
            
      std::cout << "Testing: " <<
        String::Printf("xmlfiles/call%02d.xml", pnIdxCallTest[n]).c_str()
        <<std::endl;

      callTester.Setup(
        String::Printf
          ("xmlfiles/call%02d.xml", pnIdxCallTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitCallTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_call_new.xml",
                            XML_IHGTEST_CALL_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_call.xml",
                              XML_IHGTEST_CALL_OUTPUT,
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

