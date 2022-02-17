#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "warrant_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern WarrantTester warrantTester;

#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !warrantTester. ##name().Test(testTags));      \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitWarrantTest : public CppUnit::TestCase
{

public:
  CppUnitWarrantTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitWarrantTest );
    
    CPPUNIT_TEST( TestStrike );
    CPPUNIT_TEST( TestResetDate ); 
    CPPUNIT_TEST( TestCapRatioIncreases );
    CPPUNIT_TEST( TestShareFactor );
    CPPUNIT_TEST( TestBaseRatio );
    
    CPPUNIT_TEST( TestGetBasicData );
    CPPUNIT_TEST( TestAddCallNoticeToBasicData );
    CPPUNIT_TEST( TestDividendOnResetDate );
    CPPUNIT_TEST( TestIsLastTriggerConditionMet );

    //CPPUNIT_TEST( TestLastTriggerConditionMetProperties );
   // CPPUNIT_TEST( TestOneDvsTwoD );
  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );
  IMPLEMENT_TEST( AddCallNoticeToBasicData );
  IMPLEMENT_TEST( DividendOnResetDate );
  IMPLEMENT_TEST( IsLastTriggerConditionMet );
 
  //Increase the cap ratio, price should not decrease
  void TestCapRatioIncreases()
  {
    warrantTester.TestCapRatioIncreases();
  }

  //as the strike increases the price should not decreases
  void TestStrike()
  {
    warrantTester.TestStrike();
  }
  
  //increase the number of sharefactor price should not decrease
  void TestShareFactor()
  {
    warrantTester.TestShareFactor();
  }
  
  //as the base ratio increase the price increase
  void TestBaseRatio()
  {
    warrantTester.TestBaseRatio();
  }
  
  //shift the maturity date from the valuation date 
  //to the actual reset date minus one
  void TestResetDate()
  {
    warrantTester.TestResetDate();
  }
  
  
  //shift the call date slowly and expect
  //to  have the difference between 1D and 2D
  //to be decreasing
  void TestOneDvsTwoD()
  {
    warrantTester.TestOneDvsTwoD();
  }

  //check different
  //properties when the is last trigger condition met
  void TestLastTriggerConditionMetProperties()
  {
    warrantTester.TestLastTriggerConditionMetProperties();
  }

};



int RunWarrantTester(bool bIsAcceptance)
{
  const size_t nNbWarrantTests = 3;
  std::vector<size_t> pnIdxWarrantTest;
  

  if(bIsAcceptance)
  {
    pnIdxWarrantTest.resize(0);
    pnIdxWarrantTest.push_back(1);
  }
  else
  {
    size_t n;
    pnIdxWarrantTest.resize(nNbWarrantTests);
    for(n = 0; n < nNbWarrantTests; n++)
      pnIdxWarrantTest[n] = n + 1;
  }
  
  try
  {
    size_t n;
    std::cout << "\nRuning tests for Warrant...\n";
    g_outputsCB.Read("xmlfiles/output_Warrant.xml");

   
    for(n = 0; n < pnIdxWarrantTest.size(); n++)
    {
      std::cout << "Testing: " <<
        String::Printf("xmlfiles/warrants%02d.xml",pnIdxWarrantTest[n]).c_str()
        <<std::endl;

      warrantTester.Setup(
        String::Printf("xmlfiles/warrants%02d.xml",
                       pnIdxWarrantTest[n]).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitWarrantTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_Warrant_new.xml",
                            XML_IHGTEST_WARRANT_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_Warrant.xml",
                              XML_IHGTEST_WARRANT_OUTPUT,
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

