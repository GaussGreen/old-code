#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "percs_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern PERCSTester percsTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !percsTester. ##name().Test(testTags));         \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitPERCSTest : public CppUnit::TestCase
{

public:
  CppUnitPERCSTest()
  {
  }

  void setUp() {}

  void tearDown() {}

 private:
  CPPUNIT_TEST_SUITE( CppUnitPERCSTest );

    CPPUNIT_TEST( TestGetInitialData );
    CPPUNIT_TEST( TestAddDividends );
    CPPUNIT_TEST( TestAddNewShare );
    CPPUNIT_TEST( TestAddCrossCurrency );
    CPPUNIT_TEST( TestAddNewShareAndCrossCurrency );

  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetInitialData );
  IMPLEMENT_TEST( AddDividends );
  IMPLEMENT_TEST( AddNewShare );
  IMPLEMENT_TEST( AddNewShareAndCrossCurrency );
  IMPLEMENT_TEST( AddCrossCurrency );
};


int RunPERCSTester(bool bIsAcceptance)
{
  const size_t nNbPERCSTests = 3;
  std::vector<size_t> pnIdxPERCSTest;

  if(bIsAcceptance)
  {
    pnIdxPERCSTest.resize(0);
    pnIdxPERCSTest.push_back(1);
  }
  else
  {
    pnIdxPERCSTest.resize(nNbPERCSTests);
    for(size_t n = 0; n < nNbPERCSTests; n++)
      pnIdxPERCSTest[n] = n + 1;
  }
  
  try
  {
    size_t n;
    std::cout << "\nRunning tests for PERCS...\n";
    g_outputsCB.Read("xmlfiles/output_PERCS.xml");

    
    ///------------------- general percs test ---------------------  
    for(n = 0; n < pnIdxPERCSTest.size(); n++)
    {
      percsTester.Setup(
        String::Printf("xmlfiles/percs%02d.xml", (pnIdxPERCSTest[n])).c_str());
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitPERCSTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_PERCS_new.xml",
                            XML_IHGTEST_RESET_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_PERCS.xml",
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

