#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/sharedptr.h"

#include "ito33/cppunit.h"

#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/convertiblebond.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"


#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/version.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"


#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/testdata.h"

// local files
#include "bondlike_reader.h"
#include "basic_features_tester.h"
#include "bondlike_test_runners.h"

using namespace ito33;
using namespace ito33::finance;

using std::vector;
using std::string;

extern TagsForIhgBondLikeTest testTags;

extern BasicFeaturesTester basicFeatureTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !basicFeatureTester. ##name().Test(testTags)); \
}                                                                \


class IHGCBBasicTest : public CppUnit::TestCase
{
public:
  IHGCBBasicTest() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( IHGCBBasicTest );

  CPPUNIT_TEST_SUITE_END();
};

class CppUnitBasicFeaturesTest : public CppUnit::TestCase
{

public:
  CppUnitBasicFeaturesTest()  {}

  void tearDown() {
  }

 private:
  CPPUNIT_TEST_SUITE( CppUnitBasicFeaturesTest );

    CPPUNIT_TEST( TestGetBasicData );

    CPPUNIT_TEST( TestAddSoftCall);

    CPPUNIT_TEST( TestAddPVCouponMakeWhole);
    CPPUNIT_TEST( TestAddNoPVCouponMakeWhole);

    CPPUNIT_TEST( TestAddPremiumMakeWhole);
    
    CPPUNIT_TEST( TestAddCallNoticeToBasicData);

    CPPUNIT_TEST( TestAddCallNoticeToAddSoftCallData);

    CPPUNIT_TEST( TestPriceIncreasesAsYTCIncreases );

    CPPUNIT_TEST( TestPriceIncreasesAsYTPIncreases );

    CPPUNIT_TEST( TestPriceIncreasesAsCallNoticePeriodIncreases );

    CPPUNIT_TEST( TestAddCallClaimTriggerAsPercentageOfToBasicData );

    CPPUNIT_TEST( TestAddCallIssuePriceTriggerAsPercentageOfToBasicData );

    CPPUNIT_TEST( TestAddCallPrincipalTriggerAsPercentageOfToBasicData );

  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );

  IMPLEMENT_TEST( AddSoftCall);

  IMPLEMENT_TEST( AddPVCouponMakeWhole);
    
  IMPLEMENT_TEST( AddNoPVCouponMakeWhole);
  
  IMPLEMENT_TEST( AddPremiumMakeWhole);
  
  IMPLEMENT_TEST( AddCallNoticeToBasicData);

  IMPLEMENT_TEST( AddCallNoticeToAddSoftCallData);

  IMPLEMENT_TEST( AddCallClaimTriggerAsPercentageOfToBasicData );

  IMPLEMENT_TEST( AddCallIssuePriceTriggerAsPercentageOfToBasicData );

  IMPLEMENT_TEST( AddCallPrincipalTriggerAsPercentageOfToBasicData );


  void TestPriceIncreasesAsYTCIncreases()
  {
    basicFeatureTester.TestPriceIncreasesAsYTCIncreases();
  }
  

  void TestPriceIncreasesAsYTPIncreases()
  {
    basicFeatureTester.TestPriceIncreasesAsYTPIncreases();
  }

  void TestPriceIncreasesAsCallNoticePeriodIncreases()
  {
    basicFeatureTester.TestPriceIncreasesAsCallNoticePeriodIncreases();
  }
};

extern ito33::RegressionOutputSet<finance::ModelOutput> g_outputsCB;

namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

int RunBasicFeatureTester(bool bIsAcceptance)
{
  const size_t nNbBasicTests = 13;
  std::vector<size_t> pnIdxBasicTest;

  if(bIsAcceptance)
  {
    pnIdxBasicTest.resize(0);
    pnIdxBasicTest.push_back(1);
  }
  else
  {
    pnIdxBasicTest.resize(nNbBasicTests);
    for(size_t n = 0; n < nNbBasicTests; n++)
      pnIdxBasicTest[n] = n + 1;
  }

  try
  {    
    std::cout << "\nRunning Basic feature tester...\n";

    g_outputsCB.Read("xmlfiles/output_basic.xml");

    size_t n;
    for(n = 0; n < pnIdxBasicTest.size(); n++)
    {
      std::cout << "Testing : " << 
        String::Printf("xmlfiles/basic%02d.xml", pnIdxBasicTest[n]).c_str()
        << std::endl;

      basicFeatureTester.Setup(
        String::Printf("xmlfiles/basic%02d.xml", pnIdxBasicTest[n]).c_str());
  //    t.Test();
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitBasicFeaturesTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_basic_new.xml",
                            XML_IHGTEST_BASIC_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_basic.xml",
                              XML_IHGTEST_BASIC_OUTPUT,
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

