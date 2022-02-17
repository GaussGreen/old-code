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

#include "ito33/finance/bondlike/bond.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"


#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/bondlikeoutput.h"
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
#include "bond_tester.h"
#include "bondlike_test_runners.h"

using namespace ito33;
using namespace ito33::finance;

using std::vector;
using std::string;

extern TagsForIhgBondLikeTest testTags;

extern BondTester bondTester;


#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !bondTester. ##name().Test(testTags)); \
}                                                                \

namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitBondTest : public CppUnit::TestCase
{

public:
  CppUnitBondTest()  {}

  void tearDown() {
  }

 private:
  CPPUNIT_TEST_SUITE( CppUnitBondTest );

    CPPUNIT_TEST( TestGetBasicData );

    CPPUNIT_TEST( TestPriceIncreasesAsYTCIncreases );

    CPPUNIT_TEST( TestPriceIncreasesAsYTPIncreases );

    CPPUNIT_TEST( TestPriceWithCouponsDefinedByFDF );

    CPPUNIT_TEST( TestFugitForBondWithNoConstraintsAndNoDefault );

    CPPUNIT_TEST( TestFugitForBondWithNoConstraintsAndNonNullFlatHR );

  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );

  void TestPriceIncreasesAsYTCIncreases()
  {
    bondTester.TestPriceIncreasesAsYTCIncreases();
  }
  
  void TestPriceIncreasesAsYTPIncreases()
  {
    bondTester.TestPriceIncreasesAsYTPIncreases();
  }
  
  void TestPriceWithCouponsDefinedByFDF()
  {
    bondTester.TestPriceWithCouponsDefinedByFDF();
  }
  
  void TestFugitForBondWithNoConstraintsAndNoDefault()
  {
    bondTester.TestFugitForBondWithNoConstraintsAndNoDefault();
  }
  
  void TestFugitForBondWithNoConstraintsAndNonNullFlatHR()
  {
    bondTester.TestFugitForBondWithNoConstraintsAndNonNullFlatHR();
  }

};

int RunBondTester(bool bIsAcceptance)
{
  const size_t nNbBondTests = 2;
  std::vector<size_t> pnIdxBondTest;

  if(bIsAcceptance)
  {
    pnIdxBondTest.resize(0);
    pnIdxBondTest.push_back(1);
  }
  else
  {
    pnIdxBondTest.resize(nNbBondTests);
    for(size_t n = 0; n < nNbBondTests; n++)
      pnIdxBondTest[n] = n + 1;
  }

  try
  {    
    std::cout << "\nRunning Bond tester...\n";

    g_outputsCB.Read("xmlfiles/output_bond.xml");

    size_t n;
    for(n = 0; n < pnIdxBondTest.size(); n++)
    {
      std::cout << "Testing : " << 
        String::Printf("xmlfiles/bond%02d.xml", pnIdxBondTest[n]).c_str()
        << std::endl;

      bondTester.Setup(
        String::Printf("xmlfiles/bond%02d.xml", pnIdxBondTest[n]).c_str());
  //    t.Test();
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitBondTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_bond_new.xml",
                            XML_IHGTEST_BOND_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_bond.xml",
                              XML_IHGTEST_BOND_OUTPUT,
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

