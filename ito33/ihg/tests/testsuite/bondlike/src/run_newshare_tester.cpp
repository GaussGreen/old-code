/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/bondlike/src/run_newshare_tester.cpp
// Purpose:     Tests for new share feature
// Author:      Nabil
// Created:     2005/04/13
// RCS-ID:      $Id: run_newshare_tester.cpp,v 1.4 2006/05/09 15:02:30 nabil Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "newshare_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern NewShareTester newshareTester;

#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !newshareTester. ##name().Test(testTags));   \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitNewShareTest : public CppUnit::TestCase
{

public:
  
  CppUnitNewShareTest() {}

  void setUp() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( CppUnitNewShareTest );
    CPPUNIT_TEST( TestGetBasicData );
    CPPUNIT_TEST( PriceDecreaseWhenNewShare );
  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );

  void PriceDecreaseWhenNewShare()
  {
    newshareTester.PriceDecreaseWhenNewShare();
  }
  
};

int RunNewShareTester(bool bIsAcceptance)
{
  const size_t nNbNewShareCBTests = 7;
  std::vector<size_t> pnIdxNewShareCBTest;

  if(bIsAcceptance)
  {
    pnIdxNewShareCBTest.resize(0);
    pnIdxNewShareCBTest.push_back(1);
  }
  else
  {
    size_t n;
    pnIdxNewShareCBTest.resize(nNbNewShareCBTests);
    for(n = 0; n < nNbNewShareCBTests; n++)
      pnIdxNewShareCBTest[n] = n + 1;
  }
  
  try
  {
    size_t n;
    std::cout << "\nRunning tests for new share...\n";
    g_outputsCB.Read("xmlfiles/output_newshare.xml");

    ///------------------- cb new share test ---------------------  
    for(n = 0; n < pnIdxNewShareCBTest.size(); n++)
    {
      newshareTester.Setup( String::Printf("xmlfiles/cbnewshare%02d.xml",
                            pnIdxNewShareCBTest[n]).c_str() );
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitNewShareTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_newshare_new.xml",
                            XML_IHGTEST_NEWSHARE_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_newshare.xml",
                              XML_IHGTEST_NEWSHARE_OUTPUT,
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
