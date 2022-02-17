/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/bondlike/src/run_cboption_tester.cpp
// Purpose:     Tests for cb option instrument
// Author:      Nabil
// Created:     2005/09/20
// RCS-ID:      $Id: run_cboption_tester.cpp,v 1.7 2006/07/20 14:50:07 nabil Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/ihg/version.h"
// local files
#include "cboption_tester.h"

using namespace ito33;

extern TagsForIhgBondLikeTest testTags;

extern CBOptionTester cboptionTester;

#define IMPLEMENT_TEST(name)                                     \
  void Test ##name()                                             \
{                                                                \
  CPPUNIT_ASSERT( !cboptionTester. ##name().Test(testTags));   \
}                                                                \


namespace ito33
{
extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

extern size_t g_nNbFailedSelfRegressionTests;
}

class CppUnitCBOptionTest : public CppUnit::TestCase
{

public:
  
  CppUnitCBOptionTest() {}

  void setUp() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( CppUnitCBOptionTest );
    CPPUNIT_TEST( TestGetBasicData );
    CPPUNIT_TEST( PriceIncreaseWithRecallSpread );
  CPPUNIT_TEST_SUITE_END();

private:

  IMPLEMENT_TEST( GetBasicData );

  void PriceIncreaseWithRecallSpread()
  {
    cboptionTester.PriceIncreaseWithRecallSpread();
  }
  
};

int RunCBOptionTester(bool bIsAcceptance)
{
  const size_t nNbCBOptionTests = 7;
  std::vector<size_t> pnIdxCBOptionTest;

  if(bIsAcceptance)
  {
    pnIdxCBOptionTest.resize(0);
    pnIdxCBOptionTest.push_back(1);
  }
  else
  {
    size_t n;
    pnIdxCBOptionTest.resize(nNbCBOptionTests);
    for(n = 0; n < nNbCBOptionTests; n++)
      pnIdxCBOptionTest[n] = n + 1;
  }
  
  try
  {
    size_t n;
    std::cout << "\nRunning tests for cb option...\n";
    g_outputsCB.Read("xmlfiles/output_cboption.xml");

    ///------------------- cb option test ---------------------  
    for(n = 0; n < pnIdxCBOptionTest.size(); n++)
    {
      cboptionTester.Setup( String::Printf("xmlfiles/cboption%02d.xml",
                            pnIdxCBOptionTest[n]).c_str() );
      
      CppUnit::TextUi::TestRunner runner;

      runner.addTest( CppUnitCBOptionTest::suite() );
      runner.run("");
    }

    if(g_outputsCB.Changed())
      g_outputsCB.Write( "xmlfiles/output_cboption_new.xml",
                            XML_IHGTEST_CBOPTION_OUTPUT,
                            ITO33_XML_ATTR_ROOT_VERSION,
                            ITO33_IHG_VERSION_DOT_STRING);
    else
    {
      if(g_outputsCB.Added())
      { // replace old output file 
        g_outputsCB.Write( "xmlfiles/output_cboption.xml",
                              XML_IHGTEST_CBOPTION_OUTPUT,
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
