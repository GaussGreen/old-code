/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testputcallparity.h
// Purpose:     test for put call parity helper functions
// Author:      ITO 33
// Created:     13/05/2005
// RCS-ID:      $Id: testputcallparity.h,v 1.5 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/sharedptr.h"

namespace ito33
{

namespace finance
{
  class Option;
  class SessionData;
}

namespace ihg
{
  class TheoreticalModel;
}

}


// ----------------------------------------------------------------------------
// Acceptance tests for put call parity functions
// ----------------------------------------------------------------------------
class PutCallParityTest : public CppUnit::TestCase
{
public:
  PutCallParityTest() { }

  virtual void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( PutCallParityTest );
   CPPUNIT_TEST( TestConvertPutToCall1 );
   CPPUNIT_TEST( TestConvertPutToCall2 );
   CPPUNIT_TEST( TestConvertPutToCall3 );
  CPPUNIT_TEST_SUITE_END();


  // hazard rate no dividends
  void TestConvertPutToCall1();

  // dividends no hazard rate
  void TestConvertPutToCall2();

    // dividend and hazard rate
  void TestConvertPutToCall3();
  

  // Helper function: Construct a model for pricing
  ito33::shared_ptr<ito33::ihg::TheoreticalModel> ConstructModel();

  // Helper function: Construct a session for pricing
  ito33::shared_ptr<ito33::finance::SessionData> 
    ConstructSessionData(int iDivType);

  // Helper function: Construct a basic put option for pricing
  ito33::shared_ptr<ito33::finance::Option> 
    ConstructPutOption(ito33::shared_ptr<ito33::finance::SessionData> pSessionData,
                       ito33::shared_ptr<ito33::ihg::TheoreticalModel> pModel);


  NO_COPY_CLASS(PutCallParityTest);
};
