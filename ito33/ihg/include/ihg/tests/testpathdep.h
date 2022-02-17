/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testpathdep.h
// Purpose:     test for path dep options
// Author:      ITO 33
// Created:     04/02/2005
// RCS-ID:      $Id: testpathdep.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/common.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"

namespace ito33
{

namespace finance
{
  class ConvertibleBond;
  class PutSchedule;
  class ConversionSchedule;
  class CallSchedule;
}

}

// ----------------------------------------------------------------------------
// Acceptance tests for PathDependent
// ----------------------------------------------------------------------------
class PathDepTest : public CppUnit::TestCase
{
public:
  PathDepTest() { }

  virtual void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( PathDepTest );
   CPPUNIT_TEST( TestPathDepNoEvents );
   CPPUNIT_TEST( TestPathDepConstraints );
   CPPUNIT_TEST( TestContinuousEvent );
   CPPUNIT_TEST( TestDiscreteEvent ) ;
  CPPUNIT_TEST_SUITE_END();


  // In this test we check that the path dependent
  // solver gives the same answer as in the regular
  // solver when no events are created
  // solve xx problems
  void TestPathDepNoEvents();
  
   // In this tests we check that the path dependent
  // constraints are created correctly and that 
  // they are not issues with the different
  // contstraints
  //  We solve a cb using the classic method and compare
  //  the solution with the pathDependent method
  //  We consider the following scenarios
  //  0. no calls, no conversion, no puts     PDE0
  //  1. calls, no conversion, no puts        PDE1
  //  2. no calls, conversions, no puts       PDE2
  //  3. no calls, no conversion, put         PDE3
  //  4. calls, conversions, puts             PDE4
  //  5. calls notice, no conversion, no puts PDE5
  void TestPathDepConstraints();


    //
  // continous call/conversion pde PDE solution  
  // is compared with  two path
  // dependent pdes
  // one with call/conversion provision PDE1, 
  // the other one without call/conversion provision PDE0
  //
  // at each timestep pde1 is copied inside pde0
  // at the end pde0 should be equal to the regular
  // method of solving the pde
  //
  void TestContinuousEvent();
  
  
  //
  // test discrete events are handled
  // properly
  //    5 discrete event dates
  //     calls/put/conversion at each of these dates
  void TestDiscreteEvent();

  //util function
  double Price(ito33::shared_ptr<ito33::finance::ConvertibleBond> &pCB,
    double dVol,double dLambda,ito33::Date issueDate);

  //util function
  void SetConstraints(size_t nPDENum,
    ito33::shared_ptr<ito33::finance::CallSchedule> &pCall,
    ito33::shared_ptr<ito33::finance::PutSchedule>  &pPut,
    ito33::shared_ptr<ito33::finance::ConversionSchedule> &pConv,
    ito33::Date startCallDate, ito33::Date endCallDate,
    ito33::Date startConvDate, ito33::Date endConvDate,
    ito33::Date putDate);

  //util function
  ito33::shared_ptr<ito33::finance::ConvertibleBond> 
    CreateCB(size_t nPDENum,ito33::Date issueDate, ito33::Date maturityDate,
             ito33::Date startCallDate, ito33::Date endCallDate,
             ito33::Date startConvDate, ito33::Date endConvDate,
             ito33::Date putDate);

  NO_COPY_CLASS(PathDepTest);
};
