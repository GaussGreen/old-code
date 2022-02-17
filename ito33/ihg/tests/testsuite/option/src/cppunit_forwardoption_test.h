
/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/cppunit_option_test.h
// Purpose:     Base class for testing IHG projects
// Author:      ITO33 Canada
// Created:     2005/06/08
// RCS-ID:      $Id: cppunit_forwardoption_test.h,v 1.2 2005/06/27 16:20:57 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/cppunit_option_test.h
**/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_FORWARDOPTION_TEST_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_FORWARDOPTION_TEST_H_

#include "ito33/cppunit.h"

namespace ito33 
{
namespace ihg
{
namespace test
{

class CppUnitForwardOptionTest: public CppUnit::TestCase
{
 
private:

public:
  CppUnitForwardOptionTest() {}

  void setup();


 private:
  CPPUNIT_TEST_SUITE( CppUnitForwardOptionTest );
    CPPUNIT_TEST( Test1 );
    CPPUNIT_TEST( Test2 );
    CPPUNIT_TEST( Test3 );
    CPPUNIT_TEST( Test4 );
    CPPUNIT_TEST( Test5 );
    CPPUNIT_TEST( Test6 );
    CPPUNIT_TEST( Test7 );
    CPPUNIT_TEST( Test8 );
  CPPUNIT_TEST_SUITE_END();

  void Test1();
  void Test2();
  void Test3();
  void Test4();
  void Test5();
  void Test6();
  void Test7();
  void Test8();

  

  NO_COPY_CLASS(CppUnitForwardOptionTest);

};

} //end namespace test
} //end namespace ihg
}//end namespace ito33

#endif //_IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_FORWARDOPTION_TEST_H_