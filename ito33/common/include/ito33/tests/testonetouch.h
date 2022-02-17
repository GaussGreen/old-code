/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/tests/testonetouch.h
// Purpose:     acceptance test for one touch
// Created:     2005/12/08
// RCS-ID:      $Id: testonetouch.h,v 1.1 2005/12/14 13:27:23 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class OneTouchTest : public CppUnit::TestCase
{

public:
  OneTouchTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( OneTouchTest );
    CPPUNIT_TEST( Delta2Barrier );
  CPPUNIT_TEST_SUITE_END();

  void Delta2Barrier();

  NO_COPY_CLASS(OneTouchTest);
};
