/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/tests/testnewton2d.h
// Purpose:     unit tests for newton 2d non linear solver
// Author:      ITO 33
// Created:     2004/12/22
// RCS-ID:      $Id: testnewton2d.h,v 1.2 2005/01/28 19:34:04 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_NEWTON2D_H_
#define _ITO33_TEST_NEWTON2D_H_

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"


// Test class
class Newton2DTest : public CppUnit::TestCase {
public:

  enum Test
  {
    LINEAR,
    QUADRATIC
  };

  Newton2DTest( ) {}

  void tearDown() {}

  /** 
    This is called by the newton 2d solver.  Each test needs to set
    what this operator will return.
  */
  void operator () (double dParam1, double dParam2, 
                    double &dFunc1,double &dFunc2,
                    double &dF);

private:

  CPPUNIT_TEST_SUITE( Newton2DTest );
    CPPUNIT_TEST ( Linear );
    CPPUNIT_TEST ( Quadratic );
  CPPUNIT_TEST_SUITE_END();

  void Linear();
  void Quadratic();

  /// which test functions are being used
  Test test;

  NO_COPY_CLASS( Newton2DTest );
};

#endif
