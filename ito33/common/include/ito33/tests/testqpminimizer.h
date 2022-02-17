/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/tests/testqpminimizer.h
// Purpose:     unit tests for non linear equation 1D solvers
// Author:      Pedro Ferreira, Laurence Gozalo
// Created:     2005/06/22
// RCS-ID:      $Id: testqpminimizer.h,v 1.1 2005/06/27 13:13:46 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TESTS_QPMINIMIZER_H_
#define _ITO33_TESTS_QPMINIMIZER_H_

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"

// Test class
class QPMinimizerTest : public CppUnit::TestCase {

public:

  QPMinimizerTest() { }


private:

  CPPUNIT_TEST_SUITE( QPMinimizerTest );

    // With an identity matrix, we should get the negative of the c vector
    CPPUNIT_TEST ( Identity );
    
    // Test the storage, see if the way it's used in hedging has no problem
    // A non degenerated matrix is used
    CPPUNIT_TEST ( Appended );

  CPPUNIT_TEST_SUITE_END();

  void Identity();
  void Appended();

  NO_COPY_CLASS( QPMinimizerTest );
};

#endif // #ifndef _ITO33_TESTS_QPMINIMIZER_H_
