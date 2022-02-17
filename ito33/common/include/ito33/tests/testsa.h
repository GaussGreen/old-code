/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/tests/testsa.h
// Purpose:     unit tests for simulated annealing
// Author:      ITO 33
// Created:     April 18, 2005
// RCS-ID:      $Id: testsa.h,v 1.2 2005/05/25 14:29:08 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_SIMULATEDANNEALING_H_
#define _ITO33_TEST_SIMULATEDANNEALING_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"


// Test class
class SATest : public CppUnit::TestCase {
public:

  enum Test {LINEAR, ROSENBROCK, SPHERE, PLATEAU, CALIBRATION};

  SATest( ) {}

  void tearDown() {}

  /** 
    This is called by solver.  Each test needs to set
    what this operator will return.
  */
  void operator () (const std::vector<double> &pdParam, double &dF);

private:

  CPPUNIT_TEST_SUITE( SATest );
    CPPUNIT_TEST ( TestLinear ); 
    CPPUNIT_TEST ( TestRosenbrock );
    CPPUNIT_TEST ( TestSphere );
    CPPUNIT_TEST( TestPlateau );
    //CPPUNIT_TEST( TestOptionCalibration );
  CPPUNIT_TEST_SUITE_END();

  void TestLinear();
  void TestRosenbrock();
  void TestSphere();
  void TestPlateau();
  void TestOptionCalibration();

  Test m_test; //indicate which test to run

  void CheckObjectiveFunction(const std::vector<double> &pdParam);

  NO_COPY_CLASS( SATest );
};


// Test class
class ASATest : public CppUnit::TestCase {
public:

  enum Test {LINEAR, ROSENBROCK, SPHERE, PLATEAU, CALIBRATION};

  ASATest( ) {}


  /** 
    This is called by solver.  Each test needs to set
    what this operator will return.
  */
  void operator () (const std::vector<double> &pdParam, double &dF);

private:

  CPPUNIT_TEST_SUITE( ASATest );
    CPPUNIT_TEST ( TestRosenbrock );
    CPPUNIT_TEST ( TestSphere );
    CPPUNIT_TEST( TestPlateau );
    //CPPUNIT_TEST( TestOptionCalibration );
  CPPUNIT_TEST_SUITE_END();

  void TestRosenbrock();
  void TestSphere();
  void TestPlateau();
  void TestOptionCalibration();

  Test m_test; //indicate which test to run

  void CheckObjectiveFunction(const std::vector<double> &pdParam);

  NO_COPY_CLASS( ASATest );
};

#endif
