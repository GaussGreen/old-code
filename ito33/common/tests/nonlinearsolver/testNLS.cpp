/////////////////////////////////////////////////////////////////////////////
// Name:        tests/nonlinearsolver/testNLS.cpp
// Purpose:     unit tests for non linear equation 1D solvers
// Author:      Pedro Ferreira, Laurence Gozalo
// Created:     04.12.03
// RCS-ID:      $Id: testNLS.cpp,v 1.14 2006/03/01 14:43:36 yann Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/numeric/nonlinearsolver.h"
#include "ito33/cppunit.h"

#include "ito33/tests/testNLS.h"

#include <cmath>
#include <iostream>

using namespace std;
using namespace ito33;

using ito33::numeric::RegulaFalsi;
using ito33::numeric::Secant;

static double linear(double x)
{
  return -3 * x + 4;
}
// Global test functions
static double polynom(double x)
{
  return ( x - 7.0) * x + 12.0;
}

// this function is monotonuous and so can be used with
// RegulaFalsi::FindMatchingRootBracket()
static double polynom3(double x)
{
  return ((x - 1.) * x + 1.) * x - 1.;
}

static double polynom0(double /* x */)
{
  return 17.;
}

static double expfunc(double x)
{
  return exp(x) - 3;
}

void NLSTest::Sin()
{
  Secant solver;

  // Giving a small initial guess we shoul find the root zero
  CPPUNIT_ASSERT_DOUBLES_EQUAL( solver(sin, 0.1, 0.15), 0.0, 1e-5);
}


void NLSTest::IterPol()
{
  Secant solver;

  solver(polynom, 0.1, 0.5);
  
  solver.GetNbIter();

  CPPUNIT_ASSERT(solver.GetNbIter() > 0);
}

void NLSTest::Polynom()
{
  Secant solver;

  // The roots of the polynom are 3 and 4 so we must find one of them
  double dResult = solver(polynom,0.1,0.5);
  CPPUNIT_ASSERT( fabs(dResult-3.0) < 1e-5 || fabs(dResult-4.0) < 1e-5 );
}


void NLSTest::RFSin()
{
  RegulaFalsi solver;

  // Giving a small initial guess we shoul find the root zero
  CPPUNIT_ASSERT_DOUBLES_EQUAL( solver(sin, -0.1, 0.15), 0.0, 1e-5 );
}


void NLSTest::RFIterPol()
{
  RegulaFalsi solver;

  solver(polynom, 2.5, 3.2);

  CPPUNIT_ASSERT(solver.GetNbIter() > 0);
}

void NLSTest::RFPolynom()
{
  RegulaFalsi solver;
  
  // The roots of the polynom are 3 and 4 so we must find one of them
  double dResult = solver(polynom,2.5,3.2);
  CPPUNIT_ASSERT( fabs(dResult-3.0) < 1e-5 || fabs(dResult-4.0) < 1e-5 );
}

void NLSTest::RFFindMatch()
{
  double dValue2 = RegulaFalsi::FindMatchingRootBracket(polynom3, 0.);
  CPPUNIT_ASSERT( polynom3(dValue2) > 0. );

  dValue2 = RegulaFalsi::FindMatchingRootBracket(polynom3, 1.);

  RegulaFalsi solver;

  double dResult = solver(polynom3, 1., dValue2);

  ITO33_ASSERT_DOUBLES_EQUAL( dResult, 1. );

  dValue2 = RegulaFalsi::FindMatchingRootBracket(polynom3, 2.);
  CPPUNIT_ASSERT( polynom3(dValue2) < 0. );
}

void NLSTest::RFFindMatchBad()
{
  RegulaFalsi::FindMatchingRootBracket(polynom0, 0.);
}

void NLSTest::RFPolynom3()
{
  RegulaFalsi solver;

  // The roots of the polynom is 1 so we must find it
  double dResult = solver(polynom3,0.,10);
  CPPUNIT_ASSERT( fabs(dResult - 1.0) < 1e-5 );
}

void NLSTest::RFLinearPolynom()
{
  RegulaFalsi solver;

  double dResult = solver(linear, 0. , 2.);

 
  CPPUNIT_ASSERT( fabs( dResult - 4./3.) < 1.e-16);
}


void NLSTest::RFExp()
{
  RegulaFalsi solver;

  double dResult = solver(expfunc, 0. , 2.);

 
  CPPUNIT_ASSERT( fabs( dResult - log(3.) ) < 1.e-5);
}

void NLSTest::Tolerance()
{
  // Check that tolerance cannot be negative
  Secant solver(-0.01);
}

void NLSTest::Div0()
{
  Secant solver(1e-5);

  // Test divide by zero
  solver(polynom,0.1, 0.1);
}

void NLSTest::MaxIterations()
{
  Secant solver(1e-15, 3);

  // Test the max iterations
  solver(sin, 0.2, 0.3);
}
