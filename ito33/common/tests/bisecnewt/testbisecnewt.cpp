#include "ito33/beforestd.h"
#include <cmath>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/numeric/bisecnewton.h"
#include "ito33/cppunit.h"

#include "ito33/tests/testbisecnewt.h"

using namespace std;
using namespace ito33;

// Global test function
static void polynom(double x, double &f, double &df)
{
  f = (x - 7.0) * x + 12.0;
  df = 2. * x - 7.;
}

static void mysin(double x, double &f, double &df)
{
  f = sin(x);
  df = cos(x);
}

void BisecNewtTest::Sin()
{
  numeric::BisecNewton solver;

  // Giving a small initial guess we shoul find the root zero
  CPPUNIT_ASSERT_DOUBLES_EQUAL( solver(mysin, -0.1, 0.1), 0.0, 1e-5 );
}

void BisecNewtTest::Polynom()
{
  numeric::BisecNewton solver;

  // The roots of the polynom are 3 and 4 so we must find one of them
  double dResult = solver(polynom,2.,3.);
  CPPUNIT_ASSERT( fabs(dResult-3.0) < 1e-5 || fabs(dResult-4.0) < 1e-5 );
}

void BisecNewtTest::NegativeTolerance()
{
  // Check that tolerance cannot be negative
  numeric::BisecNewton solver(-1.e-4);
}

void BisecNewtTest::InitEqual()
{
  // Test that the bracket interval is only a single point
  numeric::BisecNewton solver;
  
  solver(polynom, 0.5,0.5);
}

void BisecNewtTest::MaxIterations()
{
  numeric::BisecNewton solver(1e-15, 1);

  // Test the max iterations
  solver(mysin, -0.1, 0.15);
}
