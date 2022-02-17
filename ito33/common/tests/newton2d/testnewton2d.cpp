/////////////////////////////////////////////////////////////////////////////
// Name:        tests/newton2d/testnewton2d.cpp
// Purpose:     unit tests for non linear Newton 2d solvers
// Author:      ITO 33
// Created:     2002/12/22
// RCS-ID:      $Id: testnewton2d.cpp,v 1.3 2005/02/07 12:47:15 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"

#include "ito33/numeric/newton2d.h"
#include "ito33/numeric/numericerror.h"

#include "ito33/tests/testnewton2d.h"

using namespace ito33;
using namespace ito33::numeric;


// Global test functions
void linear(double x, double y, double& dFunc1, double& dFunc2,
            double& dF)
{
  // equation 1: f1(x,y) = x - y
  // equation 2: f2(x,y) = 2x + 3y - 2
  // solution at (0.4, 0.4)
  dFunc1 = x - y;
  dFunc2 = 2.0*x + 3.0*y - 2.0;

  dF = 0.5 * (dFunc1 * dFunc1 + dFunc2 * dFunc2);
}

void quadratic(double x, double y, double& dFunc1, double& dFunc2,
            double& dF)
{
  // equation 1: f1(x,y) = x*x - y*y
  // equation 2: f2(x,y) = 2x + 3y*y - 4
  // solution at (0.4, 0.4)
  dFunc1 = x*x - y*y;
  dFunc2 = 2.0*x + 3.0*y*y - 4.0;

  dF = 0.5 * (dFunc1 * dFunc1 + dFunc2 * dFunc2);
}


void Newton2DTest::operator () 
   (double dX, double dY, double &dFunc1, double &dFunc2,
    double &dF)
{

  switch(test){
  case LINEAR:
    linear(dX, dY, dFunc1, dFunc2, dF);
    break;
  case QUADRATIC:
    quadratic(dX, dY, dFunc1, dFunc2, dF);
    break;
  default:
    CPPUNIT_ASSERT(false);
  }

}



void Newton2DTest::Linear()
{

  test = LINEAR;

  double dParam1Lower = -10.0;
  double dParam1Upper =  10.0;
  double dParam2Lower = -10.0;
  double dParam2Upper =  10.0;

  // linear problem only needs 1 iteration
  Newton2D solver(dParam1Lower, dParam1Upper, dParam2Lower, dParam2Upper, 1.e-6, 1);

  double dX = 3.0;
  double dY = -2.1;
  numeric::NumericError err = solver(*this, dX, dY);

  // exact answer is (0.4, 0.4)
  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, dX, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, dY, 1e-6);
}


void Newton2DTest::Quadratic()
{

  test = QUADRATIC;

  double dParam1Lower = -10.0;
  double dParam1Upper =  10.0;
  double dParam2Lower = -10.0;
  double dParam2Upper =  10.0;

  // quadrartic problem should converge
  Newton2D solver(dParam1Lower, dParam1Upper, dParam2Lower, dParam2Upper, 1.e-14, 15);

  double dX = 1.2;
  double dY = 1.1;
  numeric::NumericError err = solver(*this, dX, dY);
  
  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.868517091, dX, 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.868517091, dY, 1e-6);
}
