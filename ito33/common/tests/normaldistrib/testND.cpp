/////////////////////////////////////////////////////////////////////////////
// Name:        test/normaldistrib/testND.cpp
// Purpose:     test file for normal distribution function
// Author:      Laurence
// Created:     22/09/2003
// RCS-ID:      $Id: testND.cpp,v 1.7 2004/10/05 09:13:51 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_ND_H_
#define _ITO33_TEST_ND_H_
// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <cmath>

#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/numeric/normaldistribution.h"
#include "ito33/cppunit.h"

#include "ito33/tests/testND.h"

using namespace ito33;
using namespace ito33::numeric;

// Global test function
//double polynom(double x)
//{
//  return (x-7.0)*x+12.0;
//}


void NormalDistTest::normalof0()
{
  // Verifies if the value of normal(0) is \f$ \frac{1}{\sqrt{2 \pi}} $\f
  CPPUNIT_ASSERT_DOUBLES_EQUAL( Normal(0.) , INVSQRT2PI , 1e-7 );
}

void NormalDistTest::normalofd1()
{
CPPUNIT_ASSERT_DOUBLES_EQUAL( Normal(-0.0500033335333477) , 0.398443847681252 , 1e-7 );
}
void NormalDistTest::normalof1_96()
{
  // Verifies if the value of normal(1,96) corresponds to the one given in 
  // Abramovitz/Stegun p966
  CPPUNIT_ASSERT_DOUBLES_EQUAL( Normal(1.96) , 0.0584409443 , 1e-7 );
}

void NormalDistTest::cumulatednormalof0()
{
  // Verifies if the value of cumulatednormal(0) is 0.5
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CumulatedNormal(0.) , 0.5 , 1e-7 );
}

void NormalDistTest::cumulatednormalof1_96()
{
  // Verifies if the value of cumulatednormal(1.96) corresponds to the one  
  // given in Abramovitz/Stegun p968
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CumulatedNormal(1.96) , 0.97500210485178 , 1e-7 );
}


void NormalDistTest::cumulnormsymetry()
{
  // Verifies the symetry of the cumulatednormal : we should have
  // P(-x)=1-P(x), where P(x) is the cumulated normal of x, for x > 0.
  // Here we choose x=2.
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CumulatedNormal(-2.) + CumulatedNormal(2.) , 
    1. , 1e-12 );
}

#endif
