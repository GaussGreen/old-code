/////////////////////////////////////////////////////////////////////////////
// Name:        tests/hazardrates/main.cpp
// Purpose:     main file for testing hazard rates
// Author:      David Pooley
// Created:     16/06/2004
// RCS-ID:      $Id: main.cpp,v 1.3 2005/10/27 22:27:24 dave Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/numeric/mesh/specialtimes.h"

#include "ihg/tests/testhazardrates.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(HazardRateFlatTest::suite());
  runner.addTest(HazardRateTimeOnlyTest::suite());
  runner.addTest(HazardRatePowerTest::suite());
  runner.addTest(HazardRateComboTest::suite());
  runner.addTest(HRSpotComponentPowerTest::suite());

  return runner.run("") ? 0 : 1;
}

