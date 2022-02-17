/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testacceptance/main.cpp
// Purpose:     cppunit testing for ihg functionnalitilites
// Author:      Yann d'Halluin
// Created:     11/07/2004
// RCS-ID:      $Id: main.cpp,v 1.5 2005/12/28 09:38:00 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ihg/tests/testhazardrates.h"
#include "ihg/tests/testvolatility.h"
#include "ihg/tests/testpathdep.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCB);
ITO33_FORCE_LINK_MODULE(IHGPriceCDS);

int main(int ,char **)
{

 
  /**************** UNIT TESTING FOR IHG ******************/
  CppUnit::TextUi::TestRunner runnerUnitTest;
  
  runnerUnitTest.addTest(HazardRateFlatTest::suite());
  runnerUnitTest.addTest(HazardRateTimeOnlyTest::suite());
  runnerUnitTest.addTest(HazardRatePowerTest::suite());
  runnerUnitTest.addTest(HazardRateComboTest::suite());
  runnerUnitTest.addTest(HRSpotComponentPowerTest::suite());

  runnerUnitTest.addTest(VolatilityFlatTest::suite());
  runnerUnitTest.addTest(VolatilityPowerTest::suite());
  runnerUnitTest.addTest(VolatilityTanhTest::suite());
  runnerUnitTest.addTest(VolatilityTimeOnlyTest::suite());
	 
  runnerUnitTest.addTest( PathDepTest::suite() );

  return runnerUnitTest.run("") ? 0 : 1;;
}

