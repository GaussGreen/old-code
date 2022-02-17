/////////////////////////////////////////////////////////////////////////////
// Name:        tests/spotfxrates/main.cpp
// Purpose:     main file of SpotFXRates test program
// Author:      Wang
// Created:     2004/09/02
// RCS-ID:      $Id: main.cpp,v 1.2 2004/10/05 09:13:52 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/tests/testspotfxrates.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(SpotFXRatesTest::suite());

  return runner.run("") ? 0 : 1;
}
