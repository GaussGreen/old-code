/////////////////////////////////////////////////////////////////////////////
// Name:        tests/floatingrates/main.cpp
// Purpose:     main file of floating rates test
// Author:      Nabil Ouachani
// Created:     2005/01/10
// RCS-ID:      $Id: main.cpp,v 1.2 2005/09/27 16:32:20 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testfloatingrates.h"

using namespace ito33;

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(FloatingRatesTest::suite());

  return runner.run() ? 0 : 1;
}

