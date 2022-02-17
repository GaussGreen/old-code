/////////////////////////////////////////////////////////////////////////////
// Name:        tests/simulatedannealing/main.cpp
// Purpose:     main file for unit tests for Simulated Annealing
// Author:      ITO 33
// Created:     April 18, 2005
// RCS-ID:      $Id: main.cpp,v 1.2 2005/05/25 14:53:51 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testsa.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(SATest::suite());
  runner.addTest(ASATest::suite());

  return runner.run("") ? 0 : 1;
}
