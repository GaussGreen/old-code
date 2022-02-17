/////////////////////////////////////////////////////////////////////////////
// Name:        common/tests/calibrator/main.cpp
// Purpose:     main file for unit tests for qp minimizer
// Created:     2005/06/22
// RCS-ID:      $Id: main.cpp,v 1.1 2005/06/23 08:29:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testqpminimizer.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(QPMinimizerTest::suite());

  return runner.run("") ? 0 : 1;
}
