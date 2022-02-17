/////////////////////////////////////////////////////////////////////////////
// Name:        tests/newton2d/main.cpp
// Purpose:     main file for unit tests for nonlinear newton 2d solver
// Author:      ITo 33
// Created:     2004/12/22
// RCS-ID:      $Id: main.cpp,v 1.2 2005/01/28 19:34:05 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testnewton2d.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(Newton2DTest::suite());

  return runner.run("") ? 0 : 1;
}
