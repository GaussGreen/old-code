/////////////////////////////////////////////////////////////////////////////
// Name:        tests/surfaces/main.cpp
// Purpose:     main file for testing the time mesh generation
// Author:      David
// Created:     25.05.04
// RCS-ID:      $Id: main.cpp,v 1.3 2004/10/05 09:13:52 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testtimemesh.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(TimeMeshTest::suite());

  return runner.run("") ? 0 : 1;
}

