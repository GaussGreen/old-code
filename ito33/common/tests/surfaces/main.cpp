/////////////////////////////////////////////////////////////////////////////
// Name:        tests/surfaces/main.cpp
// Purpose:     main file of domain/surface test programs
// Author:      David Pooley
// Created:     17.05.04
// RCS-ID:      $Id: main.cpp,v 1.6 2004/10/14 16:03:06 zhang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testsurfaces.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(DomainFixedSpaceMeshTest::suite());
  runner.addTest(SurfaceGeneralTest::suite());
  runner.addTest(DomainGeneralWithSurfaceTest::suite());

  return runner.run("") ? 0 : 1;
}
