/////////////////////////////////////////////////////////////////////////////
// Name:        test/error/main.cpp
// Purpose:     main file of error example/test program
// Author:      Vadim Zeitlin
// Created:     12.03.04
// RCS-ID:      $Id: main.cpp,v 1.3 2004/10/05 09:13:50 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testerror.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(ErrorCodeTestCase::suite());

  return runner.run("") ? 0 : 1;
}

