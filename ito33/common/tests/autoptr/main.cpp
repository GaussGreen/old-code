/////////////////////////////////////////////////////////////////////////////
// Name:        tests/autoptr/main.cpp
// Purpose:     main file of AutoPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: main.cpp,v 1.4 2004/10/05 09:13:48 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testautoptr.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(AutoPtrTestCase::suite());

  return runner.run("") ? 0 : 1;
}

