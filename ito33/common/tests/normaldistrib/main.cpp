/////////////////////////////////////////////////////////////////////////////
// Name:        test/normaldistrib/testND.cpp
// Purpose:     main test file for normal distribution function
// Author:      Laurence
// Created:     22/09/2003
// RCS-ID:      $Id: main.cpp,v 1.2 2004/10/05 09:13:51 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testND.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(NormalDistTest::suite());

  return runner.run("") ? 0 : 1;
}
