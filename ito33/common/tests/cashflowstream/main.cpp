/////////////////////////////////////////////////////////////////////////////
// Name:        tests/cashflowstream/main.cpp
// Purpose:     main file of cashflowstream test
// Author:      Zhang
// Created:     24.06.04
// RCS-ID:      $Id: main.cpp,v 1.3 2005/04/22 13:18:20 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testcashflowstream.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(CashFlowStreamTest::suite());

  return runner.run("") ? 0 : 1;
}

