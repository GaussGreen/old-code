/////////////////////////////////////////////////////////////////////////////
// Name:        test/dividends/main.cpp
// Purpose:     Unit test for Dividends class
// Author:      Vadim Zeitlin
// Created:     30.07.03
// RCS-ID:      $Id: main.cpp,v 1.9 2004/10/05 09:13:50 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testdividends.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(DivsTestCase::suite());

  return runner.run("") ? 0 : 1;
}

