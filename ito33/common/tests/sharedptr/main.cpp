/////////////////////////////////////////////////////////////////////////////
// Name:        test/sharedptr/main.cpp
// Purpose:     main file of SharedPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: main.cpp,v 1.4 2004/10/05 09:13:52 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/tests/testsharedptr.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(SharedPtrTestCase::suite());

  return runner.run("") ? 0 : 1;
}
