/////////////////////////////////////////////////////////////////////////////
// Name:        test/datadump/main.cpp
// Purpose:     main file of XMLWrite test program
// Author:      Vadim Zeitlin
// Created:     22.03.04
// RCS-ID:      $Id: main.cpp,v 1.7 2004/10/05 09:13:53 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/tests/testxmlwrite.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(XMLWriteTestCase::suite());

  return runner.run("") ? 0 : 1;
}

