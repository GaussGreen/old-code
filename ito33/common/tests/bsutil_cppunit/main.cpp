/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bsutil_cppunit/main.cpp
// Purpose:     main file for black-scholes utility function
// Author:      Laurence
// Created:     22/09/2003
// RCS-ID:      $Id: main.cpp,v 1.2 2004/10/05 09:13:49 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testBSut.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(BSutileTest::suite());

  return runner.run("") ? 0 : 1;
}
