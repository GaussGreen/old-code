/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/main.cpp
// Purpose:     main file of bondlike test
// Author:      Zhang
// Created:     24.06.04
// RCS-ID:      $Id: main.cpp,v 1.12 2005/05/24 10:36:00 zhang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testbondlike.h"
#include "ito33/tests/testyieldtomaturity.h"
#include "ito33/tests/testyieldtoPut.h"
#include "ito33/tests/testnewshare.h"


int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(BondLikeTest::suite());
  runner.addTest(YieldToMaturityTest::suite());
  runner.addTest(YieldToPutTest::suite());
  runner.addTest(NewShareTest::suite());

  return runner.run("") ? 0 : 1;
}

