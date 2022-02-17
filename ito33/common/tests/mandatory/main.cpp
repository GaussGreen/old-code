/////////////////////////////////////////////////////////////////////////////
// Name:        tests/mandatory/main.cpp
// Purpose:     main file of mandatory test
// Author:      Ito33 Canada
// Created:     May 11, 2005
// RCS-ID:      $Id: main.cpp,v 1.2 2005/06/27 13:20:47 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testmandatory.h"


int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(PepsAveragingPeriodTest::suite());

  return runner.run("") ? 0 : 1;
}

