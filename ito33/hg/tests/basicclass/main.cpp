/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/basicclass/main.cpp
// Purpose:     main file for testing HG basic classes
// Created:     2005/04/26
// RCS-ID:      $Id: main.cpp,v 1.2 2005/05/02 12:44:41 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "translator.h"
#include "theoreticalmodel.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(TranslatorTest::suite());

  runner.addTest(TheoreticalModelTest::suite());

  return runner.run("") ? 0 : 1;
}
