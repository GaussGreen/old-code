/////////////////////////////////////////////////////////////////////////////
// Name:        tests/nonlinearsolver/testNLS.cpp
// Purpose:     main file for unit tests for non linear equation 1D solvers
// Author:      Pedro Ferreira, Laurence Gozalo
// Created:     04.12.03
// RCS-ID:      $Id: main.cpp,v 1.2 2004/10/05 09:13:51 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testNLS.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(NLSTest::suite());

  return runner.run("") ? 0 : 1;
}
