/////////////////////////////////////////////////////////////////////////////
// Name:        tests/crosscurrency/main.cpp
// Purpose:     main file for testing cross-currency
// Author:      Nabil Ouachani
// Created:     2006/03/02
// RCS-ID:      $Id: main.cpp,v 1.2 2006/03/08 13:08:56 nabil Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ihg/tests/testcrosscurrency.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(CrossCurrencyTest::suite());

  return runner.run("") ? 0 : 1;

  
}
