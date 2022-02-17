/////////////////////////////////////////////////////////////////////////////
// Name:        tests/volatility/main.cpp
// Purpose:     main file for testing volatility
// Author:      Yann d'Halluin
// Created:     11/07/2004
// RCS-ID:      $Id: main.cpp,v 1.2 2005/06/27 16:23:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ihg/tests/testvolatility.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(VolatilityFlatTest::suite());

  return runner.run("") ? 0 : 1;

  
}

