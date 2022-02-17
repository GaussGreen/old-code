/////////////////////////////////////////////////////////////////////////////
// Name:        test/tocdouble/main.cpp
// Purpose:     main file of ToCDouble/FromCDouble test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: main.cpp,v 1.1 2005/11/29 15:39:53 vaclav Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(CppUnit::TestFactoryRegistry::getRegistry().makeTest());
  return runner.run("") ? 0 : 1;
}
