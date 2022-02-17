/////////////////////////////////////////////////////////////////////////////
// Name:        test/ycurve/main.cpp
// Purpose:     main file of RWLock test program
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: main.cpp,v 1.3 2004/10/05 09:13:49 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/tests/const_constraint_case.h"

#include "ito33/cppunit.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(ConstConstraintCase::suite());

  return runner.run("") ? 0 : 1;
}
