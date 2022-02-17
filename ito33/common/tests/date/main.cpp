/////////////////////////////////////////////////////////////////////////////
// Name:        test/date/main.cpp
// Purpose:     main file of Date test program
// Author:      Vadim Zeitlin
// Created:     2004-07-08
// RCS-ID:      $Id: main.cpp,v 1.5 2006/05/11 10:38:02 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/date.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testdate.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// helper macro
// ----------------------------------------------------------------------------

#define ASSERT_EQUAL(n, s1, s2)                                               \
        CPPUNIT_ASSERT_MESSAGE( Str::Printf("test %lu: \"%s\" != \"%s\"",     \
                                            (unsigned long)n,                 \
                                            s1.c_str(),                       \
                                            s2),                              \
                                s1 == s2 )

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
    Date dt = Date::Today();
    printf("Today is %s\n", dt.Format().c_str());

    CppUnit::TextUi::TestRunner runner;
    runner.addTest(DateTestCase::suite());

    return runner.run("") ? 0 : 1;
}


