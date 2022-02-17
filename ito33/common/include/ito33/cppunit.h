/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/cppunit.h
// Purpose:     wrapper header for CppUnit headers
// Author:      Vadim Zeitlin
// Created:     26.06.03
// RCS-ID:      $Id: cppunit.h,v 1.5 2006/05/15 15:18:14 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/cppunit.h
    @brief Central header for all CppUnit headers which our tests use.

    Normally you need to include at least 4 headers to use CppUnit and, worse,
    all of them result in warnings at level 4 with VC++ so you also have to
    include before/afterstd.h before/after including them which makes 6
    lines just for the headers. With this helper header just one line is
    enough.
 */

#ifndef _ITO33_CPPUNIT_H_
#define _ITO33_CPPUNIT_H_

// using CPPUNIT_TEST() macro results in this warning, disable it as there is
// no other way to get rid of it and it's not very useful anyhow
#ifdef _MSC_VER
  // typedef-name 'foo' used as synonym for class-name 'bar'
  #pragma warning(disable:4097)

  // unreachable code: we don't care about warnings in CppUnit headers
  #pragma warning(disable:4702)
#endif // _MSC_VER

#include "ito33/beforestd.h"
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "ito33/afterstd.h"

#ifdef _MSC_VER
  #pragma warning(default:4702)
#endif // _MSC_VER

namespace ito33
{

/// Precision used for comparing double values with the expected result.
static const double ITO33_DOUBLE_CMP_PRECISION = 1e-10;

/**
    Use this macro to compare two doubles with the default precision.

    @param dValue the real value
    @param dExpected the expected value
 */
#define ITO33_ASSERT_DOUBLES_EQUAL(dValue, dExpected) \
    CPPUNIT_ASSERT_DOUBLES_EQUAL((dValue), (dExpected), \
                                 ITO33_DOUBLE_CMP_PRECISION)

} // namespace ito33

// for VC++ automatically link in cppunit library
#ifdef _MSC_VER
  #ifdef NDEBUG
    #pragma comment(lib, "cppunit.lib")
  #else // Debug
    #pragma comment(lib, "cppunitd.lib")
  #endif // Release/Debug
#endif

/**
    A trivial main function implementation for a simple test.

    For the simple unit tests, the main function just runs the single test case
    and this macro implements such main. To use it the test must be registered
    using CPPUNIT_TEST_SUITE_REGISTRATION() macro.
 */
#define ITO33_TEST_MAIN()                                                     \
  int main()                                                                  \
  {                                                                           \
    CppUnit::TextUi::TestRunner runner;                                       \
    runner.addTest(CppUnit::TestFactoryRegistry::getRegistry().makeTest());   \
    return runner.run("") ? 0 : 1;                                            \
  }

#endif // _ITO33_CPPUNIT_H_

