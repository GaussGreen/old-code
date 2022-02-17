//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLAddin.hpp
//
//   Description : Class for creating and running regression tests for addins
//
//   Author      : Mark A Robson
//
//   Date        : 25 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ADDINTEST_HPP
#define EDG_ADDINTEST_HPP

#include "edginc/Object.hpp"
#include "edginc/RegressionTest.hpp"
#include <stdarg.h>

DRLIB_BEGIN_NAMESPACE

class ADDINS_DLL AddinTest: public CObject,
                 public IRegressionTest {
public:
    friend class AddinTestHelper;
    static CClassConstSP const TYPE;

    /** create an instance of an AddinTest - provides ability to write out
        regression file. Note does not take a copy of params */
    AddinTest(const string&          addinName,
              const IObjectSP&       params);

    /** Run regression test for this addin */
    IObjectSP runTest() const;

private:
    string          addinName;
    IObjectSP       params;

    AddinTest();
    AddinTest(const AddinTest &rhs);
    AddinTest& operator=(const AddinTest& rhs);
    
};

DRLIB_END_NAMESPACE
#endif
