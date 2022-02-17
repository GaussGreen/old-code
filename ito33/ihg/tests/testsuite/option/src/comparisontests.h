/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/src/option/comparisontests.h
// Purpose:     Base class for testing IHG projects
//              against the black Sholes solution 
// Author:      Ito33Canada
// Created:     2005/06/10
// RCS-ID:      $Id: comparisontests.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/src/option/comparisontests.h
    @brief Base class for testing IHG projects
    Results should be compared with the Black-Scholes formulas for the
    American and European options with flat rates and no discrete dividends.
*/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_COMPARISONTEST_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_COMPARISONTEST_H_

#include "ito33/sharedptr.h"

#include "ito33/xml/write.h"

#include "optioninterface.h"
#include "testparam.h"
#include "utiltest.h"

namespace ito33 
{

namespace ihg
{

namespace test
{
  

  /**
    Generic test function to compare the black scholes
    solution when different input parameters are moving

    @param pOptionInterface option object
    @param TestParam control parameters for the different tests
    @param RootTag to create the xml output
    @param testtitle xml testtitle
    @param testcomment paramtestcomment
    @param testType testtype
    @param OptionType put call
    @return true success || false failure
  */
  bool BlackScholesComparison(shared_ptr<OptionInterface> pOptionInterface,
    ito33::XML::RootTag &tag, TestParam testParam,
    std::string sTestTitle, std::string sTestComment,
    TestType testType);

 
 } //end test
} //end ihg
} //end ito33

#endif // _IHG_TESTS_TESTSUITE_OPTION_SRC_COMPARISONTEST_H_
