/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/convergencetests.h
// Purpose:     test that the numerical convergence property are valid
// Author:      Ito33Canada
// Created:     2005/06/13
// RCS-ID:      $Id: convergencetests.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/option/testsuite/numericaltests.h
    @brief Define a set of test to validate  the numerical
           convergence property 
*/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_NUMERICALTESTS_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_NUMERICALTESTS_H_

#include "ito33/xml/write.h"

#include "testparam.h"

namespace ito33 
{

namespace finance
{
  class Derivative;
}
namespace ihg
{
  class TheoreticalModel;

namespace test
{

 bool DerivativeCheckConvergence(shared_ptr<TheoreticalModel> pModel, 
   shared_ptr<finance::Derivative> pDeriv, TestParam testParam, 
   ito33::XML::Tag &tag, size_t nNbTests);

 
 } //end test
} //end ihg
} //end ito33

#endif

