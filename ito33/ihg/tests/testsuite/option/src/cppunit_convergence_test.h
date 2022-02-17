
/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/cppunit_convergence_test.h
// Purpose:     Base class for convergence testing IHG projects
// Author:      ITO33 Canada
// Created:     2005/06/08
// RCS-ID:      $Id: cppunit_convergence_test.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/cppunit_option_test.h
**/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_CONVERGENCE_TEST_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_CPPUNIT_CONVERGENCE_TEST_H_


#include "ito33/exception.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/derivative.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33 
{
namespace ihg
{
namespace test
{

class CppUnitConvergenceTest: public CppUnit::TestCase
{
 
private:
  shared_ptr<TheoreticalModel> m_pModel;
  shared_ptr<finance::Derivative> m_pDerivative;
  size_t m_nNbTests;

public:
  CppUnitConvergenceTest():m_nNbTests(4) {}

  void setup();


 private:
  CPPUNIT_TEST_SUITE( CppUnitConvergenceTest );
  
   CPPUNIT_TEST( CheckConvergence ); 

  CPPUNIT_TEST_SUITE_END();

  
  void CheckConvergence();

  NO_COPY_CLASS(CppUnitConvergenceTest);

};

} //end namespace test
} //end namespace ihg
}//end namespace ito33

#endif
