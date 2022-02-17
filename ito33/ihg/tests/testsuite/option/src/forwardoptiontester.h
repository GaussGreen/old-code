/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/forwardoptiontester.h
// Purpose:     Test class for forward pricing of European option
// Author:      Wang
// Created:     2004/03/10
// RCS-ID:      $Id: forwardoptiontester.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/forwardoptiontester.h
    @brief Test class for forward pricing of European option

    Test class for forward pricing
 */

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_FORWARDOPTIONTESTER_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_FORWARDOPTIONTESTER_H_

#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/array.h"

#include "ito33/finance/option.h"

#include "ihg/basetester.h"

namespace ito33 
{
  
namespace finance
{
  enum OptionType;
  enum ExerciseType;
}

namespace ihg 
{

/**
  @brief Class for testing forward pricing

  This class provides functions for testing forward pricing.

*/

class ForwardOptionTester : public ihg::BaseTester
{

public: 
  ForwardOptionTester() 
  {
    m_nNbRuns = 0;
  }

  // Default dtor is ok

  /// Run the forward option pricer. Pure virtual in base
  void RunPricer();

  /// Read forward option specific data
  void ReadContractParams(std::ifstream& sIn);
  void ReadContractParams(std::string& sContractFile);
  
  /// Report the converegence test results (prices and Greeks)
  void ReportConvergenceTest();


protected:

  /// The list of options to price via the forward equations
  std::list< shared_ptr<finance::Option> > m_listOfOptions;

}; // class ForwardOptionTester


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_FORWARDOPTIONTESTER_H_
