/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/optiontester.h
// Purpose:     Test class for regular options
// Author:      David
// Created:     2005/06/15
// RCS-ID:      $Id: optiontester.h,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/optiontester.h
    @brief Test class for regular options

    Test class for regular options
*/

#ifndef _IHG_OPTIONTESTER_H_
#define _IHG_OPTIONTESTER_H_

#include "ito33/array.h"

#include "ihg/basetester.h"

#include "ito33/finance/option.h"

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
  @brief Class for testing vanilla option pricing

  This class provides functions for testing vanilla option (European, 
  American) pricing.

*/
class OptionTester : public ihg::BaseTester
{

public:
  OptionTester() 
  {
    m_nNbRuns = 0;
  }

  virtual ~OptionTester() { }

  /// Run the option pricer. Pure virtual in base
  void RunPricer();

  /// Read option specific data
  void ReadContractParams(std::ifstream& sIn);
  void ReadContractParams(std::string& sContractFile);

  /// Report the converegence test results (prices and Greeks)
  void ReportConvergenceTest();

  /// Run the code and write out data for plotting graphs of the solution  
  void OutputPlotData();

  /// Get the final delta
  double GetDelta();

  /// Get the final gamma
  double GetGamma();

  /// Get the final vega
  double GetVega();

  /// Run tests for the value/convergence of the price, delta, gamma and vega
  bool RunFullTests(double dPrice, double dDelta, double dGamma, 
                    double dVega, size_t nNbTests = 4);

protected:

  /// The base option object, containing the strike, maturity, etc.
  shared_ptr<finance::Option> m_pOptionFinance;

  /// The payoff type.  Needs to be stored to handle constant payoffs 
  int m_iType;


}; // class OptionTester


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_OPTIONTESTER_H_
