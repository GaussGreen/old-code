/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/internalcoherencetest.h
// Purpose:     Base class for testing IHG projects
// Author:      Ito33 Canada
// Created:     2006/06/13
// RCS-ID:      $Id: internalcoherencetests.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/internalcoherencetest.h
    @brief Define a set of test to validate financial
    properties of the solver 
       * The option definition works
       * Put Call Parity
       * price increases when spot increases
       * price decreases when spot decreases
       * Call price increases when volatility increases
       * Put price decreases when volatility increases
       *  price decreases with risk free rate
       * Call price decreases when foreign rate increases
       * Put price increases when foreign rate increases
       * price increases when hazard rate increases
       * Delta for in the money option is 1
       * Delta for out of the money option is 0
       * Gamma for in the money and out the money is 0
       * American  option >= European 
    Base class for testing IHG projects
*/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_INTERNALCOHERENCETESTS_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_INTERNALCOHERENCETESTS_H_

#include "ito33/xml/write.h"

#include "optioninterface.h"
#include "testparam.h"

namespace ito33 
{

namespace ihg
{


namespace test
{
 
  /**
    @param TestParam control parameters for the different tests
    @param pOptionInterface option object
    @param RootTag to create the xml output
    @return true success || false failure
   */
   bool OptionPutCallParity(shared_ptr<OptionInterface> pOptionInterface, 
                         TestParam testParam, ito33::XML::RootTag &tag);


  /**
  Delta for in the money option is equal to One

  @param TestParam control parameters for the different tests
  @param pOptionInterface option object
  @param RootTag to create the xml output
  @return true success || false failure
  */
  bool OptionDeltaInTheMoneyEqualToOne(shared_ptr<OptionInterface> pOptionInterface, 
                 TestParam testParam, ito33::XML::RootTag &tag);

  /**
  Delta for out the money option is equal to zero

  @param TestParam control parameters for the different tests
  @param pOptionInterface option object
  @param RootTag to create the xml output
  @return true success || false failure
  */
  bool OptionDeltaOutOfTheMoneyEqualToZero(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag);

  /**
  Gamma for in the money and out the money option is equal to one

  @param TestParam control parameters for the different tests
  @param pOptionInterface option object
  @param RootTag to create the xml output
  @return true success || false failure
  */
  bool OptionGammaOutInOfTheMoneyEqualToZero(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag);

  /**
  Theta is negative if risk free rate is zero

  @param TestParam control parameters for the different tests
  @param pOptionInterface option object
  @param RootTag to create the xml output
  @return true success || false failure
  */
  bool OptionThetaNegativeRiskFreeRateZero(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag);


  /**
  Call goes to S as K goes to zero, no dividend

  @param TestParam control parameters for the different tests
  @param pOptionInterface option object
  @param RootTag to create the xml output
  @return true success || false failure
  */
  bool OptionCallGoesToSWhenStrikeDecreases(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag);

    /**
    Generic testing for when a parameter is increased or decreased
    @param TestParam control parameters for the different tests
    @param myOption option object
    @param RootTag to create the xml output
    @param sTestTitle title of the test
    @param stestComment comment of the test xml
    @param testType Volatility, Hazard rate, etc...
    @param testDirection way the parameter to be tested should move, up/down
    @param testMovment way the option price is supposed to go up/down
  */
  bool CoherenceTest(shared_ptr<OptionInterface> pOptionInterface, 
                     ito33::XML::RootTag &tag,  
                     std::string sTestTitle, 
                     std::string sTestComment, 
                     TestType testtype, 
                     TestDirection testDirection, 
                     TestDirection testMovement, 
                     TestParam testParam);


} //end namespace test

} //end namespace ihg

} //end namspace ito33


#endif // _IHG_TESTS_TESTSUITE_OPTION_SRC_INTERNALCOHERENCETESTS_H_
