/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/internatcoherencetests.cpp
// Purpose:     Base class for internal coherence tests
// Author:      Ito33 Canada
// Created:     2006/06/13
// RCS-ID:      $Id: internalcoherencetests.cpp,v 1.8 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/internalcoherence.cpp
    @brief Base class for testing IHG projects

    Base class for testing IHG projects
*/


#include "ito33/beforestd.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"
#include "ito33/constants.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/exercisetype.h"
#include "ito33/finance/optiontype.h"

#include "internalcoherencetests.h"
#include "testparam.h"
#include "utiltest.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33 
{

namespace ihg 
{

namespace test
{




//-----------------------------------------------------------------------------
// Put Call Parity Check
//-----------------------------------------------------------------------------
bool OptionPutCallParity(shared_ptr<OptionInterface> pOptionInterface, 
                         TestParam testParam, ito33::XML::RootTag &tag)
{
 
  std::string sTestTitle("putcallparity");

  //create test tags
  XML::Tag test("test", tag);
  test.Comment("Check that Put Call parity " 
          "is verified: C - P = S exp(- r *(Dividends) ) - K e^(-r T)");
  test.Element("name")("option_put_call_parity");
  test.Element("file")(testParam.GetFileName());

  bool bResult = false;

  pOptionInterface->SetExerciseType( finance::ExerciseType_European );

  //remove dividends
  pOptionInterface->SetDividends( shared_ptr<finance::Dividends>(new finance::Dividends()) );

  //compute call price
  pOptionInterface->SetOptionType( finance::Option_Call );
  pOptionInterface->Solve();
  double dPriceCall = pOptionInterface->GetPrice();

  //compute put price
  pOptionInterface->SetOptionType( finance::Option_Put );
  pOptionInterface->Solve();
  double dPricePut = pOptionInterface->GetPrice();

  //PutCallParity :  C + K exp(-rT) = P + S
  double dS      = pOptionInterface->GetSpotSharePrice();
  double dStrike = pOptionInterface->GetStrike();

  //Get numerical yield curve
  shared_ptr<finance::YieldCurve> 
    pYieldCurve = pOptionInterface->GetYieldCurve();

  //Get numerical foreign curve
  shared_ptr<finance::YieldCurve> 
    pForeignCurve = pOptionInterface->GetForeignCurve();

  //get the discounted yield rate
  double dMaturityDate 
    = pOptionInterface->GetMaturityDate().GetExcel()*ONEDAY;

  double dValuationDate 
    = pOptionInterface->GetValuationDate().GetExcel()*ONEDAY;

  Array<double> whatever(2);
  whatever[0] = dValuationDate;
  whatever[1] = dMaturityDate;

  Array<double> values(2);

  pYieldCurve->GetDiscountFactor(whatever.Get(),values.Get(),2);
  double dYieldRate = values[1] / values[0];

  //get the discounter foreign rate
  pForeignCurve->GetDiscountFactor(whatever.Get(),values.Get(),2);
  double dForeignRate = values[1] / values[0];

  //do the actual test C - P = S - PV(Dividends) - K e^(-r T),
  double dPutCallPrice  = dPriceCall - dPricePut - (dS*dForeignRate - dStrike*dYieldRate);
 
  ito33::XML::ParameterOneValue paramvalue;
  paramvalue.m_pdParameter.push_back(0);
  paramvalue.m_pdValue.push_back(dPutCallPrice);
  
  double dTolerance = 1.e-3 * pOptionInterface->GetSpotSharePrice();

  if ( fabs( dPutCallPrice ) < dTolerance ) //weaker parameters 
  {
    test.Element("result")("pass");
    bResult = true;
  } 
  else 
  {
    pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(), sTestTitle);
    test.Element("result")("fail");
    bResult = false;
  }
  
  test.Element("summary",paramvalue);

  return bResult;
}

//-----------------------------------------------------------------------------
// Generic testing system
//-----------------------------------------------------------------------------
bool CoherenceTest(shared_ptr<OptionInterface> pOptionInterface, 
   ito33::XML::RootTag &tag,  std::string sTestTitle, std::string sTestComment, 
    TestType testType, TestDirection testDirection, 
    TestDirection testMovement, TestParam testParam)
 {
   XML::Tag localtest("coherence_test", tag);
   localtest.Comment(sTestComment);
   localtest.Element("name")(sTestTitle);
   localtest.Element("file")(testParam.GetFileName());

   //Test Parameters
   int nDirTest = 0; //indicate which the variable should increase or decrease
   int nDirMov  = 0; // indicate which way the 
                     // price is expected to go increase or decrease

   double dTestBound; //indicate max or minimun value to be reached before stopping
   double dStep         = 0;
   double dMax          = 0;
   double dMin          = 0;
   double dCurrentParam = 0;
   double dPriceOld     = 0 ;
   bool bResult         = false;

   ito33::XML::ParameterTwoValues paramValue;
  
   pOptionInterface->GetTestParameters(dStep, dMax, dMin, testParam, testType);
   
  //direction of the test determine the starting value
  ( testDirection == INCREASE )? dCurrentParam = dMin : dCurrentParam = dMax;
  
  //direction of the test increase or decrease in the parameter
  (testDirection == INCREASE) ? nDirTest = 1: nDirTest = -1;

  //maximum value of test to be reached
  (testDirection == INCREASE) ? dTestBound = dMax : dTestBound = dMin;

  //expected change in value determines the first value
  (testMovement == INCREASE ) ? dPriceOld    = 0 : dPriceOld = 1.e8;

  //expected change in the value of the price increase or decrease
  (testMovement == INCREASE) ? nDirMov = 1: nDirMov = -1;

   
   for (;;)     
   {

     if ( nDirTest>0 && dCurrentParam > dTestBound) 
     {
       break;
     }
     else if ( nDirTest < 0 &&  dCurrentParam < dTestBound )
     {
       break;
     }

     pOptionInterface->SetTestParameter(dCurrentParam, testType);

     //solve
     pOptionInterface->Solve();
     double dPrice = pOptionInterface->GetPrice();

     /*
     std::cout << "ndirtest = " << nDirTest << " dCurrentParam " << dCurrentParam
               << " dTestBound " << dTestBound << " dPrice " << dPrice 
               << " dPriceOld " << dPriceOld  << std::endl;
     */
    
     
     double dError = fabs(dPrice - dPriceOld);
     dError = dPrice > 1. ? dError/dPrice: dError; 

     paramValue.m_pdParameter.push_back(  dCurrentParam );
     paramValue.m_pdValue1.push_back( dPrice ) ;
     paramValue.m_pdValue2.push_back( dPriceOld ) ;
    

     if ( (nDirMov > 0 && dPrice < dPriceOld && dError > DOUBLETOLERANCE ) || 
         (nDirMov < 0 && dPrice > dPriceOld  && dError > DOUBLETOLERANCE ) ) 
      {           
        pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(), sTestTitle);
        localtest.Element("result")("fail");
        bResult = false;
        break;
      } 
    else
      bResult = true;
        
    //do the next step
    dCurrentParam = dCurrentParam + nDirTest*dStep;
    dPriceOld = dPrice;
   } //end while
    
    if ( bResult == true ) 
      localtest.Element("result")("pass");

    //output details
    localtest.Element("summary",paramValue);

    return bResult;

} //end generic testing

//-----------------------------------------------------------------------------
// Delta for in the money option is equal to One
//-----------------------------------------------------------------------------
bool OptionDeltaInTheMoneyEqualToOne(shared_ptr<OptionInterface> pOptionInterface, 
                  TestParam testParam, ito33::XML::RootTag &tag)
{
  // right now, we don't need check this for exotic options
  // as this is not always true
  XML::Tag localtest("test", tag);
  localtest.Comment("Check that the delta for in the money  is" 
                            " equal to one");
  localtest.Element("name")("delta_in_the_money_equal_to_one");
  localtest.Element("file")(testParam.GetFileName());

  //add no dividend as dividend may erode the value
  pOptionInterface->SetDividends( shared_ptr<finance::Dividends>(new finance::Dividends()) );

  //remove the borrow/foreign curve 
  pOptionInterface->SetForeignCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(0.0)) );

  double dStrike = pOptionInterface->GetStrike();
  double dSpot   = 0.0;
 
  if ( pOptionInterface->GetOptionType() == finance::Option_Call )
   {
     dSpot = 50.*dStrike;
   }
   else if ( pOptionInterface->GetOptionType() == finance::Option_Put )
   {
     dSpot = 0.01*dStrike;
   }

   pOptionInterface->SetSpotSharePrice(dSpot);
  
   pOptionInterface->Solve();
   double dDelta = pOptionInterface->GetDelta();

   ito33::XML::ParameterOneValue paramValue;
   paramValue.m_pdParameter.push_back(0);
   paramValue.m_pdValue.push_back(dDelta);

   localtest.Element("summary",paramValue);

   if ( fabs(fabs(dDelta) - 1.0)/fabs(dDelta) > testParam.m_dEPSCheck  )
   {
     
     pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
       "delta_in_the_money_equal_to_one");
     localtest.Element("result")("fail");
     return false;
   }

   //close the XML output
   localtest.Element("result")("pass");
 
   return true;
} // Testinternalcoherence::DeltaInTheMoneyEqualToOne()


//-----------------------------------------------------------------------------
// Delta for out  the money option is equal to zero
//-----------------------------------------------------------------------------
  bool OptionDeltaOutOfTheMoneyEqualToZero(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag)
{

  std::string testname("deltaoutofthemoneyequaltozero");
  XML::Tag localtest("test", tag);
  localtest.Comment("Check that the delta for out the money option is close to Zero");
  localtest.Element("name")("deltaoutofthemoneyequaltozero");
  localtest.Element("file")(testParam.GetFileName());

  //add no dividend as dividend may erode the value
  pOptionInterface->SetDividends( shared_ptr<finance::Dividends>(new finance::Dividends()) );

  //remove the borrow/foreign curve 
  pOptionInterface->SetForeignCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(0.0)) );

   double dStrike = pOptionInterface->GetStrike();
   double dSpot   = 0; // want to be realy out of the money 
   if ( pOptionInterface->GetOptionType() == finance::Option_Call )
   {
     dSpot = .01*dStrike;
   }
   else if ( pOptionInterface->GetOptionType() == finance::Option_Put )
   {
     dSpot = 50.*dStrike;
   }

   pOptionInterface->SetSpotSharePrice(dSpot);
   
   pOptionInterface->Solve();
   double dDelta = pOptionInterface->GetDelta();

   ito33::XML::ParameterOneValue paramValue;
   paramValue.m_pdParameter.push_back(0);
   paramValue.m_pdValue.push_back(dDelta);

   localtest.Element("summary",paramValue);

   if ( fabs(dDelta) > testParam.m_dEPSCheck  )
   {
     pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
       "delta_out_of_the_money_equal_to_zero");
     localtest.Element("result")("fail");
     return false;
   }


   //close the XML output
   localtest.Element("result")("pass");
 
   return true;
} // Testinternalcoherence::DeltaOutOfTheMoneyEqualToZero()


//-----------------------------------------------------------------------------
// Gamma for in and out  the money option is equal to zero
//-----------------------------------------------------------------------------
bool OptionGammaOutInOfTheMoneyEqualToZero(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag)
{
  XML::Tag localtest("test", tag);
  localtest.Comment("Check that the gamma for in and out of the money put" 
                     "and call is equal to zero");
  localtest.Element("name")("gammaoutinthemoneyequaltozero");
  localtest.Element("file")(testParam.GetFileName());
  
   double dStrike = pOptionInterface->GetStrike();
   double dSpot   = 0.0; 

   if ( pOptionInterface->GetOptionType() == finance::Option_Call )
   {
     dSpot = 50.*dStrike;
   }
   else if ( pOptionInterface->GetOptionType() == finance::Option_Put )
   {
     dSpot = 0.1*dStrike;
   }

  pOptionInterface->SetSpotSharePrice(dSpot);

  pOptionInterface->Solve();
  double  dGamma = pOptionInterface->GetGamma();

  ito33::XML::ParameterOneValue paramValue;
  paramValue.m_pdParameter.push_back(0);
  paramValue.m_pdValue.push_back(dGamma);

  if ( fabs(dGamma) >= testParam.m_dEPSCheck )
  {
    pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
      "gamma_in_the_money_equal_to_zero");
    localtest.Element("result")("fail");
    localtest.Element("summary",paramValue);
    return false;
  }
  
   //out of the money
   if ( pOptionInterface->GetOptionType() == ito33::finance::Option_Call )
   {
     dSpot = .1*dStrike;
   }
   else if ( pOptionInterface->GetOptionType() == ito33::finance::Option_Put )
   {
     dSpot = 50.*dStrike;
   }

   pOptionInterface->SetSpotSharePrice(dSpot);

   pOptionInterface->Solve();
   dGamma = pOptionInterface->GetGamma();
   paramValue.m_pdParameter.push_back(0);
   paramValue.m_pdValue.push_back(dGamma);

   if ( fabs(dGamma) >= testParam.m_dEPSCheck )
   { 
     pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
      "gamma_out_of_the_money_equal_to_zero");
     localtest.Element("result")("fail");
     localtest.Element("summary",paramValue);
     return false;
   }

   localtest.Element("result")("pass");
   localtest.Element("summary",paramValue);

   return true;

} // GammaOutInOfTheMoneyEqualToZero()


//-----------------------------------------------------------------------------
// Theta is negative when risk free rate is zero
//-----------------------------------------------------------------------------
bool OptionThetaNegativeRiskFreeRateZero(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag)
{
  XML::Tag localtest("test", tag);
  localtest.Comment("Theta negative when risk free rate is zero");
  localtest.Element("name")("thetanegativewhenriskfreeratezero");
  localtest.Element("file")(testParam.GetFileName());

  //Set the borrow curve, yield curve to zero
  pOptionInterface->SetYieldCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(0.0)) );

  //set the spot and strike to be the same
 pOptionInterface->SetStrike( pOptionInterface->GetSpotSharePrice() );

  // output pricing
 pOptionInterface->Solve();
 double dTheta = pOptionInterface->GetTheta();
  
 ito33::XML::ParameterOneValue paramValue;
 paramValue.m_pdParameter.push_back(0);
 paramValue.m_pdValue.push_back(dTheta);


  if ( dTheta > 0 )
  {    
    pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
      "theta_negative_risk_free_rate_zero");
    localtest.Element("result")("fail");
    return false;
  }
 
  localtest.Element("result")("pass");
  localtest.Close();
  return true;

} // end ThetaNegativeRiskFreeRateZero


//-----------------------------------------------------------------------------
// Call goes to S as Strike goes to zero
// no dividends
//-----------------------------------------------------------------------------
bool OptionCallGoesToSWhenStrikeDecreases(shared_ptr<OptionInterface> pOptionInterface, 
              TestParam testParam, ito33::XML::RootTag &tag)
{
  XML::Tag localtest("coherence_test", tag);
  localtest.Comment("Call goes to S when  strike decreases");
  localtest.Element("name")("callgoestoswhenstrikedecreases");
  localtest.Element("file")(testParam.GetFileName());

  //Set it to be a call
  pOptionInterface->SetOptionType( finance::Option_Call );

  //Remove the dividends
  pOptionInterface->SetDividends( shared_ptr<finance::Dividends>(new finance::Dividends()) );

  double dStep = testParam.m_dStrikeStep * pOptionInterface->GetStrike();
  double dMax  = testParam.m_dStrikeMax  * pOptionInterface->GetStrike();
  double dMin  = testParam.m_dStrikeMin  * pOptionInterface->GetStrike();
  bool bResult = true;

  double dSpot = pOptionInterface->GetSpotSharePrice();
  double dDistancePricetoSOld = dSpot;
  double dDistancePricetoS = 0;
  double dCurrentParam = dMax;
 
  ito33::XML::ParameterTwoValues paramValue;
 
  while ( dCurrentParam > dMin)
  {   
    pOptionInterface->SetStrike(dCurrentParam);

    // output pricing
    pOptionInterface->Solve();
    double dPrice = pOptionInterface->GetPrice();

    paramValue.m_pdParameter.push_back(dCurrentParam);
    paramValue.m_pdValue1.push_back(dPrice);
    paramValue.m_pdValue2.push_back(dSpot);

    dDistancePricetoS = dSpot - dPrice;

    if ( dDistancePricetoS > dDistancePricetoSOld )
    {
      pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
      "call_goes_to_when_strike_goes_to_zero");
      localtest.Element("result")("fail");
      bResult = false;
      break;
    }

    //do the next step
    dCurrentParam     = dCurrentParam - dStep;
    dDistancePricetoSOld = dDistancePricetoS;
  } //end while

      
    if ( bResult) 
    {
      localtest.Element("result")("pass");
    }

    //output details
    localtest.Element("summary",paramValue);

    return bResult;

} //end CallGoesToSWhenStrikeDecrease




  } //end namespace test
 } //end namespace ito33
} //namespace ihg
