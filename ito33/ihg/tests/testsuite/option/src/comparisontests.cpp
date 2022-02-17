/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/src/comparison.cpp
// Purpose:     Test for Comparison
// Author:      Ito33Canada
// Created:     2005/06/10
// RCS-ID:      $Id: comparisontests.cpp,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/src/comparisontests.cpp
    @brief Base class for testing IHG projects
    Results should be compared with the Black-Scholes formulas for the
    American and European options with flat rates and no discrete dividends.
*/

#include "ito33/beforestd.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"
#include "ito33/string.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"

#include "comparisontests.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33 
{

namespace ihg 
{

namespace test
{

//-----------------------------------------------------------------------------
// Generic test
//-----------------------------------------------------------------------------
bool BlackScholesComparison(shared_ptr<OptionInterface> pOptionInterface,
                            XML::RootTag &tag, TestParam testParam,
                            std::string sTestTitle, std::string sTestComment,
                            TestType testType)
{

  XML::Tag localtest("test_blackscholes",tag);
  localtest.Comment(sTestComment);
  localtest.Element("name")(sTestTitle);
  localtest.Element("file")(testParam.GetFileName());

  //override the exercise type to be european
  pOptionInterface->SetExerciseType(  finance::ExerciseType_European );

  //ensure that there are not hazard rate
  pOptionInterface->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(0.0)) );

  //ensure that there are no discrete dividends
  pOptionInterface->SetDividends( shared_ptr<finance::Dividends>(new finance::Dividends()) );

  //ensure that volatility is constant
  pOptionInterface->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat( .2 )) );
  
  //ensure that yield rate are constant
  pOptionInterface->SetYieldCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat( .05)) );

  //ensure that foreign rate are constant
  pOptionInterface->SetForeignCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat( .0 )) );

  //Test Parameters
  double dTestBound; //indicate max or minimun value to be reached before stopping
  double dStep = 0;
  double dMax  = 0;
  double dMin  = 0;
  double dCurrentParam = 0;
  bool bHasSucceeded   = false;

  XML::ParameterTwoValues paramValue;
 
  pOptionInterface->GetTestParameters(dStep, dMax, dMin, testParam, testType);
 
  /*std::cout << "dStep : " << dStep << std::endl;
  std::cout << "dMax  : " << dMax << std::endl;
  std::cout << "dMin  : " << dMin << std::endl;*/

  //tests are going from dMin to dMax
  dCurrentParam = dMin ;
  dTestBound    = dMax;
  double dTolerance  = 1.e-3;

  while ( dCurrentParam < dTestBound )
  {   
    pOptionInterface->SetTestParameter(dCurrentParam,testType);
      
    pOptionInterface->Solve();
    double dPrice = pOptionInterface->GetPrice();

    //compute the Comparison Price
    double dBlsPrice = pOptionInterface->AnalyticalSolve(); 

    paramValue.m_pdParameter.push_back(  dCurrentParam );
    paramValue.m_pdValue1.push_back( dPrice ) ;
    paramValue.m_pdValue2.push_back( dBlsPrice ) ;

    double dError = 0.0;

    if (  dPrice < 1.0 )
      dError = fabs( dPrice - dBlsPrice  );
    else
      dError = fabs( dPrice/dBlsPrice - 1.);

   
    if (  dBlsPrice > 0 && 
       dError > dTolerance * pOptionInterface->GetSpotSharePrice() ) 
    {
      bHasSucceeded = false;
      pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(), sTestTitle);
      break;
    } 
    else
    {
      bHasSucceeded = true;
    }


    //do the next step
    dCurrentParam = dCurrentParam + dStep;
  } //end while
  
  if ( bHasSucceeded ) 
    localtest.Element("result")("pass");
  else
    localtest.Element("result")("fail");


  //output details
  localtest.Element("summary",paramValue);

  return bHasSucceeded;

} //end generic testing

} //end test namespace
} //end ihg namespace
}//end ito33 namespace
