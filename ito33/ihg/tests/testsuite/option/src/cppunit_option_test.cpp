
/////////////////////////////////////////////////////////////////////////////
// Name:       ihg/tests/testsuite/option/src/cppunit_option_test.cpp
// Purpose:     Base class for testing IHG projects
// Author:      Yann d'Halluin
// Created:     2005/06/13
// RCS-ID:      $Id: cppunit_option_test.cpp,v 1.9 2006/08/22 15:26:31 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/cppunit_option_test.cpp
**/
#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include <string>
#include "ito33/afterstd.h"


#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/exercisetype.h"

#include "ihg/xml/pricingreader.h"

#include "ito33/ihg/hazardratecallback.h"
#include "ito33/ihg/volatilitycallback.h"

#include "optiontester.h"
#include "cppunit_option_test.h"
#include "testparam.h"
#include "optioninterface.h"
#include "comparisontests.h"
#include "internalcoherencetests.h"

extern ito33::XML::RootTag root;
extern ito33::ihg::test::TestParam testParam;

extern const ito33::Error ITO33_UNEXPECTED;

extern void __stdcall parametricHR(double dTime, const double *pdS, double *pdValues, 
                            size_t nNbS, int i);

extern void __stdcall parametricVol(double dTime, const double *pdS, double *pdValues, 
                            size_t nNbS, int i);

namespace ito33 
{
  

namespace ihg
{

namespace test
{

void CppUnitOptionTest::setup()
{
  XML::PricingReader reader( testParam.GetFileName() );
  shared_ptr<finance::SessionData> pSessionData(reader.ReadSessionData());

  finance::DerivativeVisitorGoodType visitor;
  reader.ReadDerivatives(visitor);
  shared_ptr<finance::Option> pOption = visitor.GetOption();
  
  if(!pOption)
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no option in input xml file");

  shared_ptr<ito33::ihg::TheoreticalModel> 
      pModel(new ito33::ihg::TheoreticalModel);
  reader.ReadTheoreticalModel(pModel);

  m_pOptionInterface = make_ptr( new OptionInterface
                                     (pSessionData, pOption, pModel) );

}
  
void CppUnitOptionTest::PriceAfterMaturityGeneratesAnError()  
{
  setup();

  Date maturityDate(2003, Date::Jan, 1);
  Date valuationDate(2004, Date::Jan, 1);

  m_pOptionInterface->SetMaturityDate( maturityDate );
  m_pOptionInterface->SetValuationDate( valuationDate );

  m_pOptionInterface->Solve();
}

void CppUnitOptionTest::PriceAlwaysPositive()
{
  setup();

  m_pOptionInterface->Solve();

  if ( m_pOptionInterface->GetPrice() < 0 )
    m_pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
    "price_always_positive");

  CPPUNIT_ASSERT( m_pOptionInterface->GetPrice() >= 0 );
}

void CppUnitOptionTest::GammaAlwaysPositive()
{
  setup();

  m_pOptionInterface->Solve();

  const double dCheck = -1.e-12;
  if ( m_pOptionInterface->GetGamma() < dCheck )
    m_pOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
    "price_always_positive");

  CPPUNIT_ASSERT( m_pOptionInterface->GetGamma() >= dCheck );
}

void CppUnitOptionTest::CallIncreasesWhenSpotIncreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Call );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "call_increase_when_spot_increase", "", SPOT, INCREASE, INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
} // end CallIncreaseWhenSpotIncrease()

void CppUnitOptionTest::CallDecreasesWhenSpotDecreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Call );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "call_decrease_when_spot_decrease", "", SPOT, DECREASE, DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
} // end CallIncreaseWhenSpotIncrease()


void CppUnitOptionTest::PutDecreasesWhenSpotIncreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Put);

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "put_decrease_when_spot_increase", "", SPOT, INCREASE, DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::PutIncreasesWhenSpotDecreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Put);

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "put_increase_when_spot_decrease", "", SPOT, DECREASE, INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}


void CppUnitOptionTest::PriceIncreasesWhenVolatilityIncreases()
{
  setup();

  m_pOptionInterface->SetStrike( m_pOptionInterface->GetSpotSharePrice() );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "price_increase_when_volatility_increase", "", VOL, INCREASE, INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}
  
void CppUnitOptionTest::PriceDecreasesWhenVolatilityDecreases()
{
  setup();

  m_pOptionInterface->SetStrike( m_pOptionInterface->GetSpotSharePrice() );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "price_decrease_when_vol_decrease", "", VOL, DECREASE, DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}



void CppUnitOptionTest::AmericanGreaterEuropean()
{
  setup();

  m_pOptionInterface->SetExerciseType( finance::ExerciseType_European);

  m_pOptionInterface->Solve();

  double dPriceEuropean = m_pOptionInterface->GetPrice();
 
  m_pOptionInterface->SetExerciseType( finance::ExerciseType_American);

  m_pOptionInterface->Solve();

  double dPriceAmerican = m_pOptionInterface->GetPrice();

  CPPUNIT_ASSERT( dPriceEuropean <= dPriceAmerican );
}

void CppUnitOptionTest::PriceIncreasesWhenHazardRateIncreases()
{
  setup();

  m_pOptionInterface->SetStrike( m_pOptionInterface->GetSpotSharePrice() );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "price_increase_when_hazardrate_increase", "", HAZARDRATE, INCREASE, INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::PriceIncreasesWhenMaturityIncreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Call );
  m_pOptionInterface->SetStrike( m_pOptionInterface->GetSpotSharePrice() );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "price_increase_when_maturity_increase", "", MATURITY, INCREASE, INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::PriceDecreasesWhenRiskFreeRateDecreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Call );
  m_pOptionInterface->SetStrike( m_pOptionInterface->GetSpotSharePrice() );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "price_decrease_when_risk_free_rate_decrease", "", YIELDRATE, DECREASE,
    DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}
 
void CppUnitOptionTest::CallDecreasesWhenForeignRateIncreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Call );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "call_decrease_when_foreign_rate_decrease", "", FOREIGNRATE, INCREASE,
    DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::PutIncreasesWhenForeignRateIncreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Put );
  m_pOptionInterface->SetStrike( m_pOptionInterface->GetSpotSharePrice() );

  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "put_increase_when_foreign_rate_increase", "", FOREIGNRATE, INCREASE,
    INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::DeltaInTheMoneyEqualToOne()
{
 setup();

 bool bResult = OptionDeltaInTheMoneyEqualToOne(m_pOptionInterface, testParam, root);

 CPPUNIT_ASSERT( bResult );
}
  
void CppUnitOptionTest::DeltaOutOfTheMoneyEqualToZero()
{
  setup();

  bool bResult = OptionDeltaOutOfTheMoneyEqualToZero(m_pOptionInterface, testParam,
    root);

  CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::GammaOutInOfTheMoneyEqualToZero()
{
  setup();

  bool bResult = OptionGammaOutInOfTheMoneyEqualToZero(m_pOptionInterface, testParam,
    root);

 CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::ThetaNegativeRiskFreeRateZero()
{
  setup();

  bool bResult = OptionThetaNegativeRiskFreeRateZero(m_pOptionInterface, testParam,
    root);

 CPPUNIT_ASSERT( bResult );
}

void CppUnitOptionTest::PutDecreasesWhenStrikeDecreases()
{
  setup();

  m_pOptionInterface->SetOptionType( finance::Option_Put );


  bool bResult = CoherenceTest(m_pOptionInterface, root, 
    "put_decrease_when_strike_decrease", "", STRIKE, DECREASE,
    DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

 
void CppUnitOptionTest::CallGoesToSWhenStrikeDecreases()
{
  setup();
  
  bool bResult = OptionCallGoesToSWhenStrikeDecreases(m_pOptionInterface, testParam,
    root);

  CPPUNIT_ASSERT( bResult );
}

//____________________________________________________________________________________//

void CppUnitOptionLegacyTests::Test3()
{
  /*
  std::cout<< "Test 3" << std::endl;
  std::cout<< "  Purpose: Basic European call option" <<std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, zero hr" << std::endl;
  std::cout<< "           non-uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 13.7163503;
  double dDelta = 0.61054344;
  double dGamma = 0.012784239;
  double dVega = 38.3527;
  //double dRho = 0.4734;
  //double dThetaOneDay = -0.0209;
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );
}


void CppUnitOptionLegacyTests::Test4()
{
  /*
  std::cout<< "Test 4" << std::endl;
  std::cout<< "  Purpose: Basic European put option" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, zero hr" <<std::endl;
  std::cout<< "           non-uniform meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile   =  "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 9.870196496;
  double dDelta = -0.389456558;
  double dGamma = 0.012784239;
  double dVega = 38.3527;
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}


void CppUnitOptionLegacyTests::Test5()
{
  /*
  std::cout<< "Test 5" << std::endl;
  std::cout<< "  Purpose: Basic European call option with foreign yield" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, zero dividends, medium vol, zero hr" << std::endl;
  std::cout<< "           non-uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";
  
  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 12.544088388;
  double dDelta = 0.57353911;
  double dGamma = 0.012740207;
  double dVega = 38.2206;
  
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}


void CppUnitOptionLegacyTests::Test6()
{
  /*
  std::cout<< "Test 6" << std::endl;
  std::cout<< "  Purpose: Basic European put option with foreign yield" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, zero dividends, medium vol, zero hr" << std::endl;
  std::cout<< "           non-uniform meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 10.658718856;
  double dDelta = -0.40685304;
  double dGamma = 0.012740207;
  double dVega = 38.2206;
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}

void CppUnitOptionLegacyTests::Test7()
{
  /*
  std::cout<< "Test 7" << std::endl;
  std::cout<< "  Purpose: European call option with foreign yield and 2 dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, two mixed dividends, medium vol," << std::endl;
  std::cout<< "           zero hr, non-uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 9.4709471;
  double dDelta = 0.47360002;
  double dGamma = 0.012719321;
  double dVega = 37.123866;
  
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}


void CppUnitOptionLegacyTests::Test8()
{
  /*
  std::cout<< "Test 8" << std::endl;
  std::cout<< "  Purpose: European put option with foreign yield and two mixed dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, two mixed dividends, medium vol," << std::endl;
  std::cout<< "           zero hr, uniform meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 13.266954;
  double dDelta = -0.45900693;
  double dGamma = 0.012244083;
  double dVega = 36.733276;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);
  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 8.e-4) );


}

void CppUnitOptionLegacyTests::Test9()
{
  /*
  std::cout<< "Test 9" << std::endl;
  std::cout<< "  Purpose: Basic American call option" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, zero hr" << std::endl;
  std::cout<< "           non-uniform meshes, call payoff, American" << std::endl;
  std::cout<< std::endl << std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\amercall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 13.7163503;
  double dDelta = 0.61054344;
  double dGamma = 0.012784239;
  double dVega = 38.3527;
  //double dRho = 0.4734;
  //double dThetaOneDay = -0.0209;
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}


void CppUnitOptionLegacyTests::Test10()
{
  /*
  std::cout<< "Test 10" << std::endl;
  std::cout<< "  Purpose: Basic American put option" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, zero hr" << std::endl;
  std::cout<< "           non-uniform meshes, put payoff, American" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\amerput.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 10.257353;
  double dDelta = -0.41152001;
  double dGamma = 0.014036219;
  double dVega  = 38.35527;

  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}

void CppUnitOptionLegacyTests::Test11()
{
  /*
  std::cout<< "Test 11" << std::endl;
  std::cout<< "  Purpose: American put option with foreign yield and dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, two mixed dividends, medium vol," << std::endl;
  std::cout<< "           zero hr, non-uniform meshes, put payoff, American" << std::endl;
  std::cout<< std::endl;
 */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\amerput.txt";
 

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
    
  double dPrice = 13.548943;
  double dDelta = -0.47219695;
  double dGamma = 0.012754025;
  double dVega = 36.7223;

  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}

void CppUnitOptionLegacyTests::Test12()
{
  /*
  std::cout<< "Test 12" << std::endl;
  std::cout<< "  Purpose: Basic European call option with flat hazard rate" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, medium hr" << std::endl;
  std::cout<< "           non-uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrmedium.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 18.815157;
  double dDelta = 0.73041516;
  double dGamma = 0.011013006;
  double dVega = 33.039085;

  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );
}


void CppUnitOptionLegacyTests::Test13()
{
  /*
  std::cout<< "Test 13" << std::endl;
  std::cout<< "  Purpose: Basic European put option with flat hazard rate" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, medium hr" << std::endl;
  std::cout<< "           non-uniform meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrmedium.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";
  

  OptionTester optionTester;
  
  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 14.969003;
  double dDelta = -0.26958485;
  double dGamma = 0.011013195;
  double dVega  = 33.041; //42.189333;

   CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega) );

}

void CppUnitOptionLegacyTests::Test14()
{
  /*
  std::cout<< "Test 14" << std::endl;
  std::cout<< "  Purpose: Basic European call option with space dependent hazard rate" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, medium hr" << std::endl;
  std::cout<< "           uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";
  

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 23.625661;
  double dDelta = 0.77601324;
  double dGamma = 0.00801165;
  double dVega = 25.326749;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);

  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 6.e-4) );

}

void CppUnitOptionLegacyTests::Test15()
{
  /*
  std::cout<< "Test 15" << std::endl;
  std::cout<< "  Purpose: Basic European put option with space dependent hazard rate" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, medium vol, medium hr" << std::endl;
  std::cout<< "           *uniform* meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 19.779507;
  double dDelta = -0.22398676;
  double dGamma = 0.0080118;
  double dVega = 42.10209;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);

  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 8.e-4) );
}

void CppUnitOptionLegacyTests::Test16()
{
  /*
  std::cout<< "Test 16" << std::endl;
  std::cout<< "  Purpose: Basic European call option with vol surface" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, parametric vol, zero hr" << std::endl;
  std::cout<< "           uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  double dPrice = 8.117522;
  double dDelta = 0.63029509;
  double dGamma = 0.02501627;
  double dVega = 37.45192;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega, 4);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);

  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 1.e-3) );

}


void CppUnitOptionLegacyTests::Test17()
{
  /*
  std::cout<< "Test 17" << std::endl;
  std::cout<< "  Purpose: Basic European put option with vol surface" << std::endl;
  std::cout<< "  Details: medium yield, zero foreign, zero dividends, parametric vol, zero hr" << std::endl;
  std::cout<< "           uniform meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  double dPrice = 4.271367;
  double dDelta = -0.3697050;
  double dGamma = 0.02501645;
  double dVega = 37.45199;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega, 4);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);

  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 1.e-3) );

}


void CppUnitOptionLegacyTests::Test18()
{
  /*
  std::cout<< "Test 18" << std::endl;
  std::cout<< "  Purpose: European call with surface vol, surface hr, mixed dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, mixed dividends, parametric vol," << std::endl;
  std::cout<< "           parametric hr, uniform meshes, call payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 14.075379;
  double dDelta = 0.75419804;
  double dGamma = 0.012779212;
  double dVega = 19.46705;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega, 4);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 1.e-3) );

}


void CppUnitOptionLegacyTests::Test19()
{
  /*
  std::cout<< "Test 19" << std::endl;
  std::cout<< "  Purpose: European put with surface vol, surface hr, mixed dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, mixed dividends, parametric vol," << std::endl;
  std::cout<< "           parametric hr, uniform meshes, put payoff, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\europut.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 17.856460;
  double dDelta = -0.18187765;
  double dGamma = 0.0121788;
  double dVega = 36.08077;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega, 4);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);
 
  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 1.e-3) );

}


void CppUnitOptionLegacyTests::Test20()
{
  /*
  std::cout<< "Test 20" << std::endl;
  std::cout<< "  Purpose: American call with surface vol, surface hr, mixed dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, mixed dividends, parametric vol," << std::endl;
  std::cout<< "           parametric hr, uniform meshes, call payoff, American" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\amercall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 14.075836;
  double dDelta = 0.7543199;
  double dGamma = 0.01280839;
  double dVega = 19.49316;
  optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega, 4);

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 1.e-3) );

}


void CppUnitOptionLegacyTests::Test21()
{
  /*
  std::cout<< "Test 21" << std::endl;
  std::cout<< "  Purpose: American put with surface vol, surface hr, mixed dividends" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, mixed dividends, parametric vol," << std::endl;
  std::cout<< "           parametric hr, uniform meshes, put payoff, American" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\amerput.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(100.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 18.2503;
  double dDelta = -0.18951;
  double dGamma = 0.0129;
  double dVega  = 19.62;
  CPPUNIT_ASSERT( optionTester.RunFullTests(dPrice, dDelta, dGamma, dVega, 4) );

  //std::cout<< "  With non-uniform meshes" << std::endl;

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  optionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( optionTester.RunPriceTests(dPrice, 5.e-4) );
}

void CppUnitOptionLegacyTests::Test22()
{
  /*
  std::cout<< "Test 22" << std::endl;
  std::cout<< "  Purpose: Check put-call parity with surface vol, surface hr" << std::endl;
  std::cout<< "  Details: medium yield, small foreign, zero dividends, parametric vol," << std::endl;
  std::cout<< "           parametric hr, non-uniform meshes, put/call payoffs, European" << std::endl;
  std::cout<< std::endl;
  */

  std::string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  std::string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  std::string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  std::string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  std::string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  std::string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  std::string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  std::string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  std::string sOptionFile        = "txtfiles\\eurocall.txt";

  OptionTester optionTester;

  optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  optionTester.SetSpotSharePrice(80.0);
  
  optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  optionTester.RunConvergenceTest(4);
  double dCallPrice = optionTester.GetPrice();

  sOptionFile = "txtfiles\\europut.txt";
  optionTester.ReadContractParams(sOptionFile);
  optionTester.RunConvergenceTest(4);
  double dPutPrice = optionTester.GetPrice();

  double dR = 0.04;
  double dRf = 0.02;
  double dS0 = 80.0;
  double dK = 100.0;
  //std::cout<< "  Call = " << dCallPrice << ", Put = " << dPutPrice << std::endl;
  double dCMinusP = dCallPrice - dPutPrice;
  //std::cout<< "  Call - Put  = " << dCMinusP << std::endl;
  double r = log(1 + dR);  
  double dOther = dS0 - dK * exp(-r * 1.0) - (dS0*dRf)/(1+dRf);
  //std::cout<< "  S - Ke^(-rT) - PV(foreign) = " << dOther << std::endl;
  bool bResult = optionTester.ReportPass(0.0, dCMinusP - dOther, 1.e-4);

//  std::cout<< std::endl;

  CPPUNIT_ASSERT( bResult );
}
} //end namespace test
}//end namespace ihg
}//end namespace ito33

