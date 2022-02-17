
/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/cppunit_asianoption_test.cpp
// Purpose:     Base class for testing asian option
// Author:      Ito 33 Canada
// Created:     2005/06/13
// RCS-ID:      $Id: cppunit_asianoption_test.cpp,v 1.10 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"
#include "ito33/dateutils.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/exoticoption/curran.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/optiontype.h"

#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/volatilityflat.h"

#include "ihg/xml/pricingreader.h"

#include "cppunit_asianoption_test.h"
#include "testparam.h"
#include "comparisontests.h"
#include "internalcoherencetests.h"

extern ito33::XML::RootTag root;
extern ito33::ihg::test::TestParam testParam;

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33 
{
  

namespace ihg
{

namespace test
{

void CppUnitAsianOptionTest::setup()
{
 
  XML::PricingReader reader( testParam.GetFileName() );
  shared_ptr<finance::SessionData> pSessionData(reader.ReadSessionData());

  finance::DerivativeVisitorGoodType visitor;
  reader.ReadDerivatives(visitor);
  shared_ptr<finance::AsianOption> pAsianOption = visitor.GetAsianOption();
  
  if(!pAsianOption)
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no asian option in input xml file");

  shared_ptr<ito33::ihg::TheoreticalModel> 
      pModel(new ito33::ihg::TheoreticalModel);
  reader.ReadTheoreticalModel(pModel);

  m_pAsianOptionInterface = make_ptr
      ( new AsianOptionInterface(pSessionData, pAsianOption, pModel) );

}
  
void CppUnitAsianOptionTest::PriceAfterMaturityGeneratesAnError()  
{
   setup();

  Date maturityDate(2003, Date::Jan, 1);
  Date valuationDate(2004, Date::Jan, 1);

  m_pAsianOptionInterface->SetMaturityDate( maturityDate );
  m_pAsianOptionInterface->SetValuationDate( valuationDate );

  m_pAsianOptionInterface->Solve();
}


void CppUnitAsianOptionTest::PriceAlwaysPositive()
{
  setup();

  m_pAsianOptionInterface->Solve();

  if ( m_pAsianOptionInterface->GetPrice() < 0 )
    m_pAsianOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
    "price_always_positive");

  CPPUNIT_ASSERT( m_pAsianOptionInterface->GetPrice() >= 0 );
}


void CppUnitAsianOptionTest::GammaAlwaysPositive()
{
  setup();

  m_pAsianOptionInterface->Solve();

  double dGamma = m_pAsianOptionInterface->GetGamma();

  if ( dGamma < 0 && fabs(dGamma) > 1.e-3 )
    m_pAsianOptionInterface->CreateDebugOutputFile(testParam.GetFileName(),
    "price_always_positive");

  if ( fabs(dGamma) < 1.e-3 )
    dGamma = fabs(dGamma);

  CPPUNIT_ASSERT( dGamma >= 0 );
}

void CppUnitAsianOptionTest::PriceIncreasesWhenVolatilityIncreases()
{
  setup();

  if ( m_pAsianOptionInterface->GetStrike() > 0 )
    m_pAsianOptionInterface->SetStrike( m_pAsianOptionInterface->GetSpotSharePrice() );
  else
    return ;

  testParam.m_dVolStep = .2;
  testParam.m_dVolMax  = 1.;

  bool bResult = CoherenceTest(m_pAsianOptionInterface, root, 
    "price_increase_when_volatility_increase", "", VOL, INCREASE, INCREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}
  
void CppUnitAsianOptionTest::PriceDecreasesWhenVolatilityDecreases()
{
  setup();

  if ( m_pAsianOptionInterface->GetStrike() > 0 )
    m_pAsianOptionInterface->SetStrike( m_pAsianOptionInterface->GetSpotSharePrice() );
  else
    return ;

  testParam.m_dVolStep = .2;
  testParam.m_dVolMax  = 1.;

  bool bResult = CoherenceTest(m_pAsianOptionInterface, root, 
    "price_decrease_when_vol_decrease", "", VOL, DECREASE, DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}


void CppUnitAsianOptionTest::AmericanGreaterEuropean()
{
  setup();

  m_pAsianOptionInterface->SetExerciseType( finance::ExerciseType_European);

  m_pAsianOptionInterface->Solve();

  double dPriceEuropean = m_pAsianOptionInterface->GetPrice();
 
  m_pAsianOptionInterface->SetExerciseType( finance::ExerciseType_American);

  m_pAsianOptionInterface->Solve();

  double dPriceAmerican = m_pAsianOptionInterface->GetPrice();

  CPPUNIT_ASSERT( dPriceEuropean <= dPriceAmerican );
}

void CppUnitAsianOptionTest::PriceIncreasesWhenHazardRateIncreases()
{
  setup();

  if ( m_pAsianOptionInterface->GetStrike() > 0 )
  {
    m_pAsianOptionInterface->SetStrike( m_pAsianOptionInterface->GetSpotSharePrice() );
  
    testParam.m_dHazardRateStep = .2;

    bool bResult = CoherenceTest(m_pAsianOptionInterface, root, 
      "price_increase_when_hazardrate_increase", "", HAZARDRATE, 
      INCREASE, INCREASE, testParam);

    CPPUNIT_ASSERT( bResult );
  }
}


void CppUnitAsianOptionTest::PutDecreasesWhenStrikeDecreases()
{
  setup();

  m_pAsianOptionInterface->SetOptionType( finance::Option_Put );

 if ( m_pAsianOptionInterface->GetStrike() < 0 )
  return;

 testParam.m_dStrikeStep = .2;
 testParam.m_dStrikeMax = 1.;

  bool bResult = CoherenceTest(m_pAsianOptionInterface, root, 
    "put_decrease_when_strike_decrease", "", STRIKE, DECREASE,
    DECREASE, testParam);

  CPPUNIT_ASSERT( bResult );
}

///////////////////////////////////////////////////////////////////////////////

void CppUnitCurranAsianOptionTest::CheckCurranAtAveragingPeriod()
{
  double dSpot           = 100.0;
  double dContinuousRate = 0.09;
  double dAnnualRate     = exp(dContinuousRate) - 1.0;

  double ONEWEEK       = 1./52.;
  double dMaturityTime = 52.*ONEWEEK;

  size_t nNbAveraging = 53;

  Date valuationDate   = Date(2006, Date::Jan, 1);
  Date maturityDate    = GetDateFrom( dMaturityTime + GetDoubleFrom( valuationDate ) );
  Date avgStartDate    = GetDateFrom( GetDoubleFrom(valuationDate)  - 1./52);

  finance::ExerciseType exerType   = finance::ExerciseType_European;
    
  shared_ptr<finance::Numeraire> pNumeraire( new finance::Numeraire("EUR"));
  shared_ptr<finance::Equity> pEquity( new finance::Equity(dSpot, pNumeraire) );
  shared_ptr<finance::YieldCurve> pyf(new finance::YieldCurveFlat(0.0));
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<finance::RateData> pRateData( new finance::RateData );
  shared_ptr<finance::YieldCurve> pyc(new finance::YieldCurveFlat(dAnnualRate));
  pRateData->SetYieldCurve(pNumeraire, pyc);

  //At averaging period
  //start average date == valuation date
  std::vector<Curran> atAveragingPeriod;
  atAveragingPeriod.push_back( Curran(.05, 95.,  8.81) );
  atAveragingPeriod.push_back( Curran(.05, 100., 4.31) );
  atAveragingPeriod.push_back( Curran(.05, 105., 0.95) );
  
  atAveragingPeriod.push_back( Curran(.10, 95., 8.91) );
  atAveragingPeriod.push_back( Curran(.10, 100., 4.91) );
  atAveragingPeriod.push_back( Curran(.10, 105., 2.06) );

  atAveragingPeriod.push_back( Curran(.30, 90., 14.96) );
  atAveragingPeriod.push_back( Curran(.30, 100., 8.80) );
  atAveragingPeriod.push_back( Curran(.30, 110., 4.67) );

  atAveragingPeriod.push_back( Curran(.50, 90., 18.14) );
  atAveragingPeriod.push_back( Curran(.50, 100., 12.98) );
  atAveragingPeriod.push_back( Curran(.50, 110., 9.07) );

  size_t nSize = atAveragingPeriod.size();
  size_t nIdx;

  
  shared_ptr<finance::SessionData> 
    pSessionData(new finance::SessionData(pRateData, pEquity, valuationDate));

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++)
  {
   
    double dStrike = atAveragingPeriod[nIdx].m_dStrike;

    shared_ptr<finance::AsianOption> asianOptCall(new 
      finance::AsianOption(dStrike, maturityDate, finance::Option_Call, exerType, 
                 avgStartDate, nNbAveraging));
 
    asianOptCall->SetCurrentAverage( dSpot, 1);

    shared_ptr<finance::AsianOption> asianOptPut(new 
      finance::AsianOption(dStrike, maturityDate, finance::Option_Put, exerType,
                   avgStartDate, nNbAveraging) );

    asianOptPut->SetCurrentAverage( dSpot, 1);

    asianOptCall->SetSessionData(pSessionData);
    asianOptPut->SetSessionData(pSessionData);

    // setup the model
    double dVolatility = atAveragingPeriod[nIdx].m_dVol;
    double dHazardRate = 0.0;
    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
    pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

    // price
    shared_ptr<finance::ModelOutput> outputCall = pModel->Compute(*asianOptCall);
    shared_ptr<finance::ModelOutput> outputPut = pModel->Compute(*asianOptPut);

    double dPriceCall = finance::CurranCall(dSpot, dStrike, 1.0, 53, 
      0.0, dVolatility, dContinuousRate, 0);

     double dPricePut = finance::CurranPut(dSpot, dStrike, 1.0, 53, 
      0.0, dVolatility, dContinuousRate, 0);

     double dErrorCall = fabs(outputCall->GetPrice()/dPriceCall-1.0);
     if ( dPriceCall < 1.0 )
       dErrorCall = fabs(outputCall->GetPrice() - dPriceCall);

     double dErrorPut = fabs(outputPut->GetPrice()/dPricePut-1.0);
     if ( dPricePut < 1.0 )
       dErrorPut = fabs(outputPut->GetPrice() - dPricePut);

     /*
     std::cout << std::endl;
     std::cout << "vol=" << dVolatility << ", K =" << dStrike << std::endl;
     std::cout <<"--->PDE call              =" << outputCall->GetPrice() << std::endl;
     std::cout <<"--->Call curran analytical=" << dPriceCall << std::endl;
     std::cout <<"------>Call relative error="<< dErrorCall <<std::endl;
     std::cout <<"--->PDE put               =" << outputPut->GetPrice()<< std::endl;
     std::cout <<"--->Put curran analytical ="<< dPricePut << std::endl;
     std::cout <<"------>Put relative error="<<  dErrorPut<<std::endl;
     */

     CPPUNIT_ASSERT( dErrorCall < 1.e-3 * dSpot );
     CPPUNIT_ASSERT( dErrorPut < 1.e-3 * dSpot );
  }

}

void CppUnitCurranAsianOptionTest::CheckCurranBeforeAveragingPeriod()
{

  double dSpot           = 100.0;
  double dContinuousRate = 0.09;
  double dAnnualRate     = exp(dContinuousRate) - 1.0;

  double ONEWEEK       = 1./52.;
  double dMaturityTime = 72.*ONEWEEK;
  double dAvgStartTime      = 20.*ONEWEEK;
  double dDividend     = 0.;
  size_t nNbAveraging  = 53;

  Date valuationDate   = Date(2006, Date::Jan, 1);
  Date maturityDate    = GetDateFrom( dMaturityTime + GetDoubleFrom( valuationDate ) );
  Date avgStartDate    = GetDateFrom( dAvgStartTime + GetDoubleFrom( valuationDate )  - 1/52.);

  finance::ExerciseType exerType   = finance::ExerciseType_European;
    
  shared_ptr<finance::Numeraire> pNumeraire( new finance::Numeraire("EUR"));
  shared_ptr<finance::Equity> pEquity( new finance::Equity(dSpot, pNumeraire) );
  shared_ptr<finance::YieldCurve> pyf(new finance::YieldCurveFlat(0.0));
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<finance::YieldCurve> pyc(new finance::YieldCurveFlat(dAnnualRate));
  shared_ptr<finance::RateData> pRateData(new finance::RateData);
  pRateData->SetYieldCurve(pNumeraire, pyc);

  //At averaging period
  //start average date == valuation date
  std::vector<Curran> atAveragingPeriod;
  atAveragingPeriod.push_back( Curran(.05, 95.,  11.76) );
  atAveragingPeriod.push_back( Curran(.05, 100., 7.39) );
  atAveragingPeriod.push_back( Curran(.05, 105., 3.47) );
  
  atAveragingPeriod.push_back( Curran(.10, 95., 11.96) );
  atAveragingPeriod.push_back( Curran(.10, 100., 8.07) );
  atAveragingPeriod.push_back( Curran(.10, 105., 4.87) );

  atAveragingPeriod.push_back( Curran(.30, 90., 19.71) );
  atAveragingPeriod.push_back( Curran(.30, 100., 13.43) );
  atAveragingPeriod.push_back( Curran(.30, 110., 9.05) );

  atAveragingPeriod.push_back( Curran(.50, 90., 24.10) );
  atAveragingPeriod.push_back( Curran(.50, 100., 19.37) );
  atAveragingPeriod.push_back( Curran(.50, 110., 15.47) );

  size_t nSize = atAveragingPeriod.size();
  size_t nIdx;

  
  shared_ptr<finance::SessionData> 
    pSessionData(new finance::SessionData(pRateData, pEquity, valuationDate));

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++)
  {
   
    double dStrike = atAveragingPeriod[nIdx].m_dStrike;

    shared_ptr<finance::AsianOption> asianOptCall(new 
      finance::AsianOption(dStrike, maturityDate, finance::Option_Call, 
                  exerType, avgStartDate, nNbAveraging - 1));

        
    shared_ptr<finance::AsianOption> asianOptPut(new 
      finance::AsianOption(dStrike, maturityDate, finance::Option_Put, 
                  exerType, avgStartDate, nNbAveraging - 1) );

    asianOptCall->SetSessionData(pSessionData);
    asianOptPut->SetSessionData(pSessionData);

    // setup the model
    double dVolatility = atAveragingPeriod[nIdx].m_dVol;
    double dHazardRate = 0.0;
    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
    pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

    // price
    shared_ptr<finance::ModelOutput> outputCall = pModel->Compute(*asianOptCall);
    shared_ptr<finance::ModelOutput> outputPut = pModel->Compute(*asianOptPut);

    double dPriceCall = finance::CurranCall(dSpot, dStrike, dMaturityTime, nNbAveraging,
                 dAvgStartTime, dVolatility, dContinuousRate, dDividend);

     double dPricePut = finance::CurranPut(dSpot, dStrike, dMaturityTime, nNbAveraging,
                 dAvgStartTime, dVolatility, dContinuousRate, dDividend);

     
     double dErrorCall = fabs(outputCall->GetPrice()/dPriceCall-1.0);
     if ( dPriceCall < 1.0 )
       dErrorCall = fabs(outputCall->GetPrice() - dPriceCall);

     double dErrorPut = fabs(outputPut->GetPrice()/dPricePut-1.0);
     if ( dPricePut < 1.0 )
       dErrorPut = fabs(outputPut->GetPrice() - dPricePut);

     
     CPPUNIT_ASSERT( dErrorCall < 1.e-4*dSpot );
     CPPUNIT_ASSERT( dErrorPut < 1.e-4*dSpot );

     /*
     std::cout << std::endl;
     std::cout << "vol=" << dVolatility << ", K =" << dStrike << std::endl;
     std::cout <<"--->PDE call              =" << outputCall->GetPrice() << std::endl;
     std::cout <<"--->Call curran analytical=" << dPriceCall << std::endl;
     std::cout <<"------>Call relative error="<< dErrorCall <<std::endl;
     std::cout <<"--->PDE put               =" << outputPut->GetPrice()<< std::endl;
     std::cout <<"--->Put curran analytical ="<< dPricePut << std::endl;
     std::cout <<"------>Put relative error="<< dErrorPut <<std::endl;
     */

  }

  
}

void CppUnitCurranAsianOptionTest::CheckCurranIntoAveragingPeriod()
{

  double dSpot           = 100.0;
  double dContinuousRate = 0.09;
  double dAnnualRate     = exp(dContinuousRate) - 1.0;

  double ONEWEEK        = 1./52.;
  double dMaturityTime  = 52.* ONEWEEK;
  double dValuationTime = 20.* ONEWEEK;
  double dAvgStartTime  = 0;

  double dDividend        = 0.;
  double dA               = 100.0;
  size_t nNbIntoAveraging = 20;

  size_t nNbAveraging = 53;

  Date avgStartDate    = Date(2006, Date::Jan, 1);
  Date maturityDate    = GetDateFrom( GetDoubleFrom(avgStartDate) + dMaturityTime  + 1./52 );
  Date valuationDate   = GetDateFrom( GetDoubleFrom(avgStartDate) + dValuationTime + 1./52);

  finance::ExerciseType exerType   = finance::ExerciseType_European;
    
  shared_ptr<finance::Numeraire> pNumeraire( new finance::Numeraire("EUR"));
  shared_ptr<finance::Equity> pEquity( new finance::Equity(dSpot, pNumeraire) );
  shared_ptr<finance::YieldCurve> pyf(new finance::YieldCurveFlat(0.0));
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<finance::RateData> pRateData(new finance::RateData);
  shared_ptr<finance::YieldCurve> pyc(new finance::YieldCurveFlat(dAnnualRate));
  pRateData->SetYieldCurve(pNumeraire, pyc);

  //At averaging period
  //start average date == valuation date
  std::vector<Curran> atAveragingPeriod;
  atAveragingPeriod.push_back( Curran(.05, 95.,  6.39) );
  atAveragingPeriod.push_back( Curran(.05, 100., 1.73) );
  atAveragingPeriod.push_back( Curran(.05, 105., 0.01) );
  
  atAveragingPeriod.push_back( Curran(.10, 95., 6.40) );
  atAveragingPeriod.push_back( Curran(.10, 100., 2.11) );
  atAveragingPeriod.push_back( Curran(.10, 105., 0.20) );

  atAveragingPeriod.push_back( Curran(.30, 90., 11.32) );
  atAveragingPeriod.push_back( Curran(.30, 100., 4.13) );
  atAveragingPeriod.push_back( Curran(.30, 110., 0.92) );

  atAveragingPeriod.push_back( Curran(.50, 90., 12.30) );
  atAveragingPeriod.push_back( Curran(.50, 100., 6.24) );
  atAveragingPeriod.push_back( Curran(.50, 110., 2.78) );

  size_t nSize = atAveragingPeriod.size();
  size_t nIdx;

  shared_ptr<finance::SessionData> 
    pSessionData(new finance::SessionData(pRateData, pEquity, valuationDate));

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++)
  {
   
    double dStrike = atAveragingPeriod[nIdx].m_dStrike;

    shared_ptr<finance::AsianOption> asianOptCall(new 
      finance::AsianOption(dStrike, maturityDate, finance::Option_Call, 
                  exerType, avgStartDate, nNbAveraging));

    asianOptCall->SetCurrentAverage(dA, nNbIntoAveraging + 1);
        
    shared_ptr<finance::AsianOption> asianOptPut(new 
      finance::AsianOption(dStrike, maturityDate, finance::Option_Put, 
                  exerType, avgStartDate, nNbAveraging) );

    asianOptPut->SetCurrentAverage(dA, nNbIntoAveraging + 1);

    asianOptCall->SetSessionData(pSessionData);
    asianOptPut->SetSessionData(pSessionData);

    // setup the model
    double dVolatility = atAveragingPeriod[nIdx].m_dVol;
    double dHazardRate = 0.0;
    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
    pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

    // price
    shared_ptr<finance::ModelOutput> outputCall = pModel->Compute(*asianOptCall);
    shared_ptr<finance::ModelOutput> outputPut = pModel->Compute(*asianOptPut);

    double dPriceCall = finance::CurranCall(dSpot, dStrike, 
                 dMaturityTime - dValuationTime, nNbAveraging,
                 dAvgStartTime, dVolatility, dContinuousRate, dDividend, 
                 dA, nNbIntoAveraging);

     double dPricePut = finance::CurranPut(dSpot, dStrike, 
                 dMaturityTime - dValuationTime, nNbAveraging,
                 dAvgStartTime, dVolatility, dContinuousRate, dDividend, 
                 dA, nNbIntoAveraging);
     
     double dErrorCall = fabs(outputCall->GetPrice()/dPriceCall-1.0);
     if ( dPriceCall < 1.0 )
       dErrorCall = fabs(outputCall->GetPrice() - dPriceCall);

     double dErrorPut = fabs(outputPut->GetPrice()/dPricePut-1.0);
     if ( dPricePut < 1.0 )
       dErrorPut = fabs(outputPut->GetPrice() - dPricePut);


     CPPUNIT_ASSERT( dErrorCall < 1.e-3*dSpot );
     CPPUNIT_ASSERT( dErrorPut < 1.e-3*dSpot );

     /*
     std::cout << std::endl;
     std::cout << "vol=" << dVolatility << ", K =" << dStrike << std::endl;
     std::cout <<"--->PDE call              =" << outputCall->GetPrice() << std::endl;
     std::cout <<"--->Call curran analytical=" << dPriceCall << std::endl;
     std::cout <<"------>Call relative error=" << dErrorCall <<std::endl;
     std::cout <<"--->PDE put               =" << outputPut->GetPrice()<< std::endl;
     std::cout <<"--->Put curran analytical =" << dPricePut << std::endl;
     std::cout <<"------>Put relative error=" << dErrorPut <<std::endl;
     */

  }
}

//-----------------------------------------------------------------------------//
  
void CppUnitAsianOptionSpecificTest::ChangeCurrentAverageFixedCall()
{
  //std::cout << "Change current Average for fixed call " << std::endl;
  ChangeCurrentAverage(finance::Option_Call, false);
}

void CppUnitAsianOptionSpecificTest::ChangeCurrentAverageFixedPut()
{
  //std::cout << "Change current Average for fixed Put " << std::endl;
 ChangeCurrentAverage(finance::Option_Put, false);
}

void CppUnitAsianOptionSpecificTest::ChangeCurrentAverageFloatingCall()
{
  //std::cout << "Change current Average for Floating call " << std::endl;
  ChangeCurrentAverage(finance::Option_Call, true);
}

void CppUnitAsianOptionSpecificTest::ChangeCurrentAverageFloatingPut()
{
  //std::cout << "Change current Average for Floating Put " << std::endl;
  ChangeCurrentAverage(finance::Option_Put, true);
}

void CppUnitAsianOptionSpecificTest::ChangeCurrentAverage(
  finance::OptionType optionType, bool bIsFloating)
{
  //std::cout << "Change the current average" << std::endl;
  // Set the session data
  Date valuationDate(2005, Date::Jan, 1);

  // Set the asian option data
  double dStrike = 100.0;

  if ( bIsFloating )
    dStrike = 0.0;

  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);  
  finance::ExerciseType exerType   = finance::ExerciseType_European;
  Date avgStartDate = valuationDate;

  // setup the model
  double dVolatility = .2;
  double dHazardRate = 0.02;
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
  pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

  // Setup the session
  valuationDate.AddMonths(10);
  shared_ptr<finance::SessionData> pSessionData = InitSessionData(valuationDate);
  
  //std::cout << "Issue date: " << issueDate << std::endl;
  //std::cout << "Maturity date: " << maturityDate << std::endl;
  //std::cout << "Valuation date: " << valuationDate << std::endl;
  //std::cout << std::endl;
  
  // loop over current averages
  double dCurrentAverage;

  double dPriceOld = 0.0;

  if ( optionType == finance::Option_Call && bIsFloating == false )
    dPriceOld = 0.0;
  else if  ( optionType == finance::Option_Put && bIsFloating == true )
    dPriceOld = 0.0;
  else if ( optionType == finance::Option_Call && bIsFloating == true ) 
    dPriceOld = 1.e99;
  else if(  optionType == finance::Option_Put && bIsFloating == false )
    dPriceOld = 1.e99;

  for (dCurrentAverage = 50.0; dCurrentAverage < 151.0; dCurrentAverage += 10.0)
  {
      
    // Setup the asian option
    shared_ptr<finance::AsianOption> asianOpt;
    if (dStrike > 0.0)
      asianOpt = make_ptr( new finance::AsianOption
                               (dStrike, maturityDate, optionType, 
                                exerType, avgStartDate, 52) );
    else
      asianOpt = make_ptr( new finance::AsianOption
                               (maturityDate, optionType, 
                                exerType, avgStartDate, 52) );
      
    asianOpt->SetCurrentAverage(dCurrentAverage, 1);

    asianOpt->SetSessionData(pSessionData);

    //pModel->SetDebugOutputFile("./asian.xml");

    // price
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);
    double dPrice = output->GetPrice();

    //std::cout << "current average = " << dCurrentAverage
    //         << ", price = " <<  dPrice << std::endl;

    if ( dPrice > 1.e-2 )
    {
      //fix call (A - K, 0.0)
      if ( optionType == finance::Option_Call && bIsFloating == false)
        CPPUNIT_ASSERT( dPrice >= dPriceOld );
      //floating call ( S - A, 0.0 )
      else if ( optionType == finance::Option_Call && bIsFloating == true )
        CPPUNIT_ASSERT( dPrice <= dPriceOld );
      //fixed put ( K - A , 0.0 )
      else if ( optionType == finance::Option_Put && bIsFloating == false )
        CPPUNIT_ASSERT( dPrice <= dPriceOld );
      //floating put ( A - S, 0.0 )
      else if ( optionType == finance::Option_Put && bIsFloating == true )
        CPPUNIT_ASSERT( dPrice >= dPriceOld );
    }

    dPriceOld = output->GetPrice();
  }
  
  //std::cout << std::endl << std::endl;

}

shared_ptr<finance::SessionData> 
  CppUnitAsianOptionSpecificTest::InitSessionData(Date valuationDate)
{ 

  double dSpot = 100.0;
  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;

  shared_ptr<finance::Numeraire> pNumeraire( new finance::Numeraire("EUR"));

  shared_ptr<finance::Equity> pEquity( new finance::Equity(dSpot, pNumeraire) );

  shared_ptr<finance::YieldCurve> pyf(new finance::YieldCurveFlat(0.0));
  
  pEquity->SetBorrowCurve(pyf);

  //if ( bFakeDividend )
  //{
  //  shared_ptr<finance::Dividends> pDiv = new Dividends();
  //  pDiv->Add(finance::Dividend::Type::Cash, Date(1900,Date::Jun,1),1.);
  //  pEquity->SetDividends(pDiv);
  //}

  shared_ptr<finance::YieldCurve> pyc(new finance::YieldCurveFlat(dAnnualRate));

  shared_ptr<finance::RateData> pRateData( new finance::RateData );
  pRateData->SetYieldCurve(pNumeraire, pyc);
  
  shared_ptr<finance::SessionData> 
    pSessionData(new finance::SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}


} //end namespace test
}//end namespace ihg
}//end namespace ito33

