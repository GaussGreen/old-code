/////////////////////////////////////////////////////////////////////////////
// Name:        tests/hazardrates/main.cpp
// Purpose:     main file for testing hazard rates
// Author:      David Pooley
// Created:     16/06/2004
// RCS-ID:      $Id: testhazardrates.cpp,v 1.11 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <vector>
#include "ito33/afterstd.h"

#include "math.h"

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/array.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"

#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/parametrization.h"

#include "ihg/xml/hazardrate.h"

#include "ito33/tests/utilexml.h"

#include "ito33/xml/write.h"
#include "ihg/xml/spotcomponent.h"

#include "ihg/tests/testhazardrates.h"
#include "ihg/tests/testparametrization.h"


using namespace ito33;
using namespace ito33::ihg;

// ----------------------------------------------------------------------------
// HazardRateFlat tests
// ----------------------------------------------------------------------------

void HazardRateFlatTest::GetValue()
{
  HazardRateFlat hr(0.2);
  
  double dValue = hr.GetValue();
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, 0.2);

}


void HazardRateFlatTest::GetHazardRates()
{
  double dRate = 0.3;
  HazardRateFlat hr(dRate);
  
  // Check the GetHazardRate function for various times
  // and array sizes
  double dValue;
  double dSpot = 100.0;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, dRate);

  hr.GetHazardRates(1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, dRate);

  double pdSpots[2] = {100.0, 200.0};
  double pdValues[2];
  
  hr.GetHazardRates(0.0, pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], dRate);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], dRate);

  hr.GetHazardRates(100.1, pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], dRate);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], dRate);

}

void HazardRateFlatTest::IsTimeOnly()
{
  HazardRateFlat hr(1.0);

  bool bIsTimeOnly = hr.IsTimeOnly();
  CPPUNIT_ASSERT(bIsTimeOnly == true);
}

void HazardRateFlatTest::Dump()
{
  std::ostringstream oss;

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<hazard_rate_flat>\n"
                "<flat>0.2</flat>\n"
                "</hazard_rate_flat>\n"
                "</root>"
              );

  double dRate = 0.2;
  HazardRateFlat hr(dRate);

  ito33::XML::RootTag root("root",oss);

  hr.Dump(root);

}


void HazardRateFlatTest::RateNegative()
{
  // Cannot have negative rate
  HazardRateFlat hr(-0.1);
}

void HazardRateFlatTest::RateTooLarge()
{
  // This one should throw
  HazardRateFlat hr2(10.0 + 1.e-14);

}


// ----------------------------------------------------------------------------
// HazardRateTimeOnly tests
// ----------------------------------------------------------------------------

shared_ptr<HazardRateTimeOnly> HazardRateTimeOnlyTest::CreateHRObject(size_t nSize)
{
  // Create a hazard rate time only class with nSize entries
  Array<double> pdRates(nSize);
  Array<Date> pdDates(nSize);
  for (size_t nIdx = 0; nIdx < nSize; nIdx++)
  {
    pdRates[nIdx] = 0.1 * (1 + nIdx);
    pdDates[nIdx] = Date(2000 + nIdx, Date::Jan, 1);
  }

  return shared_ptr<HazardRateTimeOnly>
    ( new HazardRateTimeOnly(pdDates.Get(), pdRates.Get(), nSize) );
}


void HazardRateTimeOnlyTest::GetHazardRates()
{
  double dValue;
  double dSpot = 100.0;

  // Check the GetHazardRate function for various times
  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(1);
  
  double dTime0 = GetDoubleFrom( Date(2000, Date::Jan, 1) );
  
  hr->GetHazardRates(dTime0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, dValue);

  hr->GetHazardRates(dTime0 - 1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, dValue);

  hr->GetHazardRates(dTime0 + 1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, dValue);

  // Use two times, instead of just one
  hr = CreateHRObject(2);

  double dTime1 = GetDoubleFrom( Date(2001, Date::Jan, 1) );

  hr->GetHazardRates(dTime0 - 1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, dValue);

  hr->GetHazardRates(dTime0 - 1.e-8, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, dValue);

  hr->GetHazardRates(dTime0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

  hr->GetHazardRates( (dTime0 + dTime1)/2.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

  hr->GetHazardRates(dTime1 - 1.e-8, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

  hr->GetHazardRates(dTime1, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

  hr->GetHazardRates(dTime1 + 1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

}


void HazardRateTimeOnlyTest::IsTimeOnly()
{
  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(1);
  bool bIsTimeOnly = hr->IsTimeOnly();
  CPPUNIT_ASSERT(bIsTimeOnly == true);
}


void HazardRateTimeOnlyTest::Dump()
{
   std::ostringstream oss;

  Array<double> pdRates(1);
  Array<Date> pdDates(1);

  pdRates[0] = 0.1;
  pdDates[0] = Date(2001, Date::Jan, 1);
 

 shared_ptr<HazardRateTimeOnly> pHr( 
    new HazardRateTimeOnly(pdDates.Get(), pdRates.Get(), 1) ); 
 
   ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<hazard_rate_time_only>\n"
                "<dates>\n"
                "<date>2001-01-01</date>\n"
                "</dates>\n"
                "<values>\n"
                "<value>0.1</value>\n"
                "</values>\n"
                "</hazard_rate_time_only>\n"
                "</root>"
              );

  ito33::XML::RootTag root("root",oss);

  pHr->Dump(root);

}


void HazardRateTimeOnlyTest::GetValueAtTime()
{
  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(3);
  double dTime0 = GetDoubleFrom( Date(2000, Date::Jan, 1) );
  double dTime1 = GetDoubleFrom( Date(2001, Date::Jan, 1) );
  double dTime2 = GetDoubleFrom( Date(2002, Date::Jan, 1) );

  double dValue;
  
  dValue = hr->GetValueAtTime(dTime0 - 1.e-8);
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, dValue);

  dValue = hr->GetValueAtTime(dTime0);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

  dValue = hr->GetValueAtTime(dTime1 - 1.e-8);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, dValue);

  dValue = hr->GetValueAtTime(dTime1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, dValue);

  dValue = hr->GetValueAtTime(dTime2 - 1.e-8);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, dValue);

  dValue = hr->GetValueAtTime(dTime2);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, dValue);

  dValue = hr->GetValueAtTime(dTime2 + 1.e-8);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, dValue);

}


void HazardRateTimeOnlyTest::GetValues()
{

  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(3);

  std::vector<double> pdValues = hr->GetValues();
  CPPUNIT_ASSERT( pdValues.size() == 3 );
  
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, pdValues[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, pdValues[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, pdValues[2]);

}


void HazardRateTimeOnlyTest::GetTimeComponentValues()
{
  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(3);

  std::vector<double> pdValues = hr->GetTimeComponentValues();
  CPPUNIT_ASSERT( pdValues.size() == 3 );
  
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, pdValues[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, pdValues[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, pdValues[2]);

}


void HazardRateTimeOnlyTest::GetTimes()
{

  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(3);

  std::vector<Date> pdDates = hr->GetDates();
  CPPUNIT_ASSERT( pdDates.size() == 3 );
  
  CPPUNIT_ASSERT( pdDates[0] == Date(2000, Date::Jan, 1) );
  CPPUNIT_ASSERT( pdDates[1] == Date(2001, Date::Jan, 1) );
  CPPUNIT_ASSERT( pdDates[2] == Date(2002, Date::Jan, 1) );

}


void HazardRateTimeOnlyTest::ResetTimeComponent()
{

  shared_ptr<HazardRateTimeOnly> hr = CreateHRObject(3);

  Array<double> pdRatesNew(4);
  Array<Date> pdDatesNew(4);
  for (size_t nIdx = 0; nIdx < 4; nIdx++)
  {
    pdRatesNew[nIdx] = 1.0 + 0.1 * nIdx;
    pdDatesNew[nIdx] = Date(2010 + nIdx, Date::Feb, 2);
  }

  hr->ResetTimeComponent(pdDatesNew.Get(), pdRatesNew.Get(), 4);

  // Check the the change was made
  std::vector<Date> pdDates = hr->GetDates();
  CPPUNIT_ASSERT( pdDates.size() == 4 );
  
  CPPUNIT_ASSERT( pdDates[0] == Date(2010, Date::Feb, 2) );
  CPPUNIT_ASSERT( pdDates[1] == Date(2011, Date::Feb, 2) );
  CPPUNIT_ASSERT( pdDates[2] == Date(2012, Date::Feb, 2) );
  CPPUNIT_ASSERT( pdDates[3] == Date(2013, Date::Feb, 2) );

  std::vector<double> pdValues = hr->GetTimeComponentValues();
  CPPUNIT_ASSERT( pdValues.size() == 4 );
  
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, pdValues[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(1.1, pdValues[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(1.2, pdValues[2]);
  ITO33_ASSERT_DOUBLES_EQUAL(1.3, pdValues[3]);

  // Check if a function can access the new values
  double dTime1 = GetDoubleFrom( pdDates[1] );

  double dValue = hr->GetValueAtTime(dTime1 + 1.e-8);
  ITO33_ASSERT_DOUBLES_EQUAL(1.2, dValue);

  dValue = hr->GetValueAtTime(dTime1 - 1.e-8);
  ITO33_ASSERT_DOUBLES_EQUAL(1.1, dValue);

  dValue = hr->GetValueAtTime(dTime1);
  ITO33_ASSERT_DOUBLES_EQUAL(1.2, dValue);

}


void HazardRateTimeOnlyTest::RateNegative1()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= -0.1;
  pdRates[1]= 0.1;
  pdRates[2]= 0.2;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  HazardRateTimeOnly hr(pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateTimeOnlyTest::RateNegative2()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 0.1;
  pdRates[2]= -0.2;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  HazardRateTimeOnly hr(pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateTimeOnlyTest::RateTooLarge()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 15.0;
  pdRates[2]= 0.4;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  HazardRateTimeOnly hr(pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateTimeOnlyTest::NonincreasingDates1()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 0.2;
  pdRates[2]= 0.3;
  pdDates[0] = Date(2001, Date::Jan, 1);
  pdDates[1] = Date(2000, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  HazardRateTimeOnly hr(pdDates.Get(), pdRates.Get(), 3);

}


void HazardRateTimeOnlyTest::NonincreasingDates2()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 0.2;
  pdRates[2]= 0.3;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2000, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  HazardRateTimeOnly hr(pdDates.Get(), pdRates.Get(), 3);

}


void HazardRateTimeOnlyTest::OutputSerialization()
{
  std::vector<Date> pDates;
  std::vector<double> pdValues;

  pDates.push_back( Date(2005, Date::Jan, 1) );
  pDates.push_back( Date(2005, Date::Jan, 2) );
  pDates.push_back( Date(2005, Date::Jan, 3) );
  pDates.push_back( Date(2005, Date::Jan, 4) );
  pDates.push_back( Date(2005, Date::Jan, 5) );
  pDates.push_back( Date(2005, Date::Jan, 6) );
  pDates.push_back( Date(2005, Date::Jan, 7) );
  pDates.push_back( Date(2005, Date::Jan, 8) );
  pDates.push_back( Date(2005, Date::Jan, 9) );
  pDates.push_back( Date(2005, Date::Jan, 10) );
  pdValues.push_back( .1 );
  pdValues.push_back( .2 );
  pdValues.push_back( .3 );
  pdValues.push_back( .4 );
  pdValues.push_back( .5 );
  pdValues.push_back( .6 );
  pdValues.push_back( .7 );
  pdValues.push_back( .8 );
  pdValues.push_back( .9 );
  pdValues.push_back( 1. );


  ihg::HazardRateTimeOnly hr(pDates, pdValues);

  finance::ModelParametersConsumerTest visitor;

  hr.GetModelParameters(visitor);

  CPPUNIT_ASSERT( visitor.m_pdVal.size() == 10 );

  CPPUNIT_ASSERT( visitor.m_sCategoryName == MODEL_PARAM_NAME_HR_TIMEONLY);

  size_t nIdx;

  for ( nIdx = 0 ; nIdx < visitor.m_pDates.size(); nIdx++ )
  {
    CPPUNIT_ASSERT( visitor.m_pDates[nIdx] == Date(2005, Date::Jan, nIdx+1 ) );
    ITO33_ASSERT_DOUBLES_EQUAL( visitor.m_pdVal[nIdx], (nIdx+1)/10. );
  }

}

// ----------------------------------------------------------------------------
// HazardRatePower tests
// ----------------------------------------------------------------------------

void HazardRatePowerTest::GetHazardRates()
{
  m_dAlpha = 0.3;
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);
  
  // Check the GetHazardRate function for various times
  // and array sizes
  double dValue;
  double dSpot = 100.0;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, m_dAlpha);

  hr.GetHazardRates(1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, m_dAlpha);

  dSpot = 90.0;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(dSpot) );

  double pdSpots[2] = {80.0, 120.0};
  double pdValues[2];
  
  hr.GetHazardRates(0.0, pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL( pdValues[0], Eval(pdSpots[0]) );
  ITO33_ASSERT_DOUBLES_EQUAL( pdValues[1], Eval(pdSpots[1]) );

  hr.GetHazardRates(100.1, pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], Eval(pdSpots[0]) );
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], Eval(pdSpots[1]) );

}

void HazardRatePowerTest::IsTimeOnly()
{
  m_dAlpha = 0.3;
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);

  bool bIsTimeOnly = hr.IsTimeOnly();
  CPPUNIT_ASSERT(bIsTimeOnly == false);

  m_dBeta = 0.0;
  HazardRatePower hr2(m_dAlpha, m_dBeta, m_dS0);

  bIsTimeOnly = hr2.IsTimeOnly();
  CPPUNIT_ASSERT(bIsTimeOnly == true);

}

void HazardRatePowerTest::Dump()
{
  std::ostringstream oss;

  m_dAlpha = 0.3;
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<hazard_rate_power>\n"
                "<alpha>0.3</alpha>\n"
                "<beta>0.5</beta>\n"
                "<S0>100</S0>\n"
                "</hazard_rate_power>\n"
                "</root>"
              );

  ito33::XML::RootTag root("root",oss);

 
  hr.Dump(root);

}


void HazardRatePowerTest::AvoidZeroDivide()
{
  m_dAlpha = 0.3;
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);
  
  // Check that the function does not divide by zero.
  // Anything less than 1.e-16 is set to 1.e-16
  double dValue;
  double dSpot = 0.0;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-16) );

  dSpot = 1.e-32;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-16) );

  dSpot = 1.e-16;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-16) );

  dSpot = 1.e-8;
  hr.GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-8) );

}


void HazardRatePowerTest::AlphaNegative()
{
  m_dAlpha = -0.3;
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);
}


void HazardRatePowerTest::S0Negative()
{
  m_dAlpha = 0.3;
  m_dBeta = 0.5;
  m_dS0 = -100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);
}


void HazardRatePowerTest::BetaNegative()
{
  m_dAlpha = 0.3;
  m_dBeta = -0.5;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);
}

void HazardRatePowerTest::BetaTooLarge()
{
  m_dAlpha = 0.3;
  m_dBeta = 2.1;
  m_dS0 = 100.0;
  HazardRatePower hr(m_dAlpha, m_dBeta, m_dS0);
}

void HazardRatePowerTest::SerializationOutput()
{
   double dAlpha        = .2;
  double dBeta         = 1.5;
  double dS0           = 52.3;

  HazardRatePower hr(dAlpha, dBeta, dS0);

  finance::ModelParametersConsumerTest visitor;

  hr.GetModelParameters(visitor);


  CPPUNIT_ASSERT( visitor.m_pParameters.size() == 3);

  CPPUNIT_ASSERT( visitor.m_sCategoryName == MODEL_PARAM_NAME_HR_POWER);

  CPPUNIT_ASSERT( visitor.m_pParameters[0].name == SCALAR_MODEL_PARAM_NAME_ALPHA);
  CPPUNIT_ASSERT( visitor.m_pParameters[1].name == SCALAR_MODEL_PARAM_NAME_BETA); 
  CPPUNIT_ASSERT( visitor.m_pParameters[2].name == SCALAR_MODEL_PARAM_NAME_S0);

  CPPUNIT_ASSERT( visitor.m_pParameters[0].value == .2 );
  CPPUNIT_ASSERT( visitor.m_pParameters[1].value == 1.5 );
  CPPUNIT_ASSERT( visitor.m_pParameters[2].value == 52.3 );
}
// ----------------------------------------------------------------------------
// HazardRateCombo tests
// ----------------------------------------------------------------------------

shared_ptr<HazardRateCombo> HazardRateComboTest::CreateHRObject(size_t nSize)
{
  // Create a hazard rate combo class with nSize entries, and combined
  // using the specified operator

  // The return structure
  shared_ptr<HazardRateCombo> pHRCombo;

  // The spot component
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  if (nSize == 0)
  {
    pHRCombo = make_ptr( new HazardRateCombo(pSpotComponent) );
  }
  else
  {

    Array<double> pdRates(nSize);
    Array<Date> pdDates(nSize);
    for (size_t nIdx = 0; nIdx < nSize; nIdx++)
    {
      pdRates[nIdx] = 0.1 * (1 + nIdx);
      pdDates[nIdx] = Date(2000 + nIdx, Date::Jan, 1);
    }
    pHRCombo = make_ptr( new HazardRateCombo
        (pSpotComponent, pdDates.Get(), pdRates.Get(), nSize) );
  }

    
  return pHRCombo;
}


void HazardRateComboTest::GetHazardRatesCtor1()
{
  double dValue;
  double dSpot = 100.0;

  // If there is no time component, then it should be multiplying by zero
  shared_ptr<HazardRateCombo> hr = CreateHRObject(0);
  
  hr->GetHazardRates(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, dValue);

  hr->GetHazardRates(1000.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, dValue);

}

void HazardRateComboTest::GetHazardRatesCtor2()
{
  double dValue;
  double dSpot = 100.0;

  double dTime0 = GetDoubleFrom( Date(2000, Date::Jan, 1) );
  double dTime1 = GetDoubleFrom( Date(2001, Date::Jan, 1) );
  double dTime2 = GetDoubleFrom( Date(2002, Date::Jan, 1) );

  // Now test multiplication.  
  shared_ptr<HazardRateCombo> hr = CreateHRObject(3);
  
  hr->GetHazardRates(dTime0 - 1.e-4, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(Eval(dSpot) * 0.1, dValue);

  hr->GetHazardRates(dTime1 - 1.e-4, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(Eval(dSpot) * 0.2, dValue);

  hr->GetHazardRates(dTime2 - 1.e-4, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(Eval(dSpot) * 0.3, dValue);

  hr->GetHazardRates(dTime2 + 1.e-4, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(Eval(dSpot) * 0.3, dValue);

}


void HazardRateComboTest::IsTimeOnly()
{
  shared_ptr<HazardRateCombo> hr = CreateHRObject(1);
  bool bIsTimeOnly = hr->IsTimeOnly();
  CPPUNIT_ASSERT(bIsTimeOnly == false);
}


void HazardRateComboTest::Dump()
{
  shared_ptr<HazardRateCombo> hr = CreateHRObject(1);

  
  std::ostringstream oss;

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<hazard_rate_combo>\n"
                "<spot_component>\n"
                "<spot_component_power>\n"
                  "<beta>0.5</beta>\n"
                  "<S0>100</S0>\n"
                "</spot_component_power>\n"
                "</spot_component>\n"
                "<time_component>\n"
                  "<dates>\n"
                  "<date>2000-01-01</date>\n"
                  "</dates>\n"
                  "<values>\n"
                  "<value>0.1</value>\n"
                  "</values>\n"
                "</time_component>\n"
                "</hazard_rate_combo>\n"
                "</root>"
              );

 
  ito33::XML::RootTag root("root",oss);

 
  hr->Dump(root);
  
}



void HazardRateComboTest::GetTimeComponentValues()
{
  shared_ptr<HazardRateCombo> hr = CreateHRObject(3);

  std::vector<double> pdValues = hr->GetTimeComponentValues();
  CPPUNIT_ASSERT( pdValues.size() == 3 );
  
  ITO33_ASSERT_DOUBLES_EQUAL(0.1, pdValues[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.2, pdValues[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.3, pdValues[2]);

}


void HazardRateComboTest::GetTimes()
{

  shared_ptr<HazardRateCombo> hr = CreateHRObject(3);

  std::vector<Date> pdDates = hr->GetDates();
  CPPUNIT_ASSERT( pdDates.size() == 3 );
  
  CPPUNIT_ASSERT( pdDates[0] == Date(2000, Date::Jan, 1) );
  CPPUNIT_ASSERT( pdDates[1] == Date(2001, Date::Jan, 1) );
  CPPUNIT_ASSERT( pdDates[2] == Date(2002, Date::Jan, 1) );

}


void HazardRateComboTest::ResetTimeComponent()
{

  shared_ptr<HazardRateCombo> hr = CreateHRObject(3);

  Array<double> pdRatesNew(4);
  Array<Date> pdDatesNew(4);
  for (size_t nIdx = 0; nIdx < 4; nIdx++)
  {
    pdRatesNew[nIdx] = 1.0 + 0.1 * nIdx;
    pdDatesNew[nIdx] = Date(2010 + nIdx, Date::Feb, 2);
  }

  hr->ResetTimeComponent(pdDatesNew.Get(), pdRatesNew.Get(), 4);

  // Check that the change was made
  std::vector<Date> pdDates = hr->GetDates();
  CPPUNIT_ASSERT( pdDates.size() == 4 );
  
  CPPUNIT_ASSERT( pdDates[0] == Date(2010, Date::Feb, 2) );
  CPPUNIT_ASSERT( pdDates[1] == Date(2011, Date::Feb, 2) );
  CPPUNIT_ASSERT( pdDates[2] == Date(2012, Date::Feb, 2) );
  CPPUNIT_ASSERT( pdDates[3] == Date(2013, Date::Feb, 2) );

  std::vector<double> pdValues = hr->GetTimeComponentValues();
  CPPUNIT_ASSERT( pdValues.size() == 4 );
  
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, pdValues[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(1.1, pdValues[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(1.2, pdValues[2]);
  ITO33_ASSERT_DOUBLES_EQUAL(1.3, pdValues[3]);

  // Check if a function can access the new values
  double dSpot = 90.0;
  double dValue;
  double dTime1 = GetDoubleFrom( pdDates[1] );
  double dTime2 = GetDoubleFrom( pdDates[2] );

  hr->GetHazardRates(dTime1 - 1.e-4, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(Eval(dSpot) * 1.1, dValue);

  hr->GetHazardRates(dTime2 - 1.e-4, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(Eval(dSpot) * 1.2, dValue);

}


void HazardRateComboTest::RateNegative1()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= -0.1;
  pdRates[1]= 0.1;
  pdRates[2]= 0.2;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);

  m_dBeta = 0.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  HazardRateCombo hr(pSpotComponent, pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateComboTest::RateNegative2()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 0.1;
  pdRates[2]= -0.2;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  HazardRateCombo hr(pSpotComponent, pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateComboTest::NonincreasingDates1()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 0.2;
  pdRates[2]= 0.3;
  pdDates[0] = Date(2001, Date::Jan, 1);
  pdDates[1] = Date(2000, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);
  
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  HazardRateCombo hr(pSpotComponent, pdDates.Get(), pdRates.Get(), 3);

}


void HazardRateComboTest::NonincreasingDates2()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 0.2;
  pdRates[2]= 0.3;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2000, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);

  m_dBeta = 0.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  HazardRateCombo hr(pSpotComponent, pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateComboTest::RateTooLarge()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= 12.;
  pdRates[2]= 0.2;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);

  m_dBeta = 0.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  HazardRateCombo hr(pSpotComponent, pdDates.Get(), pdRates.Get(), 3);

}

void HazardRateComboTest::BetaTooLarge()
{
  // Cannot have negative rate
  Array<double> pdRates(3);
  Array<Date> pdDates(3);
  pdRates[0]= 0.1;
  pdRates[1]= .2;
  pdRates[2]= 0.2;
  pdDates[0] = Date(2000, Date::Jan, 1);
  pdDates[1] = Date(2001, Date::Jan, 1);
  pdDates[2] = Date(2002, Date::Jan, 1);

  m_dBeta = 2.5;
  m_dS0 = 100.0;
  shared_ptr<SpotComponent> pSpotComponent
  ( new HRSpotComponentPower(m_dBeta, m_dS0) );

  HazardRateCombo hr(pSpotComponent, pdDates.Get(), pdRates.Get(), 3);

}

// ----------------------------------------------------------------------------
// HRSpotComponentPowerTest tests
// ----------------------------------------------------------------------------

void HRSpotComponentPowerTest::GetHazardRates()
{
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HRSpotComponentPower hr(m_dBeta, m_dS0);
  
 
  // Check the GetHazardRate function for various times
  // and array sizes
  double dValue;
  double dSpot = 100.0;
  hr.GetValues(&dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, 1);

  hr.GetValues( &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, 1);

  dSpot = 90.0;
  hr.GetValues( &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(dSpot) );

  double pdSpots[2] = {80.0, 120.0};
  double pdValues[2];
  
  hr.GetValues( pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL( pdValues[0], Eval(pdSpots[0]) );
  ITO33_ASSERT_DOUBLES_EQUAL( pdValues[1], Eval(pdSpots[1]) );

  hr.GetValues( pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], Eval(pdSpots[0]) );
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], Eval(pdSpots[1]) );

}

void HRSpotComponentPowerTest::Dump()
{
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HRSpotComponentPower hr(m_dBeta, m_dS0);

   std::ostringstream oss;

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<spot_component>\n"
                "<spot_component_power>\n"
                  "<beta>0.5</beta>\n"
                  "<S0>100</S0>\n"
                "</spot_component_power>\n"
                "</spot_component>\n"
                "</root>"
              );

 
  ito33::XML::RootTag root("root",oss);

 
  hr.Dump(XML_TAG_SPOTCOMPONENT_ROOT, root);

}


void HRSpotComponentPowerTest::AvoidZeroDivide()
{
  m_dBeta = 0.5;
  m_dS0 = 100.0;
  HRSpotComponentPower hr(m_dBeta, m_dS0);
  
  // Check that the function does not divide by zero.
  // Anything less than 1.e-16 is set to 1.e-16
  double dValue;
  double dSpot = 0.0;
  hr.GetValues( &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-16) );

  dSpot = 1.e-32;
  hr.GetValues(  &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-16) );

  dSpot = 1.e-16;
  hr.GetValues(  &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-16) );

  dSpot = 1.e-8;
  hr.GetValues(  &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, Eval(1.e-8) );

}

void HRSpotComponentPowerTest::S0Negative()
{
  // Cannot have negative rate
  m_dBeta = 0.5;
  m_dS0 = -100.0;
  HRSpotComponentPower hr(m_dBeta, m_dS0);
}


void HRSpotComponentPowerTest::BetaNegative()
{
  // Cannot have hazard rate increase with spot
  m_dBeta = -0.5;
  m_dS0 = 100.0;
  HRSpotComponentPower hr(m_dBeta, m_dS0);
}

void HRSpotComponentPowerTest::BetaTooLarge()
{
  // Cannot have hazard rate increase with spot
  m_dBeta = 2.1;
  m_dS0 = 100.0;
  HRSpotComponentPower hr(m_dBeta, m_dS0);
}
