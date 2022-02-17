/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testhazardrates.h
// Purpose:     header file for testing hazard rates
// Author:      David Pooley
// Created:     16/06/2004
// RCS-ID:      $Id: testhazardrates.h,v 1.12 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
#include <math.h>

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/common.h"
#include "ito33/sharedptr.h"

namespace ito33
{

namespace ihg
{
  class HazardRateTimeOnly;
  class HazardRateCombo;
}

}


// ----------------------------------------------------------------------------
// HazardRateFlat tests
// ----------------------------------------------------------------------------

class HazardRateFlatTest : public CppUnit::TestCase
{
public:
  HazardRateFlatTest() { }

  virtual void tearDown() 
  {
    m_oss.str().clear();
  }

private:
  CPPUNIT_TEST_SUITE( HazardRateFlatTest );
    CPPUNIT_TEST( GetValue );
    CPPUNIT_TEST( GetHazardRates );
    CPPUNIT_TEST( IsTimeOnly );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST_EXCEPTION( RateNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( RateTooLarge, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void GetValue();
  void GetHazardRates();
  void IsTimeOnly();
  void Dump();

  // Test paramater validity
  void RateNegative();
  void RateTooLarge();

  std::ostringstream m_oss;
  std::string m_expected;

  NO_COPY_CLASS(HazardRateFlatTest);
};


// ----------------------------------------------------------------------------
// HazardRateTimeOnly tests
// ----------------------------------------------------------------------------

class HazardRateTimeOnlyTest : public CppUnit::TestCase
{
public:
  HazardRateTimeOnlyTest() { }

  virtual void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( HazardRateTimeOnlyTest );
    CPPUNIT_TEST( GetHazardRates );
    CPPUNIT_TEST( IsTimeOnly );    
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( GetValueAtTime );
    CPPUNIT_TEST( GetValues );
    CPPUNIT_TEST( GetTimeComponentValues );
    CPPUNIT_TEST( GetTimes );
    CPPUNIT_TEST( ResetTimeComponent );
    CPPUNIT_TEST( OutputSerialization );
    CPPUNIT_TEST_EXCEPTION( RateNegative1, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( RateNegative2, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( RateTooLarge, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NonincreasingDates1, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NonincreasingDates2, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void GetHazardRates();
  void IsTimeOnly();
  void Dump();
  void GetValueAtTime();
  void GetValues();
  void GetTimeComponentValues();
  void GetTimes();
  void ResetTimeComponent();
  void OutputSerialization();

  // Test paramater validity
  void RateNegative1();
  void RateNegative2();
  void NonincreasingDates1();
  void NonincreasingDates2();
  void RateTooLarge();


  // Helper functions
  ito33::shared_ptr<ito33::ihg::HazardRateTimeOnly> CreateHRObject(size_t nSize);

  NO_COPY_CLASS(HazardRateTimeOnlyTest);
};


// ----------------------------------------------------------------------------
// HazardRatePower tests
// ----------------------------------------------------------------------------

class HazardRatePowerTest : public CppUnit::TestCase
{
public:
  HazardRatePowerTest() { }

  virtual void tearDown() {}


private:
  CPPUNIT_TEST_SUITE( HazardRatePowerTest );
    CPPUNIT_TEST( GetHazardRates );
    CPPUNIT_TEST( IsTimeOnly );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( AvoidZeroDivide );
    CPPUNIT_TEST( SerializationOutput );
    CPPUNIT_TEST_EXCEPTION( AlphaNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( S0Negative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( BetaNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( BetaTooLarge, ito33::Exception);
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void GetHazardRates();
  void IsTimeOnly();
  void Dump();
  void AvoidZeroDivide();
  void SerializationOutput();
  
  // Test paramater validity
  void AlphaNegative();
  void S0Negative();
  void BetaNegative();
  void BetaTooLarge();

  // Helper function and parameters
  double m_dAlpha;
  double m_dBeta;
  double m_dS0;

  double Eval(double dS)
  {
    return m_dAlpha * pow(m_dS0/dS, m_dBeta);
  }


  NO_COPY_CLASS(HazardRatePowerTest);
};


// ----------------------------------------------------------------------------
// HazardRateCombo tests
// ----------------------------------------------------------------------------

class HazardRateComboTest : public CppUnit::TestCase
{
public:
  HazardRateComboTest() { }

  virtual void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( HazardRateComboTest );
    CPPUNIT_TEST( GetHazardRatesCtor1 );
    CPPUNIT_TEST( GetHazardRatesCtor2 );
    CPPUNIT_TEST( IsTimeOnly );    
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( GetTimeComponentValues );
    CPPUNIT_TEST( GetTimes );
    CPPUNIT_TEST( ResetTimeComponent );
    CPPUNIT_TEST_EXCEPTION( RateNegative1, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( RateNegative2, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NonincreasingDates1, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NonincreasingDates2, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( RateTooLarge, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( BetaTooLarge, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void GetHazardRatesCtor1();
  void GetHazardRatesCtor2();
  void IsTimeOnly();
  void Dump();
  void GetTimeComponentValues();
  void GetTimes();
  void ResetTimeComponent();

  // Test paramater validity
  void RateNegative1();
  void RateNegative2();
  void NonincreasingDates1();
  void NonincreasingDates2();
  void RateTooLarge();
  void BetaTooLarge();

  // Helper functions and variables
  ito33::shared_ptr<ito33::ihg::HazardRateCombo> CreateHRObject(size_t nSize);

  double Eval(double dS)
  {
    return pow(m_dS0 / dS, m_dBeta);
  }

  double m_dBeta;
  double m_dS0;

  NO_COPY_CLASS(HazardRateComboTest);
};


// ----------------------------------------------------------------------------
// HRSpotComponentPowerTest tests
// ----------------------------------------------------------------------------

class HRSpotComponentPowerTest : public CppUnit::TestCase
{
public:
  HRSpotComponentPowerTest() { }

  virtual void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( HRSpotComponentPowerTest );
    CPPUNIT_TEST( GetHazardRates );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( AvoidZeroDivide );
    CPPUNIT_TEST_EXCEPTION( S0Negative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( BetaNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( BetaTooLarge, ito33::Exception);
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void GetHazardRates();
  void Dump();
  void AvoidZeroDivide();
  
  // Test paramater validity
  void S0Negative();
  void BetaNegative();
  void BetaTooLarge();

  // Helper function and parameters
  double m_dBeta;
  double m_dS0;

  double Eval(double dS)
  {
    return pow(m_dS0/dS, m_dBeta);
  }


  NO_COPY_CLASS(HRSpotComponentPowerTest);
};
