/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testvolatility.h
// Purpose:     testing volatility
// Author:      Yann d'Halluin
// Created:     11/07/2004
// RCS-ID:      $Id: testvolatility.h,v 1.8 2006/01/04 16:12:48 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/common.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace ihg
{
  class VolatilityFlat;
}
}

// ----------------------------------------------------------------------------
// VolatilityFlat tests
// ----------------------------------------------------------------------------

class VolatilityFlatTest : public CppUnit::TestCase
{
public:
  VolatilityFlatTest() { }

  virtual void tearDown() 
  {
    m_oss.str().clear();
  }

private:
  CPPUNIT_TEST_SUITE( VolatilityFlatTest );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( Value );
    CPPUNIT_TEST( GetVolsSquared );
    CPPUNIT_TEST( GetVols );
    CPPUNIT_TEST( SerializationOutput );
    CPPUNIT_TEST_EXCEPTION( VolatilityNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( VolatilityTooLarge, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void GetVolsSquared();
  void GetVols();
  void Value();
  void Dump();

  // Test parameter validity
  void VolatilityNegative();
  void VolatilityTooLarge();

  //test c++ output Serialization
  void SerializationOutput();

  std::ostringstream m_oss;
  std::string m_expected;

  NO_COPY_CLASS(VolatilityFlatTest);
};


class VolatilityPowerTest : public CppUnit::TestCase
{
public:
  VolatilityPowerTest() { }

  virtual void tearDown() 
  {
    m_oss.str().clear();
  }

private:
  CPPUNIT_TEST_SUITE( VolatilityPowerTest );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( GetVolsSquared );
    CPPUNIT_TEST( GetVols );
    CPPUNIT_TEST( SerializationOutput );
    CPPUNIT_TEST_EXCEPTION( AlphaNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( SpotNegative, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void Dump();
  void GetVolsSquared();
  void GetVols();

  // Test parameter validity
  void AlphaNegative();
  void SpotNegative();

  //test c++ output Serialization
  void SerializationOutput();

  std::ostringstream m_oss;
  std::string m_expected;

  NO_COPY_CLASS(VolatilityPowerTest);
};



class VolatilityTanhTest : public CppUnit::TestCase
{
public:
  VolatilityTanhTest() { }

  virtual void tearDown() 
  {
    m_oss.str().clear();
  }

private:
  CPPUNIT_TEST_SUITE( VolatilityTanhTest );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST( GetVolsSquared );
    CPPUNIT_TEST( GetVols );
    CPPUNIT_TEST( SerializationOutput );
    CPPUNIT_TEST_EXCEPTION( LeftNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( RightNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( ScaleNegative, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( SpotNegative, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void Dump();
  void GetVolsSquared();
  void GetVols();

  // Test parameter validity
  void LeftNegative();
  void RightNegative();
  void ScaleNegative();
  void SpotNegative();

  //test c++ output Serialization
  void SerializationOutput();

  std::ostringstream m_oss;
  std::string m_expected;

  NO_COPY_CLASS(VolatilityTanhTest);
};


class VolatilityTimeOnlyTest : public CppUnit::TestCase
{
public:
  VolatilityTimeOnlyTest() { }

  virtual void tearDown() 
  {
    m_oss.str().clear();
  }

private:
  CPPUNIT_TEST_SUITE( VolatilityTimeOnlyTest );
    CPPUNIT_TEST( SerializationOutput );
    CPPUNIT_TEST( Dump );
    CPPUNIT_TEST_EXCEPTION( VolatilityTooLarge, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  //xml dump
  void Dump();

  // Test parameter validity
  void VolatilityTooLarge();

  //test c++ output Serialization
  void SerializationOutput();

  std::ostringstream m_oss;
  std::string m_expected;

  NO_COPY_CLASS(VolatilityTimeOnlyTest);
};