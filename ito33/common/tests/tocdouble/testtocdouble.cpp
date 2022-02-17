/////////////////////////////////////////////////////////////////////////////
// Name:        test/tocdouble/main.cpp
// Purpose:     test file of ToCDouble/FromCDouble test program
// Author:      Vaclav Slavik
// Created:     2005-11-28
// RCS-ID:      $Id: testtocdouble.cpp,v 1.1 2005/11/29 15:39:53 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"

#include "ito33/cppunit.h"

#include <limits>

#ifdef _MSC_VER
  #include <float.h>
#else
  #include <math.h>
#endif

using namespace ito33;
using namespace std;

#ifdef _MSC_VER
  #define isnan(x)  _isnan(x)
  #define isinf(x)  (_fpclass(x) == _FPCLASS_PINF)
#endif


// ----------------------------------------------------------------------------
// test class
// ----------------------------------------------------------------------------

class ToCDoubleTestCase : public CppUnit::TestCase
{
public:
  ToCDoubleTestCase() {}

private:
  CPPUNIT_TEST_SUITE( ToCDoubleTestCase );
    CPPUNIT_TEST( Numbers );
    CPPUNIT_TEST( NaN );
    CPPUNIT_TEST( Infinity );
  CPPUNIT_TEST_SUITE_END();

  void Numbers();
  void NaN();
  void Infinity();

  void TestToStr(double value, const std::string& str);
  void TestFromStr(double value, const std::string& str);
  void TestNumber(double value, const std::string& str);

  NO_COPY_CLASS(ToCDoubleTestCase);
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

void ToCDoubleTestCase::TestToStr(double value, const std::string& str)
{
  // with fixed precision:
  CPPUNIT_ASSERT( str == String::FromCDouble(value, 15) );
  // with unlimited precision:
  CPPUNIT_ASSERT( str == String::FromCDouble(value, -1) );
}

void ToCDoubleTestCase::TestFromStr(double value, const std::string& str)
{
  double value2;
  CPPUNIT_ASSERT( String::ToCDouble(str, &value2) );
  CPPUNIT_ASSERT( value2 == value );
}

void ToCDoubleTestCase::TestNumber(double value, const std::string& str)
{
  TestToStr(value, str);
  TestFromStr(value, str);
}

#define TEST_NUMBER(x) TestNumber(x, #x)

void ToCDoubleTestCase::Numbers()
{
  TEST_NUMBER(0);
  TEST_NUMBER(1);
  TEST_NUMBER(12.34);
  TEST_NUMBER(12005);
  TEST_NUMBER(12);
  TEST_NUMBER(0.33);
  TEST_NUMBER(0.0036);
  TEST_NUMBER(188.2);
  TEST_NUMBER(188.2345);
  TEST_NUMBER(-0.00000000001);
}

void ToCDoubleTestCase::NaN()
{
  double value;
  CPPUNIT_ASSERT( String::ToCDouble("NaN", &value) );
  CPPUNIT_ASSERT( isnan(value) );

  TestToStr(numeric_limits<double>::quiet_NaN(), "NaN");
}

void ToCDoubleTestCase::Infinity()
{
  double value;
  CPPUNIT_ASSERT( String::ToCDouble("Infinity", &value) );
  CPPUNIT_ASSERT( isinf(value) );

  TestToStr(numeric_limits<double>::infinity(), "Infinity");
}

CPPUNIT_TEST_SUITE_REGISTRATION(ToCDoubleTestCase);
