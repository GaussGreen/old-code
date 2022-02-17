/////////////////////////////////////////////////////////////////////////////
// Name:        test/testycurve.h
// Purpose:     header file of yield curve test program
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: testycurve.h,v 1.5 2005/04/27 17:42:44 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_TESTYCURVE_H_
#define _ITO33_TEST_TESTYCURVE_H_

#include "ito33/common.h"
#include "ito33/date.h"
#include "ito33/exception.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/yieldcurve_swap.h"
#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// test classes
// ----------------------------------------------------------------------------

class YCurveFlatTestCase : public CppUnit::TestCase
{
public:

  YCurveFlatTestCase() { }

  void tearDown() { }

private:

  CPPUNIT_TEST_SUITE( YCurveFlatTestCase );

    CPPUNIT_TEST( GetFlatRate );
    CPPUNIT_TEST_EXCEPTION( NegativeRate, ito33::Exception );

  CPPUNIT_TEST_SUITE_END();

  void GetFlatRate();

  void NegativeRate();

  NO_COPY_CLASS(YCurveFlatTestCase);
};


class YCurveAnnuallyCompoundedTestCase : public CppUnit::TestCase
{
public:
  YCurveAnnuallyCompoundedTestCase() : m_yc(ito33::Date(2004, ito33::Date::Jan, 1)) { }

  void tearDown() { m_yc.Clear(); }

private:
  CPPUNIT_TEST_SUITE( YCurveAnnuallyCompoundedTestCase );
    CPPUNIT_TEST( AddLegsReverse );

    CPPUNIT_TEST_EXCEPTION( NoLeg, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( SingleLeg, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( AddLegDuplicate, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( AddInvalidLeg, ito33::Exception );

    CPPUNIT_TEST( Perturb );
    CPPUNIT_TEST( Clear );
    CPPUNIT_TEST( GetLegs );

  CPPUNIT_TEST_SUITE_END();

  void AddLegsReverse();

  void NoLeg();
  void SingleLeg();
  void AddLegDuplicate();
  void AddInvalidLeg();
  void Clear();
  void SetFlatRate();
  void GetLegs();

  void Perturb();

  ito33::finance::YieldCurveAnnuallyCompounded m_yc;

  NO_COPY_CLASS(YCurveAnnuallyCompoundedTestCase);
};



class YCurveSwapTestCase : public CppUnit::TestCase
{
  /*
    Actually, bootstrapping tests are not included
    */
public:
  YCurveSwapTestCase() : m_yc(ito33::Date(2004, ito33::Date::Jan, 1)) { }

  void tearDown() { m_yc.Clear(); }

private:
  CPPUNIT_TEST_SUITE( YCurveSwapTestCase );
    CPPUNIT_TEST_EXCEPTION( AddCashDuplicate, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( AddSwapDuplicate, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( AddSwapTermBeforeCashTerm, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( Clear, ito33::Exception );

    CPPUNIT_TEST( Perturb );

  CPPUNIT_TEST_SUITE_END();

  void AddCashDuplicate();
  void AddSwapDuplicate();
  void AddSwapTermBeforeCashTerm();
  void Clear();

  void Perturb();

  ito33::finance::YieldCurveSwap m_yc;

  NO_COPY_CLASS(YCurveSwapTestCase);
};

#endif // #ifndef _ITO33_TEST_TESTYCURVE_H_
