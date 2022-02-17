/////////////////////////////////////////////////////////////////////////////
// Name:        testmoneymarket.h
// Purpose:     input validation for moneymarket
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testmoneymarket.h,v 1.1 2004/11/19 16:16:02 yann Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class MoneyMarketTest : public CppUnit::TestCase
{

public:
  MoneyMarketTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( MoneyMarketTest );
    CPPUNIT_TEST_EXCEPTION( NoYieldCurve, ito33::Exception );
    CPPUNIT_TEST( Dump );
  CPPUNIT_TEST_SUITE_END();

  void NoYieldCurve();
  void Dump();

  NO_COPY_CLASS(MoneyMarketTest);
};
