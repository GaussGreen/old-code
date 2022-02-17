/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testyieldtomaturity.h
// Purpose:     header file for yield-to-maturity test
// Author:      Nabil
// Created:     2005/02/17
// RCS-ID:      $Id: testyieldtomaturity.h,v 1.4 2005/05/24 10:31:03 zhang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_YIELDTOMATURITY_H_
#define _ITO33_TEST_YIELDTOMATURITY_H_

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/exception.h"

class YieldToMaturityTest : public CppUnit::TestCase { 

public:

  YieldToMaturityTest() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( YieldToMaturityTest );

    CPPUNIT_TEST ( YTMAtAnnualCouponDates );
    CPPUNIT_TEST ( YTMAtSemiAnnualCouponDates );
    CPPUNIT_TEST ( YTMAtQuarterCouponDates );
    CPPUNIT_TEST ( YTMAtBiMonthCouponDates );
    CPPUNIT_TEST ( YTMAtMonthCouponDates );
    
    CPPUNIT_TEST ( YTMForAZeroCouponBond );
    CPPUNIT_TEST ( YTMForAConvertibleBond );

  CPPUNIT_TEST_SUITE_END();

  void YTMAtAnnualCouponDates();
  void YTMAtSemiAnnualCouponDates();
  void YTMAtQuarterCouponDates();
  void YTMAtBiMonthCouponDates();
  void YTMAtMonthCouponDates();

  void YTMForAZeroCouponBond();  
  void YTMForAConvertibleBond();
 
  NO_COPY_CLASS( YieldToMaturityTest );

}; // Class YieldToMaturityTest

#endif // _ITO33_TEST_YIELDTOMATURITY_H_
