/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testyieldtoput.h
// Purpose:     header file for yield-to-put test
// Author:      Nabil
// Created:     2005/02/23
// RCS-ID:      $Id: testyieldtoput.h,v 1.5 2006/08/19 22:28:08 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_YIELDTOPUT_H_
#define _ITO33_TEST_YIELDTOPUT_H_

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/sessiondata.h"

class YieldToPutTest : public CppUnit::TestCase 
{

public:

  YieldToPutTest() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( YieldToPutTest );

    CPPUNIT_TEST ( YTPAtAnnualCouponDates );
    CPPUNIT_TEST ( YTPAtSemiAnnualCouponDates );
    CPPUNIT_TEST ( YTPAtQuarterCouponDates );
    CPPUNIT_TEST ( YTPAtBiMonthCouponDates );
    CPPUNIT_TEST ( YTPAtMonthCouponDates );
    
    CPPUNIT_TEST ( YTPForAZeroCouponBond );
    CPPUNIT_TEST ( YTPForAConvertibleBond );

    CPPUNIT_TEST_EXCEPTION ( YTPForBondWithoutPutSchedule, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( YTPForConvertibleBondWithoutPutSchedule, 
                             ito33::Exception);

  CPPUNIT_TEST_SUITE_END();

  void YTPAtAnnualCouponDates();
  void YTPAtSemiAnnualCouponDates();
  void YTPAtQuarterCouponDates();
  void YTPAtBiMonthCouponDates();
  void YTPAtMonthCouponDates();

  void YTPForAZeroCouponBond();  
  void YTPForAConvertibleBond();

  void YTPForBondWithoutPutSchedule();  
  void YTPForConvertibleBondWithoutPutSchedule();

  // common functions
  ito33::shared_ptr<ito33::finance::SessionData> InitSessionData();
  void YTPAtCouponDates(ito33::finance::Frequency frequency);
 
  NO_COPY_CLASS( YieldToPutTest );

}; // Class YieldToPutTest

#endif // _ITO33_TEST_YIELDTOPUT_H_
