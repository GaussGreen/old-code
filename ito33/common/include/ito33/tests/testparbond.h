/////////////////////////////////////////////////////////////////////////////
// Name:        testparbond.h
// Purpose:     unit test for par bond
// Author:      ZHANG Yunzhi
// Created:     05/23/2005
// RCS-ID:      $Id: testparbond.h,v 1.1 2005/06/08 15:41:49 zhang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class ParBondTest : public CppUnit::TestCase
{

public:
  ParBondTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( ParBondTest );
    CPPUNIT_TEST_EXCEPTION( BadYTM, ito33::Exception );   
    CPPUNIT_TEST_EXCEPTION( BadSpread, ito33::Exception );   
    CPPUNIT_TEST_EXCEPTION( BadFrequency, ito33::Exception );    
    CPPUNIT_TEST_EXCEPTION( BadRecoveryRate, ito33::Exception );   
    CPPUNIT_TEST_EXCEPTION( BadDayCountConvention, ito33::Exception );   
    CPPUNIT_TEST_EXCEPTION( BadMaturity, ito33::Exception );    
    CPPUNIT_TEST_EXCEPTION( BadFrequencyMaturity, ito33::Exception );

    CPPUNIT_TEST( GetCouponTrivial );
    CPPUNIT_TEST( GetCouponAndCashFlowStream );
  CPPUNIT_TEST_SUITE_END();

  void BadRecoveryRate();
  void BadYTM();
  void BadSpread();
  void BadFrequency();
  void BadDayCountConvention();
  void BadFrequencyMaturity();
  void BadMaturity();

  void GetCouponTrivial();
  void GetCouponAndCashFlowStream();


  NO_COPY_CLASS(ParBondTest);
};
