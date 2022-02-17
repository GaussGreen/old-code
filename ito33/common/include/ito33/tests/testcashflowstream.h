/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testcashflowstream.h
// Purpose:     header file for cash flow stream test
// Author:      Zhang 
// Created:     24.06.2004
// RCS-ID:      $Id: testcashflowstream.h,v 1.15 2006/06/09 15:55:12 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/exception.h"

class CashFlowStreamTest : public CppUnit::TestCase { 

public: 
  CashFlowStreamTest() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( CashFlowStreamTest );

    CPPUNIT_TEST_EXCEPTION( BadPaymentDates, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NegativeAmount, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( UndefinedFrequency, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( LastDateBeforeContractingDate, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( FirstDateBeforeContractingDate, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( LeapYearAdjustException, ito33::Exception );

    CPPUNIT_TEST( Uniform_Accrued );
    CPPUNIT_TEST( Uniform_Accrued_From_SWX );    
    CPPUNIT_TEST( Uniform_Irregular_Extreme_Period_From_ISMA );
    CPPUNIT_TEST( Uniform_LastRegular );
    CPPUNIT_TEST( Uniform_LastIrregular );
    CPPUNIT_TEST( Uniform_EndOfMonth );
    CPPUNIT_TEST( Uniform_NoEndOfMonth );
    CPPUNIT_TEST( Uniform_InMonths );
    CPPUNIT_TEST( Uniform_ShortLastPayment );
    CPPUNIT_TEST( Uniform_OnePayment );
    CPPUNIT_TEST( Uniform_TwoPayments );
    CPPUNIT_TEST( Uniform_AlwaysRegularButNotEndOfMonth );
    CPPUNIT_TEST( General_Accrued );
    CPPUNIT_TEST( FirstAndLastCoupon );
    CPPUNIT_TEST( LeapYearAdjust );
    
    CPPUNIT_TEST( Uniform_Dump );

  CPPUNIT_TEST_SUITE_END();

  void BadPaymentDates();
  void NegativeAmount();
  void UndefinedFrequency();
  void LastDateBeforeContractingDate();
  void FirstDateBeforeContractingDate();

  void Uniform_Accrued();
  // Examples from the SWX documentation. (cf Mangue\ito33\docs\CashFlows)
  void Uniform_Accrued_From_SWX();
  // Examples from the ISMA documentation. (cf Mangue\ito33\docs\CashFlows)
  void Uniform_Irregular_Extreme_Period_From_ISMA();
  
  void Uniform_LastRegular();
  void Uniform_LastIrregular();
  void Uniform_EndOfMonth();
  void Uniform_NoEndOfMonth();
  void Uniform_InMonths();

  // Tests the ctor with the boolean bHasShortLastPayment
  void Uniform_ShortLastPayment();
  void Uniform_OnePayment();
  void Uniform_TwoPayments();

  void Uniform_AlwaysRegularButNotEndOfMonth();
  void General_Accrued();
  void FirstAndLastCoupon();
  void LeapYearAdjust();
  void LeapYearAdjustException();

  // Tests for the function Dump()
  void Uniform_Dump();
    
  NO_COPY_CLASS( CashFlowStreamTest );
}; // Class CashFlowStreamTest
