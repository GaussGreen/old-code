/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testbondlike.h
// Purpose:     header file for bondlike test
// Author:      Zhang (converted to cppunit by David)
// Created:     24.06.2004
// RCS-ID:      $Id: testbondlike.h,v 1.15 2006/03/08 13:07:40 nabil Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_BONDLIKE_H_
#define _ITO33_TEST_BONDLIKE_H_

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/exception.h"

class BondLikeTest : public CppUnit::TestCase { 

public: 
  BondLikeTest() {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( BondLikeTest );

    // put
    CPPUNIT_TEST_EXCEPTION ( CheckPut_DuplicateDate, ito33::Exception );
    CPPUNIT_TEST ( CheckPut_Sort );
    CPPUNIT_TEST_EXCEPTION( CheckPutStrikePositive, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( CheckYieldToPutNotTooSmall, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION( CheckYieldToPutNotTooLarge, ito33::Exception);
    CPPUNIT_TEST ( PutDump );


    // call
    CPPUNIT_TEST ( CheckCall_Sort );
    CPPUNIT_TEST_EXCEPTION ( CheckCall_WrongSchedule, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( CheckCall_HardSoft, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( CheckCall_HardSoft1, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( NoticePeriodTooLong, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( NoticePeriodTooSmall, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( SetPremiumMakeWholeNegative, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( SetTriggerCheckPeriod, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckCallPeriodTriggerPositive, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckCallPeriodStartDateBeforeEndDate
                             , ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( CheckCallPeriodStrikePositive
                             , ito33::Exception );
    
    CPPUNIT_TEST_EXCEPTION( CheckYieldToCallNotTooSmall, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION( CheckYieldToCallNotTooLarge, ito33::Exception);
    
    
    CPPUNIT_TEST( CallPeriodDump );
    CPPUNIT_TEST( CallScheduleDump );

    CPPUNIT_TEST_EXCEPTION ( CheckCallFixedShareStartDateBeforeEndDate,
                             ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckCallFixedShareRatioPositive,
                             ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckCallFixedSharedTriggerRatePositive,
                             ito33::Exception);

    CPPUNIT_TEST( CallFixedShareDump );

    // conversion
     //CPPUNIT_TEST ( CheckConversion_Sort );
    CPPUNIT_TEST_EXCEPTION ( CheckConversion_WrongSchedule, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION (CheckConversionPeriodStartDateBeforeEndDate
                            , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckConversionPeriodRatioPositive
                             , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckConversionPeriodTriggerRatePositive
                           , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckConversionPeriodExtremeTriggerRatePositive
                           , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckConversionPeriodChangeRatePositive
                           , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckConversionPeriodChangeRateNegative
                             , ito33::Exception );
    CPPUNIT_TEST( ConversionPeriodDump );
    CPPUNIT_TEST( ConversionScheduleDump );


    //bondliketerms
    CPPUNIT_TEST_EXCEPTION ( CheckBondLikeTermsIssueDateBeforeMaturityDate
                          , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckbondLikeTermsIssuePricePositive
                          , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckBondLikeTermsRecoveryTooLarge
                          , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckbondLikeTermsRecoveryRateNegative
                        , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( CheckBondLikeTermsNominalPositive
                          , ito33::Exception);

    // BondTerms
    CPPUNIT_TEST_EXCEPTION ( BondTermsFloatingRatesThenCashFlowStream
                          , ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( BondTermsCashFlowStreamThenFloatingRates
                          , ito33::Exception);

    // Compute Put price
    CPPUNIT_TEST_EXCEPTION( ComputePutPrice_Exception_NoPut,
                            ito33::Exception);
    CPPUNIT_TEST_EXCEPTION( ComputePutPrice_Exception_NotPutDate,
                            ito33::Exception);
    CPPUNIT_TEST( ComputePutPrice_Value);
    CPPUNIT_TEST( ComputePutPrice_Value_AccretingBond);

    // Compute call price
    CPPUNIT_TEST_EXCEPTION( ComputeCallPrice_Exception_NoCall,
                            ito33::Exception);
    CPPUNIT_TEST_EXCEPTION( ComputeCallPrice_Exception_NotCallDate,
                            ito33::Exception);
    CPPUNIT_TEST( ComputeCallPrice_Value);

  
    // Compute accrued
    CPPUNIT_TEST_EXCEPTION( ComputeAccrued_Exception,
                            ito33::Exception);
    CPPUNIT_TEST( ComputeAccrued_Value);

  
    // Cross-currency (CC)
    CPPUNIT_TEST_EXCEPTION( CC_InvalidFixedFXRate,
                            ito33::Exception);
  
    // Fixed quanto (FQ)
    CPPUNIT_TEST_EXCEPTION( FQ_InvalidFXVolatility_High,
                            ito33::Exception);
    CPPUNIT_TEST_EXCEPTION( FQ_InvalidFXVolatility_Low,
                            ito33::Exception);

    CPPUNIT_TEST_EXCEPTION( FQ_InvalidCorrelationBetweenShareAndFXRate_High,
                            ito33::Exception);
    CPPUNIT_TEST_EXCEPTION( FQ_InvalidCorrelationBetweenShareAndFXRate_Low,
                            ito33::Exception);

  CPPUNIT_TEST_SUITE_END();

  // put
  void CheckPut_DuplicateDate();
  void CheckPut_Sort();
  void CheckPutStrikePositive();
  void CheckYieldToPutNotTooSmall();
  void CheckYieldToPutNotTooLarge();

  void PutDump();


  // call
  void CheckCall_Sort();
  void CheckCall_WrongSchedule();
  void CheckCall_HardSoft();
  void CheckCall_HardSoft1();
  void NoticePeriodTooLong();
  void NoticePeriodTooSmall();
  void SetPremiumMakeWholeNegative();
  void SetTriggerCheckPeriod();
  void CheckCallPeriodTriggerPositive();
  void CheckCallPeriodStartDateBeforeEndDate();
  void CheckCallPeriodStrikePositive();
  void CheckYieldToCallNotTooSmall();
  void CheckYieldToCallNotTooLarge();
  void CallPeriodDump();
  void CallScheduleDump();

  void CheckCallFixedShareStartDateBeforeEndDate();
  void CheckCallFixedShareRatioPositive();
  void CheckCallFixedSharedTriggerRatePositive();
  void CallFixedShareDump();

  // conversion
  // void CheckConversion_Sort();
  void CheckConversion_WrongSchedule();
  void CheckConversionPeriodStartDateBeforeEndDate();
  void CheckConversionPeriodRatioPositive();
  void CheckConversionPeriodTriggerRatePositive();
  void CheckConversionPeriodExtremeTriggerRatePositive();
  void CheckConversionPeriodChangeRatePositive();
  void CheckConversionPeriodChangeRateNegative();
  void ConversionPeriodDump();
  void ConversionScheduleDump();

 
  //bondliketerms
  void CheckBondLikeTermsIssueDateBeforeMaturityDate();
  void CheckbondLikeTermsIssuePricePositive();
  void CheckBondLikeTermsRecoveryTooLarge();
  void CheckbondLikeTermsRecoveryRateNegative();
  void CheckBondLikeTermsNominalPositive();

  //BondTerms
  void BondTermsFloatingRatesThenCashFlowStream();
  void BondTermsCashFlowStreamThenFloatingRates();

  // compute put price
  void ComputePutPrice_Exception_NoPut();
  void ComputePutPrice_Exception_NotPutDate();
  void ComputePutPrice_Value();
  void ComputePutPrice_Value_AccretingBond();

  // compute call price
  void ComputeCallPrice_Exception_NoCall();
  void ComputeCallPrice_Exception_NotCallDate();
  void ComputeCallPrice_Value();
  
  // compute accrued
  void ComputeAccrued_Exception();
  void ComputeAccrued_Value();

  // Cross-currency (CC)
  void CC_InvalidFixedFXRate();

  // Fixed Quanto (FQ)
  void FQ_InvalidFXVolatility_High();
  void FQ_InvalidFXVolatility_Low();
  void FQ_InvalidCorrelationBetweenShareAndFXRate_High();
  void FQ_InvalidCorrelationBetweenShareAndFXRate_Low();

  NO_COPY_CLASS( BondLikeTest );
}; // Class BondLikeTest

#endif
