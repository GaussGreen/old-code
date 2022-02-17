/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testmandatory.h
// Purpose:     header file for mandatory test
// Author:      Ito33 Canada
// Created:     May 11, 2005
// RCS-ID:      $Id: testmandatory.h,v 1.3 2006/06/10 23:15:13 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_MANDATORY_H_
#define _ITO33_TEST_MANDATORY_H_

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/exception.h"

class PepsAveragingPeriodTest : public CppUnit::TestCase { 

public: 
  PepsAveragingPeriodTest() {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( PepsAveragingPeriodTest );    

  CPPUNIT_TEST_EXCEPTION ( AveragingPeriodTooLarge, ito33::Exception );
  CPPUNIT_TEST_EXCEPTION ( NegativeNumberOfSampling, ito33::Exception );
  //CPPUNIT_TEST_EXCEPTION ( TooManySamplesUsed, ito33::Exception );

  CPPUNIT_TEST_EXCEPTION ( SetConversionRatioAveragingWhenStockAveragingPeriod,
    ito33::Exception );

  CPPUNIT_TEST_EXCEPTION ( SetStockAveragingWhenConversionAveragingPeriod, 
    ito33::Exception );
  
  CPPUNIT_TEST_EXCEPTION ( GetStockAveragingWhenConversionAveragingPeriod,
    ito33::Exception );
  
  CPPUNIT_TEST_EXCEPTION ( GetConversionRatioAveragingWhenStockAveragingPeriod,
    ito33::Exception );
   
  CPPUNIT_TEST_EXCEPTION ( SetNegativeStockAverage, ito33::Exception );
  CPPUNIT_TEST_EXCEPTION ( SetNegativeConversionRatioAverage, ito33::Exception );

  CPPUNIT_TEST( DumpStockAverage );
  CPPUNIT_TEST( DumpConversionRatioAverage );

  CPPUNIT_TEST_SUITE_END();

  void AveragingPeriodTooLarge();
  void NegativeNumberOfSampling();
  void TooManySamplesUsed();

  void SetConversionRatioAveragingWhenStockAveragingPeriod();
  void SetStockAveragingWhenConversionAveragingPeriod();
  void GetStockAveragingWhenConversionAveragingPeriod();
  void GetConversionRatioAveragingWhenStockAveragingPeriod();
  void SetNegativeStockAverage();
  void SetNegativeConversionRatioAverage();

  void DumpStockAverage();
  void DumpConversionRatioAverage();

  NO_COPY_CLASS( PepsAveragingPeriodTest  );
}; // Class PepsAveragingPeriodTest 

#endif
