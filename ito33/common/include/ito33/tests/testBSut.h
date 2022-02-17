/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testBSut.
// Purpose:     test file for black-scholes utility function
// Author:      Laurence
// Created:     22/09/2003
// RCS-ID:      $Id: testBSut.h,v 1.3 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_BSUT_H_
#define _ITO33_TEST_BSUT_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"

class BSutileTest : public CppUnit::TestCase { 
public: 
  BSutileTest() {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( BSutileTest );
    CPPUNIT_TEST ( CallPriceFS );
    CPPUNIT_TEST ( CallPriceFE );
    CPPUNIT_TEST ( CallPriceFR );
    CPPUNIT_TEST ( CallPriceSpot0 );
    CPPUNIT_TEST ( CallPriceRate0 );
    CPPUNIT_TEST ( CallPriceInfty );
    CPPUNIT_TEST ( PutPriceFS );
    CPPUNIT_TEST ( PutPriceFE );
    CPPUNIT_TEST ( PutPriceFR );
    CPPUNIT_TEST ( PutPriceSpot0 );
    CPPUNIT_TEST ( PutPriceRate0 );
    CPPUNIT_TEST ( PutPriceInfty );
    CPPUNIT_TEST ( VerifVega );
    CPPUNIT_TEST ( VegaCall );
    CPPUNIT_TEST ( ImpliedVol);
    CPPUNIT_TEST ( Delta );
    CPPUNIT_TEST ( DeltaTime );
//    CPPUNIT_TEST ( VerifDeltaDeriv );
    CPPUNIT_TEST ( GetStrike );
    CPPUNIT_TEST ( GetStrikeTime );
  CPPUNIT_TEST_SUITE_END();
  
//tests for function CalculateCallPrice
  void CallPriceFS();
  void CallPriceFE();
  void CallPriceFR();
  void CallPriceSpot0();
  void CallPriceRate0();
  void CallPriceInfty();
//tests for function CalculatePutPrice
  void PutPriceFS();
  void PutPriceFE();
  void PutPriceFR();
  void PutPriceSpot0();
  void PutPriceRate0();
  void PutPriceInfty();
//tests for function CalculateVegaCall
  void VerifVega();
  void VegaCall();
//tests for function CalculateImpliedVolatility
  void ImpliedVol();
//tests for function CalculateDelta
  void Delta();
  void DeltaTime();
//test for function DeltaDeriv
//  void VerifDeltaDeriv();
//tests for function GetStrikeFromDelta
  void GetStrike();
  void GetStrikeTime(); 

  NO_COPY_CLASS( BSutileTest );
}; // Class BSutileTest

#endif
