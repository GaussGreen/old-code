/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/testbondlike_putcallprice.cpp
// Purpose:     test file for put and call prices
// Author:      Nabil
// Created:     2005/02/27
// RCS-ID:      $Id: testbondlike_putcallprice.cpp,v 1.9 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"
#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callperiod.h"

#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/yieldcurve_flat.h"

#include "ito33/tests/testbondlike.h"

//#define TEST_YTPDEBUG

const double YTP_RELATIVE_ERROR = 1.e-4;

using namespace ito33;

using namespace ito33::finance;

shared_ptr<SessionData> put_call_price_InitSessionData()
{
   Date valuationDate = Date(2003, Date::Jan, 14);

   shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  // Setup the equity, and attach to session
  shared_ptr<Equity> pEquity(new Equity(95., pCurrency));
  pEquity->SetBorrowCurve(shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
  
  // Setup the issuer, and attach to the session  
  shared_ptr<RateData> pRateData(new RateData() );
  pRateData->SetYieldCurve( pCurrency, 
    shared_ptr<YieldCurve>(new YieldCurveFlat(0.04)) );
                
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

shared_ptr<BondTerms> put_call_price_InitBondTerms()
{
  Date
    issueDate = Date(2002, Date::Jan, 14),
    firstPaymentDate = Date(2003, Date::Jan, 14),
    maturityDate = Date(2053, Date::Jan, 14);

  double
    dIssuePrice = 1,
    dNominal = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  // Create the bondterms ============================================

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                       maturityDate, dNominal, dRedemptionRate,
                       dRecoveryRate)
      );
  
  // Add cash Flow Stream ============================================

  double dCouponRate = 0.02;

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate, firstPaymentDate, maturityDate, dCouponRate,
                      Date::DayCountConvention_30360, Frequency_SemiAnnual 
                    )
              );
    
  bc->SetCashDistribution(pInterests);

  return bc;
}

shared_ptr<BondTerms> put_call_price_InitAccretingBondTerms()
{
  Date
    issueDate = Date(2002, Date::Jun, 10),
    firstPaymentDate = Date(2003, Date::Jun, 10),
    maturityDate = Date(2012, Date::Jun, 10);

  double
    dIssuePrice = 1,
    dNominal = 5000,
    dRedemptionRate = 1.16805,
    dRecoveryRate = 0.;

  // Create the bondterms ============================================

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                       maturityDate, dNominal, dRedemptionRate,
                       dRecoveryRate)
      );

  bc->SetAccretingBond(0.025);
  bc->SetYieldCompoundingFrequency(Frequency_Annual);

  // Add cash Flow Stream ============================================
  
  double dCouponRate = 0.01;

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate, firstPaymentDate, maturityDate, dCouponRate,
                      Date::DayCountConvention_Act365, Frequency_Annual 
                    )
              );
    
  bc->SetCashDistribution(pInterests);

  return bc;
}

void BondLikeTest::ComputePutPrice_Exception_NoPut()
{
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      ( shared_ptr<ConversionPeriod>
        ( new ConversionPeriod(Date(2001, Date::Jan, 1), Date(2003, Date::Feb, 2), 1) )
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));
  pCB->SetSessionData(pSessionData);

  pCB->ComputePutPrice(Date(2002, Date::Feb, 1));
}


void BondLikeTest::ComputePutPrice_Exception_NotPutDate()
{
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      ( shared_ptr<ConversionPeriod>
        ( new ConversionPeriod(Date(2001, Date::Jan, 1), Date(2003, Date::Feb, 2), 1) )
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));
  pCB->SetSessionData(pSessionData);

  shared_ptr<PutSchedule> pPuts(new PutSchedule());
  pPuts->AddPutWithStrike(Date(2004, Date::Jan, 2), 1);
  pCB->ComputePutPrice(Date(2002, Date::Feb, 1));
}


void BondLikeTest::ComputePutPrice_Value()
{
  shared_ptr<PutSchedule> pPuts;
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      (
        shared_ptr<ConversionPeriod>
        ( new ConversionPeriod(Date(2001, Date::Jan, 1), Date(2003, Date::Feb, 2), 1) )
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));


  pPuts = make_ptr( new PutSchedule() );
  pPuts->AddPutWithStrike(Date(2003, Date::Jan, 2), 1);
  pPuts->AddPutWithStrike(Date(2003, Date::Jan, 14), 1);
  pPuts->AddPutWithStrike(Date(2004, Date::Jan, 15), 1);
  pPuts->AddPutWithStrike(Date(2006, Date::Jan, 14), 1);

  pCB->SetSessionData(pSessionData);
  pCB->SetPutSchedule(pPuts);

  double dPut;

  dPut = pCB->ComputePutPrice(Date(2003, Date::Jan, 2));  
  CPPUNIT_ASSERT( fabs(dPut - 112.1266) < 1.e-3 );


  dPut = pCB->ComputePutPrice(Date(2003, Date::Jan, 14));
  CPPUNIT_ASSERT( dPut == 112.2 );

  dPut = pCB->ComputePutPrice(Date(2004, Date::Jan, 15));
  CPPUNIT_ASSERT( fabs(dPut - 110.006) < 1.e-3 );

  dPut = pCB->ComputePutPrice(Date(2006, Date::Jan, 14));
  CPPUNIT_ASSERT( dPut == 111.1 );

  /// set yield to Put
  pPuts = make_ptr( new PutSchedule() );
  pPuts->AddPutWithYield(Date(2004, Date::Jan, 14), 0.03);
  pPuts->AddPutWithYield(Date(2004, Date::Jan, 15), 0.03);

  pCB->SetPutSchedule(pPuts);
  
  dPut = pCB->ComputePutPrice(Date(2004, Date::Jan, 14));

  CPPUNIT_ASSERT( fabs(dPut - 113.367) < 1.e-3);

  dPut = pCB->ComputePutPrice(Date(2004, Date::Jan, 15));
  CPPUNIT_ASSERT( fabs(dPut - 112.276) < 1.e-3);  
}

void BondLikeTest::ComputePutPrice_Value_AccretingBond()
{
  shared_ptr<PutSchedule> pPuts;
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitAccretingBondTerms());

  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      (
        shared_ptr<ConversionPeriod>
        ( new ConversionPeriod(Date(2001, Date::Jan, 1), Date(2003, Date::Feb, 2), 1) )
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));

  pPuts = make_ptr( new PutSchedule() );
  
  pPuts->AddPutWithStrike(Date(2004, Date::Jan, 15), 1);
  pPuts->AddPutWithStrike(Date(2006, Date::Jun, 10), 1);

  pCB->SetSessionData(pSessionData);
  pCB->SetPutSchedule(pPuts);

  double dPut;

  dPut = pCB->ComputePutPrice(Date(2004, Date::Jan, 15));
  CPPUNIT_ASSERT( fabs(dPut - 5149.7380) < 1.e-3 );

  dPut = pCB->ComputePutPrice(Date(2006, Date::Jun, 10));

  CPPUNIT_ASSERT( fabs(dPut - 5360.73595) < 1.e-3);

  /// set yield to Put
  pPuts = make_ptr( new PutSchedule() );
  pPuts->AddPutWithYield(Date(2004, Date::Jan, 15), 0.025);
  pPuts->AddPutWithStrike(Date(2006, Date::Jun, 10), 1);

  pCB->SetPutSchedule(pPuts);
  
  dPut = pCB->ComputePutPrice(Date(2004, Date::Jan, 15));

  CPPUNIT_ASSERT( fabs(dPut - 5150.7488) < 1.e-3);

  dPut = pCB->ComputePutPrice(Date(2006, Date::Jun, 10));

  CPPUNIT_ASSERT( fabs(dPut - 5360.73595) < 1.e-3);
}

void BondLikeTest::ComputeCallPrice_Exception_NoCall()
{
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<Bond> pBond(new Bond(pTerms));
  pBond->SetSessionData(pSessionData);

  pBond->ComputeCallPrice(Date(2003, Date::Feb, 1));
}


void BondLikeTest::ComputeCallPrice_Exception_NotCallDate()
{
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      (
        shared_ptr<ConversionPeriod>
        ( new ConversionPeriod(Date(2001, Date::Jan, 1),
                             Date(2003, Date::Feb, 2), 1) )
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));
  pCB->SetSessionData(pSessionData);

  shared_ptr<CallSchedule> pCalls(new CallSchedule());
  pCalls->AddCallPeriod(shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike(Date(2003, Date::Jan, 2),
                                             Date(2005, Date::Jan, 2),
                                             1))
                        );
  pCalls->AddCallPeriod(shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike(
                              Date(2007, Date::Jan, 2),
                              Date(2009, Date::Jan, 2),
                              1) )
                        );
  pCB->ComputeCallPrice(Date(2006, Date::Feb, 1));
}


void BondLikeTest::ComputeCallPrice_Value()
{
  shared_ptr<CallSchedule> pCalls;
  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<Bond> pBond(new Bond(pTerms));

  pCalls = make_ptr( new CallSchedule() );
  pCalls->AddCallPeriod(shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike(
                              Date(2003, Date::Jan, 2),
                              Date(2005, Date::Jan, 2),
                              1))
                        );
  pCalls->AddCallPeriod(shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike(
                              Date(2007, Date::Jan, 2),
                              Date(2009, Date::Jan, 2),
                              1.01))
                       );

  pBond->SetSessionData(pSessionData);
  pBond->SetCallSchedule(pCalls);

  double dCall;

  std::cout.precision(16);
  dCall = pBond->ComputeCallPrice(Date(2003, Date::Jan, 14));
  CPPUNIT_ASSERT( dCall == 112.2 );


  dCall = pBond->ComputeCallPrice(Date(2003, Date::Jan, 15));
  CPPUNIT_ASSERT( fabs(dCall - 110.00611111) < 1.e-6 );

  dCall = pBond->ComputeCallPrice(Date(2007, Date::Jan, 14));
  CPPUNIT_ASSERT( fabs(dCall - 112.2) < 1.e-9 );

  // activate forfeit coupon
  pCalls->SetForfeitCoupon(true);
  dCall = pBond->ComputeCallPrice(Date(2007, Date::Jan, 14));
  CPPUNIT_ASSERT( dCall == 111.1 );
  dCall = pBond->ComputeCallPrice(Date(2008, Date::Feb, 14));
  CPPUNIT_ASSERT( fabs(dCall - 111.283333) < 1.e-6 );

  // deactivate keep accrued
  pCalls->SetKeepAccrued(false);;
  dCall = pBond->ComputeCallPrice(Date(2008, Date::Feb, 14));
  CPPUNIT_ASSERT( dCall == 111.1 );


  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      (
        shared_ptr<ConversionPeriod>
        (new ConversionPeriod(Date(2001, Date::Jan, 1), Date(2003, Date::Feb, 2), 1))
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));

  pCalls = make_ptr( new CallSchedule() );
  pCalls->AddCallPeriod(shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithYield(
                              Date(2003, Date::Jan, 2),
                              Date(2005, Date::Jan, 2),
                              .0201 ))
                       );
  pCalls->AddCallPeriod(shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithYield(
                              Date(2007, Date::Jan, 2),
                              Date(2009, Date::Jan, 2),
                              .0201))
                       );
  pCB->SetSessionData(pSessionData);
  pCB->SetCallSchedule(pCalls);

  std::cout.precision(16);
  dCall = pCB->ComputeCallPrice(Date(2003, Date::Jan, 14));
  CPPUNIT_ASSERT( fabs(dCall - 112.222) < 1.e-3 );


  dCall = pCB->ComputeCallPrice(Date(2003, Date::Jan, 15));
  CPPUNIT_ASSERT( fabs(dCall - 110.028) < 1.e-3 );

  dCall = pCB->ComputeCallPrice(Date(2007, Date::Jan, 14));
  CPPUNIT_ASSERT( fabs(dCall - 111.169) < 1.e-2 );
  
}


void BondLikeTest::ComputeAccrued_Exception()
{
  shared_ptr<BondTerms> pTerms(put_call_price_InitBondTerms());

  shared_ptr<Bond> pBond(new Bond(pTerms));

  // Should throw exception as session is not set
  pBond->GetAccruedInterestValue();
}


void BondLikeTest::ComputeAccrued_Value()
{
  double dAccrued;

  shared_ptr<SessionData> pSessionData(put_call_price_InitSessionData());

  // Zero coupon bond  
  shared_ptr<BondTerms> 
    pTerms( new BondTerms(Date(2003, Date::Jan, 14),
                           1,
                           Date(2010, Date::Feb, 14),
                           110,
                           1,
                           0.1)
          );

  pTerms->SetYieldCompoundingFrequency(Frequency_Annual);
  pTerms->SetYieldDayCountConvention(Date::DayCountConvention_ActAct);

  shared_ptr<Bond> pBond(new Bond(pTerms));
  pBond->SetSessionData(pSessionData);
  pSessionData->SetValuationDate(Date(2003, Date::Jan, 15));

  // no coupon
  dAccrued = pBond->GetAccruedInterestValue();
  CPPUNIT_ASSERT( fabs(dAccrued - 0.) < 1.e-6 );

  pTerms = put_call_price_InitBondTerms();
  pBond = make_ptr( new Bond(pTerms) );
  pBond->SetSessionData(pSessionData);

  pSessionData->SetValuationDate(Date(2003, Date::Jan, 15));
  dAccrued = pBond->GetAccruedInterestValue();
  CPPUNIT_ASSERT( fabs(dAccrued - 0.00611111) < 1.e-6 );

  shared_ptr<ConversionSchedule> pConversions(new ConversionSchedule());
  pConversions->AddConversionPeriod
      (
        shared_ptr<ConversionPeriod>
        ( new ConversionPeriod(Date(2001, Date::Jan, 1), Date(2003, Date::Feb, 2), 1) )
      );
  shared_ptr<ConvertibleBond> pCB(new ConvertibleBond(pTerms, pConversions));

  pSessionData->SetValuationDate(Date(2003, Date::Jun, 15));
  dAccrued = pBond->GetAccruedInterestValue();
  CPPUNIT_ASSERT( fabs(dAccrued - 1.1 * 151 / 180) < 1.e-2 );
}
