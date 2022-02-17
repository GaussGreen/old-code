/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/testbondlike.cpp
// Purpose:     test file for calls, puts and conversions
// Author:      Zhang (converted to cppunit by David)
// Created:     24.06.04
// RCS-ID:      $Id: testbondlike.cpp,v 1.30 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <iostream>
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

#include "ito33/finance/equity.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/call.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/callfixedshare.h"
#include "ito33/finance/bondlike/cocotype.h"
#include "ito33/finance/bondlike/bondterms.h"

#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/floatingrates.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/numeraire.h"

#include "ito33/tests/testbondlike.h"

#include "ito33/tests/utilexml.h"

using namespace ito33;

using namespace ito33::finance;

shared_ptr<SessionData> Init_CC_SessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01"); 

  shared_ptr<Numeraire> pCurrencyEUR( new Numeraire("EUR") );
  shared_ptr<Numeraire> pCurrencyUSD( new Numeraire("USD") );

  // Setup the equity, and attach to session
  shared_ptr<Equity> pEquity(new Equity(100., pCurrencyEUR));
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0) ) );
  
  // Setup the issuer, and attach to the session
  shared_ptr<RateData> pRateData( new RateData() );
  pRateData->SetYieldCurve( pCurrencyEUR, 
    shared_ptr<YieldCurve>(new YieldCurveFlat(0.04)) );
  pRateData->SetYieldCurve( pCurrencyUSD,
    shared_ptr<YieldCurve>(new YieldCurveFlat(0.02)) );
  
  // Spot FX rate
  double dSpotFXRate = 0.8;
  shared_ptr<SpotFXRates> pSpotFXRates( new SpotFXRates() );
  pSpotFXRates->SetFXRate(pCurrencyEUR, pCurrencyUSD, dSpotFXRate);
  pRateData->SetSpotFXRates(pSpotFXRates);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<finance::ConvertibleBond> 
InitCrossCurrencyCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  double
    dIssuePrice = 1,
    dParValue = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                       maturityDate, dParValue, dRedemptionRate,
                       dRecoveryRate) );
  
  

  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(false);
  
  conv->AddConversionPeriod(
    shared_ptr<ConversionPeriod>(new ConversionPeriod(issueDate, maturityDate, 1)));

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

  pConvertibleBond->SetNumeraire(shared_ptr<Numeraire>( new Numeraire("USD") ) );

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
 
  Date startCallDate(issueDate);
  Date endCallDate(maturityDate);  
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike
                             (startCallDate, endCallDate, 1) );
  pCallSchedule->AddCallPeriod( pCallPeriod );    
    
  pConvertibleBond->SetCallSchedule(pCallSchedule);

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

void BondLikeTest::CheckPut_DuplicateDate()
{
  PutSchedule puts;

  puts.AddPutWithStrike(100, 1);
  puts.AddPutWithStrike(100, 1);
}

void BondLikeTest::CheckPut_Sort()
{
  PutSchedule puts;

  puts.AddPutWithStrike(100, 1);
  puts.AddPutWithStrike(50, 1);
  puts.AddPutWithStrike(200, 1);

  const PutSchedule::Elements& p = puts.GetAll();
  PutSchedule::Elements::const_iterator iter = p.begin();

  CPPUNIT_ASSERT( iter->first == 50 );
  ++iter;
  CPPUNIT_ASSERT( iter->first == 100 );
  ++iter;
  CPPUNIT_ASSERT( iter->first == 200 );
  ++iter;
  CPPUNIT_ASSERT( iter == p.end() );
}

void BondLikeTest::CheckPutStrikePositive()
{
  PutSchedule puts;
  puts.AddPutWithStrike(100, -1);
}

void BondLikeTest::CheckYieldToPutNotTooSmall()
{
  PutSchedule puts;
  puts.AddPutWithYield(100,-.11);
}

void BondLikeTest::CheckYieldToPutNotTooLarge()
{
  PutSchedule puts;
  puts.AddPutWithYield(100,2.1);
}

void BondLikeTest::PutDump()
{

   std::ostringstream oss;

  Date   date1(2003,Date::Oct,1);
  double dStrike = 10;

  Date date2(2003,Date::Nov,1);
  double dYieldToPut = .1;

  Date date3(2003,Date::Dec,1);


  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<put_schedule>\n"
                "<keep_accrued>1</keep_accrued>\n"
                "<forfeit_coupon>0</forfeit_coupon>\n"
                "<puts>\n"
                "<put>\n"
                "<date>2003-10-01</date>\n"
                "<strike>10</strike>\n"
                "</put>\n"
                "<put>\n"
                "<date>2003-11-01</date>\n"
                "<yield>0.1</yield>\n"
                "</put>\n"
                "<put>\n"
                "<date>2003-12-01</date>\n"
                "<strike>10</strike>\n"
                "</put>\n"
                "</puts>\n"
                "</put_schedule>\n"
                "</root>"
              );


  ito33::XML::RootTag root("root",oss);

   PutSchedule puts;
   puts.AddPutWithStrike(date1, dStrike);
   puts.AddPutWithYield(date2, dYieldToPut);
   puts.AddPutWithStrike(date3, dStrike);

   puts.Dump(root);

}

void BondLikeTest::CheckCall_WrongSchedule()
{
  CallSchedule calls;

  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(100), Date(200), 100) ) );
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(1), Date(100), 100) ) );
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (50, 100, 100) ) );
}


void BondLikeTest::CheckCall_HardSoft()
{
  CallSchedule calls;

  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(100, 200, 100) );

  pCallPeriod->SetTrigger(1.);

  calls.AddCallPeriod(pCallPeriod);
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(1), Date(100), 100) ) );
}


void BondLikeTest::CheckCall_HardSoft1()
{
  CallSchedule calls;

  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   ( Date(1), Date(100), 100) ) );
  
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(100, 200, 100) );
  pCallPeriod->SetTrigger(1.);

  calls.AddCallPeriod(pCallPeriod);
}

void BondLikeTest::CheckCall_Sort()
{
  CallSchedule calls;

  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(100), Date(200), 100) ) );
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(5), Date(100), 100) ) );
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(100), Date(100), 100) ) );
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(300), Date(400), 100) ) );
  calls.AddCallPeriod( shared_ptr<CallPeriod>
                        ( CallPeriod::CreateWithStrike
                                   (Date(250), Date(300), 100) ) );

  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(1 , 5, 100) );

  pCallPeriod->SetTrigger(1.2);

  calls.AddCallPeriod(pCallPeriod);

  const CallSchedule::Elements& t(calls.GetAll());

  CallSchedule::Elements::const_iterator iter = t.begin();

  CPPUNIT_ASSERT( (*iter)->GetStartDate() == Date(1) );
  CPPUNIT_ASSERT( (*iter)->GetEndDate() == Date(5) );
  ++iter;

  CPPUNIT_ASSERT( (*iter)->GetStartDate() == Date(5) );
  CPPUNIT_ASSERT( (*iter)->GetEndDate() == Date(100) );
  ++iter;

  CPPUNIT_ASSERT( (*iter)->GetStartDate() == Date(100) );
  CPPUNIT_ASSERT( (*iter)->GetEndDate() == Date(100) );
  ++iter;

  CPPUNIT_ASSERT( (*iter)->GetStartDate() == Date(100) );
  CPPUNIT_ASSERT( (*iter)->GetEndDate() == Date(200) );
  ++iter;

  CPPUNIT_ASSERT( (*iter)->GetStartDate() == Date(250) );
  CPPUNIT_ASSERT( (*iter)->GetEndDate() == Date(300) );
  ++iter;

  CPPUNIT_ASSERT( (*iter)->GetStartDate() == Date(300) );
  CPPUNIT_ASSERT( (*iter)->GetEndDate() == Date(400) );
  ++iter;

  CPPUNIT_ASSERT( iter == t.end() );
}

//

void BondLikeTest::CheckConversion_WrongSchedule()
{
  ConversionSchedule conversions;

  conversions.AddConversionPeriod
              ( shared_ptr<ConversionPeriod>(new ConversionPeriod(Date(100), Date(200), 100) ) );

  conversions.AddConversionPeriod
              ( shared_ptr<ConversionPeriod>(new ConversionPeriod(Date(1), Date(100), 100) ) );

  conversions.AddConversionPeriod
              ( shared_ptr<ConversionPeriod>(new ConversionPeriod(50, 100, 100) ) );
}



void BondLikeTest::NoticePeriodTooLong()
{
  CallSchedule call;
  call.SetNoticePeriod(400);

}

void BondLikeTest::NoticePeriodTooSmall()
{
  CallSchedule call;
  call.SetNoticePeriod(0);
}


void BondLikeTest::SetPremiumMakeWholeNegative()
{
  CallSchedule call;
  call.SetPremiumMakeWhole(-.2);
}

void BondLikeTest::SetTriggerCheckPeriod()
{
  CallSchedule call;
  call.SetTriggerCheckPeriod(40,20);
}


void  BondLikeTest::CheckCallPeriodTriggerPositive()
{
  Date startDate(2003,Date::Oct,1);
  Date endDate(2004,Date::Oct,1);
  double dStrike = 1.0;

  shared_ptr<CallPeriod>
    pCallPeriod( CallPeriod::CreateWithStrike
                            (startDate,endDate,dStrike) );
  
  pCallPeriod->SetTrigger(-.10);
}
void BondLikeTest::CheckCallPeriodStrikePositive()
{
  Date startDate(2003,Date::Oct,1);
  Date endDate(2004,Date::Oct,1);
  double dStrike = -1.0;

  shared_ptr<CallPeriod>
    pCallPeriod( CallPeriod::CreateWithStrike
                            (startDate,endDate,dStrike) );
}

void  BondLikeTest::CheckCallPeriodStartDateBeforeEndDate()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2003,Date::Oct,1);
  double dStrike = 0;

  shared_ptr<CallPeriod>
    pCallPeriod( CallPeriod::CreateWithStrike
                            (startDate,endDate,dStrike) );
}

void BondLikeTest::CheckYieldToCallNotTooSmall()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);

  shared_ptr<CallPeriod>
    pCallPeriod( CallPeriod::CreateWithYield
                            (startDate, endDate, -.11) );
}

void BondLikeTest::CheckYieldToCallNotTooLarge()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2003,Date::Oct,1);

  shared_ptr<CallPeriod>
    pCallPeriod( CallPeriod::CreateWithYield
                             (startDate, endDate, 3.1) );
}

void BondLikeTest::CallPeriodDump()
{
  
  std::ostringstream oss;

  Date startDate(2003,Date::Oct,1);
  Date endDate(2004,Date::Oct,1);
  double dStrike = 10;
  double dTriggerRate = .3;


  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<call_period>\n"
                "<start_date>2003-10-01</start_date>\n"
                "<end_date>2004-10-01</end_date>\n"
                "<strike>10</strike>\n"
                "<trigger_rate>0.3</trigger_rate>\n"
                "</call_period>\n"
                "</root>"
              );


  ito33::XML::RootTag root("root",oss);

  shared_ptr<CallPeriod>
    pCallPeriod( CallPeriod::CreateWithStrike
                             (startDate,endDate,dStrike) );

  pCallPeriod->SetTrigger(dTriggerRate);

  pCallPeriod->Dump(root);
}


void BondLikeTest::CallScheduleDump()
{
  std::ostringstream oss;

  Date startDate(2003,Date::Oct,1);
  Date endDate(2004,Date::Oct,1);
  double dStrike = 10;
  double dTriggerRate = .3;

  shared_ptr<CallPeriod>  
    pCallPeriod1( CallPeriod::CreateWithStrike
                               (startDate, endDate, dStrike) );
  pCallPeriod1->SetTrigger(dTriggerRate);

  startDate.AddYears(1);
  endDate.AddYears(1);
  
  shared_ptr<CallPeriod>  
    pCallPeriod2( CallPeriod::CreateWithStrike
                               (startDate, endDate, dStrike) );

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<call_schedule>\n"
                "<keep_accrued>1</keep_accrued>\n"
                "<forfeit_coupon>0</forfeit_coupon>\n"
                "<trigger_as_percentage_of>principal</trigger_as_percentage_of>\n"
                "<call_periods>\n"
                "<call_period>\n"
                "<start_date>2003-10-01</start_date>\n"
                "<end_date>2004-10-01</end_date>\n"
                "<strike>10</strike>\n"
                "<trigger_rate>0.3</trigger_rate>\n"
                "</call_period>\n"
                "<call_period>\n"
                "<start_date>2004-10-01</start_date>\n"
                "<end_date>2005-10-01</end_date>\n"
                "<strike>10</strike>\n"
                "</call_period>\n"
                "</call_periods>\n"
                "</call_schedule>\n"
                "</root>"
              );


  ito33::XML::RootTag root("root",oss);

  CallSchedule callSchedule;
  callSchedule.AddCallPeriod(pCallPeriod1);
  callSchedule.AddCallPeriod(pCallPeriod2);

  callSchedule.Dump(root);

}


void BondLikeTest::CheckCallFixedShareStartDateBeforeEndDate()
{
  Date startDate = Date(2003,Date::Oct,1);
  Date endDate   = Date(2004,Date::Oct,1);
  double dRatio = 1.0;

  CallFixedShare cfs(endDate,startDate,dRatio);

}

void BondLikeTest::CheckCallFixedShareRatioPositive()
{
 Date startDate = Date(2003,Date::Oct,1);
 Date endDate   = Date(2004,Date::Oct,1);
 double dRatio  = -1.0;

 CallFixedShare cfs(startDate,endDate,dRatio);


}

void BondLikeTest::CheckCallFixedSharedTriggerRatePositive()
{
  Date startDate = Date(2003,Date::Oct,1);
  Date endDate   = Date(2004,Date::Oct,1);
  double dRatio       = 1.0;
  double dTriggerRate = -1.0;

  CallFixedShare cfs(startDate,endDate,dRatio);
  cfs.SetTrigger( dTriggerRate );
}

void BondLikeTest::CallFixedShareDump()
{

 std::ostringstream oss;
 
 Date startDate = Date(2003,Date::Oct,1);
 Date endDate   = Date(2004,Date::Oct,1);
 double dRatio       = .2;
 double dTriggerRate = .3;


  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<call_fixed_share>\n"
                "<keep_accrued>1</keep_accrued>\n"
                "<forfeit_coupon>0</forfeit_coupon>\n"
                "<trigger_as_percentage_of>principal</trigger_as_percentage_of>\n"
                "<start_date>2003-10-01</start_date>\n"
                "<end_date>2004-10-01</end_date>\n"
                "<ratio>0.2</ratio>\n"
                "<trigger_rate>0.3</trigger_rate>\n"
                "</call_fixed_share>\n"
                "</root>"
              );


  ito33::XML::RootTag root("root",oss);

  CallFixedShare cfs(startDate,endDate,dRatio);
  cfs.SetTrigger( dTriggerRate );

  cfs.Dump(root);

}


 void BondLikeTest::CheckConversionPeriodStartDateBeforeEndDate()
 {
   Date startDate(2004,Date::Oct,1);
   Date endDate(2005,Date::Oct,1);
   double dRatio = 1.0;

   ConversionPeriod convPeriod(endDate,startDate,dRatio);
 }

void BondLikeTest::CheckConversionPeriodRatioPositive()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = -1.0;

  ConversionPeriod convPeriod(startDate,endDate,dRatio);
}

void BondLikeTest::CheckConversionPeriodTriggerRatePositive()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = 1.0;
  double dTriggerRate = -1.0;
  double dChangeRate = 1.0;
  double dExtremeTriggerRate = 1.0;

  ConversionPeriod convPeriod(startDate,endDate,dRatio);

  
  convPeriod.SetCoCo(dTriggerRate,
    CoCoType_CheckAnyTimeAndConvertOnCheckDate,
    dChangeRate,dExtremeTriggerRate);

}


void BondLikeTest::CheckConversionPeriodExtremeTriggerRatePositive()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = 1.0;
  double dTriggerRate = 1.0;
  double dChangeRate  = 1.0;
  double dExtremeTriggerRate = -1.0;

  ConversionPeriod convPeriod(startDate,endDate,dRatio);


  convPeriod.SetCoCo(dTriggerRate,
    CoCoType_CheckAnyTimeAndConvertOnCheckDate,
    dChangeRate,dExtremeTriggerRate);

}


void BondLikeTest::CheckConversionPeriodChangeRatePositive()
{
  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = 1.0;
  double dTriggerRate = 1.0;
  double dChangeRate  = 1.0;
  double dExtremeTriggerRate = -1.0;

  ConversionPeriod convPeriod(startDate,endDate,dRatio);


  convPeriod.SetCoCo(dTriggerRate,
    CoCoType_CheckAnyTimeAndConvertOnCheckDate,
    dChangeRate,dExtremeTriggerRate);

}
 
void BondLikeTest::CheckConversionPeriodChangeRateNegative()
{
 Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = 1.0;
  double dTriggerRate = -1.0;
  double dChangeRate  = -1.0;
  double dExtremeTriggerRate = 1.0;

  ConversionPeriod convPeriod(startDate,endDate,dRatio);


  convPeriod.SetCoCo(dTriggerRate,
    CoCoType_CheckAnyTimeAndConvertOnCheckDate,
    dChangeRate,dExtremeTriggerRate);

 }
 

void BondLikeTest::ConversionPeriodDump()
{
  std::ostringstream oss;
 
  ExpectedXML expected(oss,
            "<?xml version=\"1.0\"?>"
            "<root>\n"
            "<conversion_period>\n"
            "<start_date>2004-10-01</start_date>\n"
            "<end_date>2005-10-01</end_date>\n"
            "<ratio>1</ratio>\n"
            "<cash_value>0</cash_value>\n"
            "<trigger_rate>1</trigger_rate>\n"
            "<change_rate_in_trigger>1</change_rate_in_trigger>\n"
            "<extreme_trigger_rate>1</extreme_trigger_rate>\n"
            "<coco_type>check_any_time_and_convert_on_check_date</coco_type>\n"
            "</conversion_period>\n"
            "</root>"
         );

  ito33::XML::RootTag root("root",oss);

  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = 1.0;
  double dTriggerRate = 1.0;
  double dChangeRate  = 1.0;
  double dExtremeTriggerRate = 1.0;

  ConversionPeriod convPeriod(startDate,endDate,dRatio);


  convPeriod.SetCoCo(dTriggerRate,
                     CoCoType_CheckAnyTimeAndConvertOnCheckDate,
                     dChangeRate,dExtremeTriggerRate);

  convPeriod.Dump(root);

}


void BondLikeTest::ConversionScheduleDump()
{
  std::ostringstream oss;
 
  ExpectedXML expected(oss,
            "<?xml version=\"1.0\"?>"
            "<root>\n"
            "<conversion_schedule>\n"
            "<keep_accrued>0</keep_accrued>\n"
            "<forfeit_coupon>0</forfeit_coupon>\n"
            "<conversion_periods>\n"
            "<conversion_period>\n"
            "<start_date>2004-10-01</start_date>\n"
            "<end_date>2005-10-01</end_date>\n"
            "<ratio>1</ratio>\n"
            "<cash_value>0</cash_value>\n"
            "<trigger_rate>1</trigger_rate>\n"
            "<change_rate_in_trigger>1</change_rate_in_trigger>\n"
            "<extreme_trigger_rate>1</extreme_trigger_rate>\n"
            "<coco_type>check_any_time_and_convert_on_check_date</coco_type>\n"
            "</conversion_period>\n"
            "</conversion_periods>\n"
            "</conversion_schedule>\n"
            "</root>"
         );

  ito33::XML::RootTag root("root",oss);

  Date startDate(2004,Date::Oct,1);
  Date endDate(2005,Date::Oct,1);
  double dRatio = 1.0;
  double dTriggerRate = 1.0;
  double dChangeRate  = 1.0;
  double dExtremeTriggerRate = 1.0;

  ConversionSchedule convSchedule;

  shared_ptr<ConversionPeriod> 
    pConvPeriod( new ConversionPeriod(startDate,endDate,dRatio) );

  pConvPeriod->SetCoCo(dTriggerRate,
                       CoCoType_CheckAnyTimeAndConvertOnCheckDate,
                       dChangeRate,dExtremeTriggerRate);

  convSchedule.AddConversionPeriod(pConvPeriod);

  convSchedule.Dump(root);
}

void BondLikeTest::CheckBondLikeTermsIssueDateBeforeMaturityDate()
{
  Date issueDate(2004,Date::Oct,1);
  Date maturityDate(2005,Date::Oct,1);
  double dIssuePrice   = 1.0;
  double dNominal      = 100.0;
  double dRecoveryRate = 1.0;

  BondLikeTerms bondLikeTerms(maturityDate, 
                              dIssuePrice, 
                              issueDate, 
                              dNominal,
                              dRecoveryRate);

}

void BondLikeTest::CheckbondLikeTermsIssuePricePositive()
{

  Date issueDate(2004,Date::Oct,1);
  Date maturityDate(2005,Date::Oct,1);
  double dIssuePrice   = -1.0;
  double dNominal      = 100.0;
  double dRecoveryRate = 1.0;

  BondLikeTerms bondLikeTerms(maturityDate, 
                              dIssuePrice, 
                              issueDate, 
                              dNominal,
                              dRecoveryRate);


}
void BondLikeTest::CheckBondLikeTermsRecoveryTooLarge()
{

  Date issueDate(2004,Date::Oct,1);
  Date maturityDate(2005,Date::Oct,1);
  double dIssuePrice   = 1.0;
  double dNominal      = 100.0;
  double dRecoveryRate = 5.0;

  BondLikeTerms bondLikeTerms(maturityDate, 
                              dIssuePrice, 
                              issueDate, 
                              dNominal,
                              dRecoveryRate);

}

void BondLikeTest::CheckbondLikeTermsRecoveryRateNegative()
{

  Date issueDate(2004,Date::Oct,1);
  Date maturityDate(2005,Date::Oct,1);
  double dIssuePrice   = 1.0;
  double dNominal      = 100.0;
  double dRecoveryRate = -1.0;

  BondLikeTerms bondLikeTerms(maturityDate, 
                              dIssuePrice, 
                              issueDate, 
                              dNominal,
                              dRecoveryRate);

}

void BondLikeTest::CheckBondLikeTermsNominalPositive()
{
  Date issueDate(2004,Date::Oct,1);
  Date maturityDate(2005,Date::Oct,1);
  double dIssuePrice   = 1.0;
  double dNominal      = -100.0;
  double dRecoveryRate = 1.0;

  BondLikeTerms bondLikeTerms(maturityDate, 
                              dIssuePrice, 
                              issueDate, 
                              dNominal,
                              dRecoveryRate);

}

void BondLikeTest::BondTermsFloatingRatesThenCashFlowStream()
{
  Date
    issueDate = Date(2002, Date::Jan, 14),
    maturityDate = Date(2004, Date::Jan, 14);

  double
    dIssuePrice = 1,
    dParValue = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                       maturityDate, dParValue, dRedemptionRate,
                       dRecoveryRate)
      );

  // Floating rates
  
  Date lastPaymentDate = maturityDate;
  Date firstPaymentDate = maturityDate;
  Frequency paymentFrequency = Frequency_Quarterly;
  const int iInterval = 12 / paymentFrequency; 

  firstPaymentDate = AddMonthsAdjustedForEndOfMonth(firstPaymentDate, 
    - 3 * iInterval);
  
  shared_ptr<FloatingRates> FloatingInterests(
    new FloatingRates(0., issueDate, firstPaymentDate, lastPaymentDate,
          paymentFrequency) );

  // Fixed coupon added
  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  Date 
    fixedCouponDate = AddMonthsAdjustedForEndOfMonth(firstPaymentDate, 
      - iInterval);
  
  double dAmount = 0.;
  
  paymentDates.push_back(fixedCouponDate);
  paymentRates.push_back(dAmount);  
  FloatingInterests->SetKnownPaymentStream(paymentDates, paymentRates);
  
  bc->SetFloatingRates(FloatingInterests);

  // Cash Flow Stream
  paymentDates.clear();
  paymentRates.clear();

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral
                    ( issueDate, paymentDates, paymentRates, 
                      Date::DayCountConvention_30360,
                      Frequency_Annual)
              );
    
  bc->SetCashDistribution(pInterests);
}

void BondLikeTest::BondTermsCashFlowStreamThenFloatingRates()
{
  Date
    issueDate = Date(2002, Date::Jan, 14),
    maturityDate = Date(2004, Date::Jan, 14);

  double
    dIssuePrice = 1,
    dParValue = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                       maturityDate, dParValue, dRedemptionRate,
                       dRecoveryRate)
      );
  
  // Cash Flow Stream

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral
                    ( issueDate, paymentDates, paymentRates,
                      Date::DayCountConvention_30360,
                      Frequency_Annual) 
              );
    
  bc->SetCashDistribution(pInterests);
  
  // Floating rates
 
  Date lastPaymentDate = maturityDate;
  Date firstPaymentDate = maturityDate;
  Frequency paymentFrequency = Frequency_Quarterly;
  const int iInterval = 12 / paymentFrequency; 

  firstPaymentDate = AddMonthsAdjustedForEndOfMonth(firstPaymentDate, 
    - 3 * iInterval);
  
  shared_ptr<FloatingRates> FloatingInterests(
    new FloatingRates(0., issueDate, firstPaymentDate, lastPaymentDate,
          paymentFrequency) );

  // Fixed coupon added
  paymentDates.clear();
  paymentRates.clear();

  Date 
    fixedCouponDate = AddMonthsAdjustedForEndOfMonth(firstPaymentDate, 
      - iInterval);
  
  double dAmount = 0.;
  
  paymentDates.push_back(fixedCouponDate);
  paymentRates.push_back(dAmount);  
  FloatingInterests->SetKnownPaymentStream(paymentDates, paymentRates);
  
  bc->SetFloatingRates(FloatingInterests);

}

void BondLikeTest::CC_InvalidFixedFXRate()
{  
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);
  
  pCB->SetFixedFXRate(-0.1);
}

void BondLikeTest::FQ_InvalidFXVolatility_High()
{  
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);
  
  pCB->SetFixedFXRate(.8);
  pCB->SetFixedQuanto(-0.2, 5.1);
}

void BondLikeTest::FQ_InvalidFXVolatility_Low()
{  
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);
  
  pCB->SetFixedFXRate(.8);
  pCB->SetFixedQuanto(-0.2, -0.1);
}

void BondLikeTest::FQ_InvalidCorrelationBetweenShareAndFXRate_High()
{  
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);
  
  pCB->SetFixedFXRate(0.8);
  pCB->SetFixedQuanto(0.2, 1.1);
}

void BondLikeTest::FQ_InvalidCorrelationBetweenShareAndFXRate_Low()
{  
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);
  
  pCB->SetFixedFXRate(0.8);
  pCB->SetFixedQuanto(0.2, -1.1);
}
