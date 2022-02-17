/////////////////////////////////////////////////////////////////////////////
// Name:        tests/ihg/testcrosscurrency.cpp
// Purpose:     testing cross-currency
// Author:      Nabil Ouachani
// Created:     2006/03/02
// RCS-ID:      $Id: testcrosscurrency.cpp,v 1.5 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
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
#include "ito33/finance/bondlike/cboption.h"
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
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/yieldcurve_flat.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/testbondlike.h"

#include "ito33/tests/utilexml.h"

#include "ihg/tests/testcrosscurrency.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCB);
ITO33_FORCE_LINK_MODULE(IHGPriceCBOption);

using namespace std;
using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;

shared_ptr<SessionData> CrossCurrencyTest::Init_CC_SessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01"); 

  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));
  
  shared_ptr<Equity> pEquity(new Equity(100, pNumeraire));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0) ) );
  
  // Setup the rate data, and attach to the session
  shared_ptr<RateData> pRateData(new RateData);
  shared_ptr<YieldCurve> pYC(new YieldCurveFlat(0.04));
  pRateData->SetYieldCurve(pNumeraire, pYC);
 
  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));
  
  // Spot FX rate
  shared_ptr<Numeraire> pNumeraireUSD(new Numeraire("USD"));
  double dSpotFXRate = 0.8;
  shared_ptr<SpotFXRates> pSpotFXRates( new SpotFXRates() );
  pSpotFXRates->SetFXRate(pNumeraire, pNumeraireUSD, dSpotFXRate);
  
  pSessionData->GetRateData()->SetSpotFXRates(pSpotFXRates);

  return pSessionData;
}

shared_ptr<finance::ConvertibleBond> 
CrossCurrencyTest::InitCrossCurrencyCB(const shared_ptr<SessionData>& pSessionData)
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

  // COCO
  shared_ptr<ConversionPeriod>
    pConPeriod(new ConversionPeriod(issueDate, maturityDate, 1));
  
  pConPeriod->SetCoCo(1., CoCoType_CheckAnyTimeAndConvertOnCheckDate, 1.,1.);

  conv->AddConversionPeriod( pConPeriod );

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

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

  // cross currency
  shared_ptr<Numeraire> pNumeraire(new Numeraire("USD"));
  shared_ptr<YieldCurve> pYC(new YieldCurveFlat(0.02));
  pSessionData->GetRateData()->SetYieldCurve(pNumeraire, pYC);

  pConvertibleBond->SetNumeraire(pNumeraire);

  return pConvertibleBond;
}

shared_ptr<finance::ConvertibleBond> 
CrossCurrencyTest::InitCB(const shared_ptr<SessionData>& pSessionData)
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

  // COCO
  shared_ptr<ConversionPeriod>
    pConPeriod(new ConversionPeriod(issueDate, maturityDate, 1));
  
  pConPeriod->SetCoCo(1., CoCoType_CheckAnyTimeAndConvertOnCheckDate, 1.,1.);

  conv->AddConversionPeriod( pConPeriod );

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

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

void CrossCurrencyTest::CC_TriggerInCurrencyOfUnderlyingButNoFixedFXRate()
{ 
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);  

  pCB->SetTriggerInCurrencyOf(TriggerInCurrencyOf_Underlying);  

  // Model
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  double dLambda = .2;
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  // Pricing
  shared_ptr<ModelOutput> output = pModel->Compute(*pCB);
}

void CrossCurrencyTest::FQ_ButNotCrossCurrencyInstrument()
{ 
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCB(pSessionData);   

  // Sets the Quanto parameters
  pCB->SetFixedQuanto( 0.2, 0.1);
  pCB->SetFixedFXRate(0.8);

  // Model
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  double dLambda = .2;
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  // Pricing
  shared_ptr<ModelOutput> output = pModel->Compute(*pCB);
}

void CrossCurrencyTest::FQ_NeededFixedFXRateMissing()
{ 
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);  

  // Sets the Quanto parameters
  pCB->SetFixedQuanto( 0.2, 0.1);  

  // Model
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  double dLambda = .2;
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  // Pricing
  shared_ptr<ModelOutput> output = pModel->Compute(*pCB);
}

void CrossCurrencyTest::CBOption_With_FQ_NotSupported()
{ 
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);  

  // Sets the Quanto parameters
  pCB->SetFixedQuanto( 0.2, 0.1);
  pCB->SetFixedFXRate(0.8);

  // Create the floating
  double 
    dRecallSpread = 0.01;

  Date 
    maturityDate(pCB->GetMaturityDate()),
    startOfAccruedDate("2003/02/01"), 
    firstUnknownCouponDate("2003/08/01"),
    lastButOneUnknownCouponDate("2004/02/01");
  
  Date::DayCountConvention dcc = Date::DayCountConvention_Act365;
  
  Frequency frequency = Frequency_Quarterly;

  shared_ptr<finance::FloatingRates>
    pFloatingRates ( new finance::FloatingRates
                                  (
                                    dRecallSpread, 
                                    startOfAccruedDate, 
                                    firstUnknownCouponDate, 
                                    maturityDate, frequency
                                  ) 
                   );
  
  pFloatingRates->SetDayCountConvention(dcc);

  vector<Date> pknownPaymentDates;
  vector<double> pknownPaymentAmounts;
  pknownPaymentDates.push_back("2003/05/01");
  pknownPaymentAmounts.push_back(0.01);

  pFloatingRates->SetKnownPaymentStream(pknownPaymentDates, 
    pknownPaymentAmounts);

  // Create the CBOption
  shared_ptr<finance::CBOption> 
    pCBOption ( new finance::CBOption(pCB,
                                      pFloatingRates,
                                      maturityDate) 
              );
  
  pCBOption->SetSessionData(pSessionData);

  // Model
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  double dLambda = .2;
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  // Pricing
  shared_ptr<ModelOutput> output = pModel->Compute(*pCBOption);
}

void CrossCurrencyTest::CBOption_And_CB_Not_Same_Currency()
{ 
  shared_ptr<SessionData> pSessionData = Init_CC_SessionData();

  shared_ptr<finance::ConvertibleBond> pCB = InitCrossCurrencyCB(pSessionData);

  // Create the floating
  double 
    dRecallSpread = 0.01;

  Date 
    maturityDate(pCB->GetMaturityDate()),
    startOfAccruedDate("2003/02/01"), 
    firstUnknownCouponDate("2003/08/01"),
    lastButOneUnknownCouponDate("2004/02/01");
  
  Date::DayCountConvention dcc = Date::DayCountConvention_Act365;
  
  Frequency frequency = Frequency_Quarterly;

  shared_ptr<finance::FloatingRates>
    pFloatingRates ( new finance::FloatingRates
                                  (
                                    dRecallSpread, 
                                    startOfAccruedDate, 
                                    firstUnknownCouponDate, 
                                    maturityDate, frequency
                                  ) 
                   );
  
  pFloatingRates->SetDayCountConvention(dcc);

  vector<Date> pknownPaymentDates;
  vector<double> pknownPaymentAmounts;
  pknownPaymentDates.push_back("2003/05/01");
  pknownPaymentAmounts.push_back(0.01);

  pFloatingRates->SetKnownPaymentStream(pknownPaymentDates, 
    pknownPaymentAmounts);

  // Create the CBOption
  shared_ptr<finance::CBOption> 
    pCBOption ( new finance::CBOption(pCB,
                                      pFloatingRates,
                                      maturityDate) 
              );
  
  pCBOption->SetSessionData(pSessionData);

  // Sets the numeraire of the cboption to the one of the equity, then, 
  // in this case, CBOption_currency != CB_currency
  pCBOption->SetNumeraire(pSessionData->GetEquity()->GetNumeraire());

  // Model
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  double dLambda = .2;
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  // Pricing
  shared_ptr<ModelOutput> output = pModel->Compute(*pCBOption);
}

