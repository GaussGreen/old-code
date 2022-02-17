
#include "utils.h"

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/frequency.h"

#include "ito33/finance/bondlike/cocotype.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

using namespace ito33::finance;

namespace ito33
{

shared_ptr<CallSchedule> GetCallSchedule(Date callStart,Date callEnd)
{
  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
    
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(callStart, callEnd, 1.0) );
  pCallSchedule->AddCallPeriod( pCallPeriod );
        
  return pCallSchedule;
}


shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01");

  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
 
  shared_ptr<YieldCurve> pyc (new YieldCurveFlat(0.04)); 
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<BondTerms> InitBondTerms()
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
                      dRecoveryRate)
      );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral(issueDate, paymentDates, 
                     paymentRates, Date::DayCountConvention_Act365,
                     Frequency_Annual) );
  
  bc->SetCashDistribution(pInterests);

  return bc;
}

shared_ptr<ConvertibleBond> InitCB(const shared_ptr<SessionData>& pSessionData)
{
  shared_ptr<BondTerms> bc = InitBondTerms();

  Date conversionStart = Date(2003, Date::Apr, 1);
  Date conversionEnd = Date(2004, Date::May, 1);

  Date callStart( bc->GetIssueDate() );
  Date callEnd( bc->GetMaturityDate() );

  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(true);
  conv->SetForfeitCoupon(false);

  conv->AddConversionPeriod(shared_ptr<ConversionPeriod>
     (new ConversionPeriod(conversionStart, conversionEnd, 1)));

  shared_ptr<ConvertibleBond> 
    pConvertibleBond( new ConvertibleBond(bc, conv) );

  pConvertibleBond->SetCallSchedule( GetCallSchedule(callStart,callEnd) );

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}


shared_ptr<AttachedWarrantConvertibleBond> 
InitAttachedWarrantCB(const shared_ptr<SessionData>& pSessionData,
                      double dShareFactor,
                      double dStrikeIn,
                      double dCap)
{

  shared_ptr<BondTerms> bc = InitBondTerms();

  Date conversionStart = Date(2003, Date::Apr, 1);
  Date conversionEnd = Date(2004, Date::May, 1);

  Date resetDate = Date(2004, Date::Jan, 1);
  double dInitialRatio = 1.0;

  // If the strike was not specified, set to the spot
  dStrikeIn = 90;//dStrikeIn;
  //if (dStrikeIn < 0.0)
  //  dStrike = pSession->GetSpotSharePrice();  

  Date callStart( bc->GetIssueDate() );
  Date callEnd( bc->GetMaturityDate() );
  
  shared_ptr<ShareDependentConversion> 
    pConv( new ShareDependentConversion(conversionStart, conversionEnd, 
                 dInitialRatio, dShareFactor) );
  pConv->SetCapRatio(dCap);
  pConv->SetResetDate(resetDate);
  
  if ( dStrikeIn >= 0 )
    pConv->SetFixedStrike(dStrikeIn);

  pConv->SetKeepAccrued(true);
  pConv->SetForfeitCoupon(false);

  shared_ptr<AttachedWarrantConvertibleBond> 
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  pAttachedWarrantCB->SetCallSchedule(GetCallSchedule(callStart,callEnd));


  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
}


shared_ptr<AttachedWarrantConvertibleBond>
InitAttachedWarrantCBNoReset(double dShareFactor, double dCap)
{
  // Construct test for Carnival_1.132_2033.pdf
  // The session variables
  Date valuationDate = Date(2005,Date::Feb,22);

  double dSpot = 53.98; //current stock price
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );

  shared_ptr<YieldCurve> pyc (new YieldCurveFlat(0.04));
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
  pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  // Setup the bond terms
  Date issueDate = Date(2004, Date::Mar, 19);
  Date maturityDate = Date(2033, Date::Dec, 1);

  double dIssuePrice     = 1.; //issue price percentage of nominal
  double dParValue       = 1000.;
  double dRedemptionRate = 1.64412;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
                      dParValue, dRedemptionRate, dRecoveryRate)
      );

  Date firstCashDate = Date(2004, Date::Jun, 1);
  Date lastCashDate = Date(2006, Date::Dec, 1);

  double dCashYield = 1.85/100.;
  shared_ptr<CashFlowStreamUniform>
      pInterests( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstCashDate,
                          lastCashDate,
                          dCashYield,
                          Date::DayCountConvention_Act365,
                          Frequency_SemiAnnual
                        )
                  );

  bc->SetCashDistribution(pInterests);

  Date conversionStart = bc->GetIssueDate();
  Date conversionEnd = bc->GetMaturityDate();

  double dInitialRatio = 1.0;

  Date callStart( bc->GetIssueDate() );
  Date callEnd( bc->GetMaturityDate() );

  Date resetDate = Date(2008, Date::Mar, 01);

  shared_ptr<ShareDependentConversion>
    pConv( new ShareDependentConversion(conversionStart, conversionEnd, 
                dInitialRatio, dShareFactor) );

  pConv->SetCapRatio(dCap);
  pConv->SetKeepAccrued(true);
  pConv->SetForfeitCoupon(false);

  shared_ptr<AttachedWarrantConvertibleBond>
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  pAttachedWarrantCB->SetCallSchedule(GetCallSchedule(callStart,callEnd));

  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
} 


shared_ptr<AttachedWarrantConvertibleBond> 
InitCarnivalTest(bool bShiftCallDate)
{

  // Construct test for Carnival_1.132_2033.pdf
  // The session variables
  Date valuationDate("2004/01/01");

  double dSpot = 646.88 / 12.18;
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));
  
/*
  // To get a price difference with coco, need dividends and small
  // incremental share factor
  shared_ptr<Dividends> pDividends(new Dividends() );
  pDividends->AddCash(Date(2006, Date::May, 01), 5.0);
  pEquity->SetDividends(pDividends);
*/

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );

  shared_ptr<YieldCurve> pyc (new YieldCurveFlat(0.04)); 
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  // Setup the bond terms
  Date issueDate = Date(2003, Date::Apr, 29);
  Date maturityDate = Date(2033, Date::Apr, 29);

  double dIssuePrice     = 646.88/1000; //issue price percentage of nominal
  double dParValue       = 1000;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
                      dParValue, dRedemptionRate, dRecoveryRate)
      );

  Date firstCashDate = Date(2003, Date::Oct, 29);
  Date lastCashDate = Date(2008, Date::Apr, 29);

  double dCashYield = 1.132/100.;
  shared_ptr<CashFlowStreamUniform> 
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate,
                      firstCashDate, 
                      lastCashDate,
                      dCashYield,
                      Date::DayCountConvention_30360,
                      Frequency_SemiAnnual
                    )
              );

  bc->SetCashDistribution(pInterests);

  double dOIDYield = 1.75/100.;
  bc->SetCashPayToZero(dOIDYield);
  bc->SetYieldCompoundingFrequency(Frequency_SemiAnnual);

  // Conversion provisions
  Date conversionStart = Date(2003,Date::Aug,31);
  Date conversionEnd = maturityDate;
  Date resetDate = Date(2008, Date::Apr, 29);

  double dInitialRatio = 12.18;
  double dShareFactor  = 11.3258;
  double dStrike = dIssuePrice * dParValue / dInitialRatio;
  
  shared_ptr<ShareDependentConversion> 
    pConv( new ShareDependentConversion
               (conversionStart, conversionEnd, dInitialRatio, dShareFactor)
         );

  pConv->SetResetDate(resetDate);
  pConv->SetFixedStrike(dStrike);

  pConv->SetKeepAccrued(false);
  pConv->SetForfeitCoupon(false);
   
  double dTriggerRate = 120./100.;
  double dChangeRate  = 0.0;
  double dExtremeTriggerRate = dTriggerRate;
  bool bIsCurrentlyActive = true;
  pConv->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertDuringNextQuarter,
                 dChangeRate, dExtremeTriggerRate, bIsCurrentlyActive);

  // Call schedule
  Date callStart = Date(2008, Date::Apr, 29);
  Date callEnd = maturityDate;

  if ( bShiftCallDate )
    callStart.AddDays(30);

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
    
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(callStart, callEnd, 1.) );
  pCallSchedule->AddCallPeriod( pCallPeriod );

  // Put provisions
  shared_ptr<PutSchedule> pPutSchedule(new PutSchedule);
  pPutSchedule->AddPutWithStrike(Date(2008, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2013, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2018, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2023, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2028, Date::Apr, 29), 1.0);

  // Construct the attached warrant cb
  shared_ptr<AttachedWarrantConvertibleBond> 
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  // Finalize the contract (call schedule, session, etc)
  pAttachedWarrantCB->SetCallSchedule(pCallSchedule);
  pAttachedWarrantCB->SetPutSchedule(pPutSchedule);
  pAttachedWarrantCB->SetConversionTriggerAsPercentageOf(TriggerAsPercentageOf_IssuePrice);
  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
}


shared_ptr<AttachedWarrantConvertibleBond> 
InitAmericanExpressTest(bool bShiftCallDate)
{
  // Construct test for Carnival_1.132_2033.pdf
  // The session variables
  Date valuationDate = Date(2005,Date::Feb,22);

  double dSpot = 53.98; //current stock price
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );

  shared_ptr<YieldCurve> pyc (new YieldCurveFlat(0.04)); 
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  // Setup the bond terms
  Date issueDate = Date(2004, Date::Mar, 19);
  Date maturityDate = Date(2033, Date::Dec, 1);

  double dIssuePrice     = 1.; //issue price percentage of nominal
  double dParValue       = 1000.;
  double dRedemptionRate = 1.64412; 
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
                      dParValue, dRedemptionRate, dRecoveryRate)
      );

  Date firstCashDate = Date(2004, Date::Jun, 1);
  Date lastCashDate = Date(2006, Date::Dec, 1);

  double dCashYield = 1.85/100.;
  shared_ptr<CashFlowStreamUniform> 
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate,
                      firstCashDate, 
                      lastCashDate,
                      dCashYield,
                      Date::DayCountConvention_Act365,
                      Frequency_SemiAnnual
                    )
              );

  bc->SetCashDistribution(pInterests);

  double dOIDYield = 1.85/100.;
  bc->SetCashPayToZero(dOIDYield);
  bc->SetYieldCompoundingFrequency(Frequency_SemiAnnual);

  // Conversion provisions
  Date conversionStart = issueDate;
  Date conversionEnd   = maturityDate;
  Date resetDate       = Date(2006, Date::Dec, 1);

  double dBaseRatio    = 14.4073;
  double dShareFactor  = 43.2219;
  double dStrike       = dIssuePrice * dParValue / dBaseRatio;
  double dCap          = 22.7635; 
  
  shared_ptr<ShareDependentConversion> 
    pConv( new ShareDependentConversion
               (conversionStart, conversionEnd, dBaseRatio, dShareFactor)
         );

  pConv->SetCapRatio(dCap);
  pConv->SetResetDate(resetDate);
  pConv->SetFixedStrike(dStrike);

  pConv->SetKeepAccrued(true);
  pConv->SetForfeitCoupon(false); 
   
  double dTriggerRate = 125./100.;
  double dChangeRate  = 0.0;
  double dExtremeTriggerRate = dTriggerRate;
  bool bIsCurrentlyActive = false;
  pConv->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertAsOfCheckDate,
                 dChangeRate, dExtremeTriggerRate, bIsCurrentlyActive);

  // Call schedule
  Date callStart = Date(2006, Date::Dec, 1);
  Date callEnd   = maturityDate;

  if ( bShiftCallDate )
    callStart.AddDays(1);

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
    
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(callStart, callEnd, 1.0) );
  pCallSchedule->AddCallPeriod( pCallPeriod );

  // Put provisions
  shared_ptr<PutSchedule> pPutSchedule(new PutSchedule);
  pPutSchedule->AddPutWithStrike(Date(2006, Date::Dec, 1), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2008, Date::Dec, 1), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2013, Date::Dec, 1), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2018, Date::Dec, 1), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2023, Date::Dec, 1), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2028, Date::Dec, 1), 1.0);

  // Construct the attached warrant cb
  shared_ptr<AttachedWarrantConvertibleBond> 
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  // Finalize the contract (call schedule, session, etc)
  pAttachedWarrantCB->SetCallSchedule(pCallSchedule);
  pAttachedWarrantCB->SetPutSchedule(pPutSchedule);
  pAttachedWarrantCB->SetConversionTriggerAsPercentageOf
                      ( TriggerAsPercentageOf_IssuePrice );
  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
}


shared_ptr<AttachedWarrantConvertibleBond> 
InitAutobacsTest()
{
  // Construct test for Autobacs 0% 2023.pdf prospectus
  // The session variables
  Date valuationDate("2004/01/01");

  // Set the spot to the base conversion price
  double dSpot = 3220;
  shared_ptr<Numeraire> pCurrency(new Numeraire("YEN"));
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));
  
  // See page 16
  shared_ptr<Dividends> pDividends(new Dividends() );
  for (size_t nYear = 2003; nYear <= 2023; nYear++)
  {
    pDividends->AddCash(Date(nYear, Date::Mar, 31), 18.0);
    pDividends->AddCash(Date(nYear, Date::Sep, 30), 18.0);
  }
  pEquity->SetDividends(pDividends);

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );

  shared_ptr<YieldCurve> pyc (new YieldCurveFlat(0.04)); 
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  // Setup the bond terms
  Date issueDate = Date(2003, Date::Sep, 22);
  Date maturityDate = Date(2023, Date::Sep, 30);

  double dParValue       = 5000000;
  double dIssuePrice     = 1.;   
  double dRedemptionRate = 1.;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
                      dParValue, dRedemptionRate, dRecoveryRate)
      );
    
  // Conversion provisions
  Date conversionStart = Date(2003,Date::Oct, 22);
  Date conversionEnd = Date(2023, Date::Sep, 15);

  double dInitialRatio = 1552.79;
  double dShareFactor  = 2173.91;
  double dStrike = 3220;
  double dCap = 1976.28; 
  
  shared_ptr<ShareDependentConversion> 
    pConv( new ShareDependentConversion
               (conversionStart, conversionEnd, dInitialRatio, dShareFactor)
         );
  pConv->SetCapRatio(dCap);
  pConv->SetFixedStrike(dStrike);

  // These should not matter, since there are no coupons. Left here
  // for testing purposes
  pConv->SetKeepAccrued(false);
  pConv->SetForfeitCoupon(false);
   
  // Coco feature
  double dTriggerRate = 110./100.;
  double dChangeRate  = 0.0;
  double dExtremeTriggerRate = dTriggerRate;
  bool bIsCurrentlyActive = true;
  pConv->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertDuringNextQuarter,
                 dChangeRate, dExtremeTriggerRate, bIsCurrentlyActive);

  // Call schedule
  Date callStart = Date(2007, Date::Sep, 30);
  Date callEnd = maturityDate;

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  
  // Again, these should not matter
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
    
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(callStart, callEnd, 1.0) );
  pCallSchedule->AddCallPeriod( pCallPeriod );

  // Put provisions
  shared_ptr<PutSchedule> pPutSchedule(new PutSchedule);
  pPutSchedule->AddPutWithStrike(Date(2007, Date::Sep, 30), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2011, Date::Sep, 30), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2015, Date::Sep, 30), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2019, Date::Sep, 30), 1.0);

  // Construct the attached warrant cb
  shared_ptr<AttachedWarrantConvertibleBond> 
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  // Finalize the contract (call schedule, session, etc)
  pAttachedWarrantCB->SetCallSchedule(pCallSchedule);
  pAttachedWarrantCB->SetPutSchedule(pPutSchedule);
  pAttachedWarrantCB->SetConversionTriggerAsPercentageOf(TriggerAsPercentageOf_IssuePrice);
  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
}

shared_ptr<finance::AttachedWarrantConvertibleBond> InitGettyTest()
{

  // Construct test for gyi.pdf
  // The session variables
  Date valuationDate("2005/07/12");

  double dSpot = 61.08;
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));  

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );

  shared_ptr<YieldCurve> pyc (new YieldCurveFlat(0.04)); 
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  // Setup the bond terms
  Date issueDate    = Date(2003, Date::Jun, 9);
  Date maturityDate = Date(2023, Date::Jun, 9);

  double dIssuePrice     = 1.0; //issue price percentage of nominal
  double dParValue       = 1000;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
                      dParValue, dRedemptionRate, dRecoveryRate)
      );

  Date firstCashDate = Date(2003, Date::Dec, 9);
  Date lastCashDate  = maturityDate;

  double dCashYield = .5/100.;
  shared_ptr<CashFlowStreamUniform> 
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate,
                      firstCashDate, 
                      lastCashDate,
                      dCashYield,
                      Date::DayCountConvention_30360,
                      Frequency_SemiAnnual
                    )
              );

  bc->SetCashDistribution(pInterests);

  // Conversion provisions
  Date conversionStart = Date(2003, Date::Jun, 9);
  Date conversionEnd   = maturityDate;
  Date resetDate       = Date(2008, Date::Jun, 9);

  double dInitialRatio = 16.372;
  double dShareFactor  = 16.372;
  double dStrike = dIssuePrice * dParValue / dInitialRatio;
  
  shared_ptr<ShareDependentConversion> 
    pConv( new ShareDependentConversion
               (conversionStart, conversionEnd, dInitialRatio, dShareFactor)
         );

  pConv->SetResetDate(resetDate);
  pConv->SetFixedStrike(dStrike);

  pConv->SetKeepAccrued(false);
  pConv->SetForfeitCoupon(false);
   
  double dTriggerRate = 120./100.;
  double dChangeRate  = 0.0;
  double dExtremeTriggerRate = dTriggerRate;
  bool bIsCurrentlyActive = true;
  pConv->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertDuringNextQuarter,
                 dChangeRate, dExtremeTriggerRate, bIsCurrentlyActive);

  // Call schedule
  Date callStart = Date(2008, Date::Jun, 13);
  Date callEnd   = maturityDate;


  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
    
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(callStart, callEnd, 1.) );
  pCallSchedule->AddCallPeriod( pCallPeriod );

  // Construct the attached warrant cb
  shared_ptr<AttachedWarrantConvertibleBond> 
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  // Finalize the contract (call schedule, session, etc)
  pAttachedWarrantCB->SetCallSchedule(pCallSchedule);
  pAttachedWarrantCB->SetConversionTriggerAsPercentageOf(TriggerAsPercentageOf_IssuePrice);
  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
}


} // ito33 namespace
