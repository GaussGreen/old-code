#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"

#include "ito33/finance/bondlike/pepsaveragingperiod.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callfixedshare.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ito33/link.h"

#include "ito33/tests/showconvergence.h"

ITO33_FORCE_LINK_MODULE(IHGPricePEPSLike);

ITO33_FORCE_LINK_MODULE(IHGPriceGeneralizedPEPSLike);

ITO33_FORCE_LINK_MODULE(IHGPricePERCSLike);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

shared_ptr<finance::PERCSLike> 
InitPERCS(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  double
    dIssuePrice = 1.,
    dNominal = 110.,
    dRecoveryRate = 0.;

  shared_ptr<BondLikeTerms> 
    pblt( new BondLikeTerms(issueDate, dIssuePrice,
                            maturityDate, dNominal,
                            dRecoveryRate)
        );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamGeneral
                    (issueDate, paymentDates, paymentRates,
                     Date::DayCountConvention_Act365,
                     Frequency_Annual)
              );

  // pblt->SetCashDistribution(pInterests);

  shared_ptr<PERCSLike> 
    pPERCS( new PERCSLike(pblt, dIssuePrice, 1.) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  shared_ptr<CallSchedule> pCall( new CallSchedule() );
  // pCall->SetKeepAccrued(false);
         
  pCall->AddCallPeriod( shared_ptr<CallPeriod>
            ( CallPeriod::CreateWithStrike(startCallDate, endCallDate, 1.) ) );//
          
  pPERCS->SetCallSchedule(pCall);

  pPERCS->SetSessionData(pSessionData);

  return pPERCS;
}

shared_ptr<finance::PEPSLike> 
InitPEPS(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2005, Date::May, 1);

  double
    dIssuePrice = 1.,
    dNominal = 110.,
    dRecoveryRate = 0.;

  shared_ptr<BondLikeTerms> 
    pblt( new BondLikeTerms(issueDate, dIssuePrice,
                            maturityDate, dNominal, dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamGeneral
                    (issueDate, paymentDates, paymentRates,
                     Date::DayCountConvention_Act365,
                     Frequency_Annual)
              );

  pblt->SetCashDistribution(pInterests);

  shared_ptr<PEPSLike> 
    pPEPS( new PEPSLike(pblt, 1.5, 0.8) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  pPEPS->EnableOptionalConversion();

  /*
  shared_ptr<CallSchedule> pCall( new CallSchedule() );
  // pCall->SetKeepAccrued(false);
    
  size_t nNoticePeriod = 20;

  //  pCall->SetNoticePeriod(nNoticePeriod);
  //  startCallDate.AddMonths(4);
  //  endCallDate.AddMonths(-4);
          
  pCall->AddCallPeriod( CallPeriod::CreateWithStrike
                                    (startCallDate, endCallDate, 1) );//
       
  pPEPS->SetCallFixedCash(pCall);

  */

  shared_ptr<CallFixedShare> 
    pCall( new CallFixedShare(issueDate, maturityDate, .8) );

  pCall->SetTrigger(1.3);
  pPEPS->SetCallFixedShare(pCall);

  pPEPS->SetSessionData(pSessionData);

  return pPEPS;
}

shared_ptr<finance::GeneralizedPEPSLike> 
InitGeneralizedPEPS(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2005, Date::May, 1);

  double
    dIssuePrice = 1.,
    dNominal = 110.,
    dRecoveryRate = 0.;

  shared_ptr<BondLikeTerms> 
    pblt( new BondLikeTerms(issueDate, dIssuePrice,
                             maturityDate, dNominal, dRecoveryRate) );
  
  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamGeneral
                    (issueDate, paymentDates, paymentRates,
                     Date::DayCountConvention_Act365,
                     Frequency_Annual)
              );

  pblt->SetCashDistribution(pInterests);

  shared_ptr<GeneralizedPEPSLike> 
    pPEPS( new GeneralizedPEPSLike(pblt, 1.5, 110/1.5, 0.8, 110/0.8) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  pPEPS->EnableOptionalConversion();

  /*
  shared_ptr<CallSchedule> pCall = new CallSchedule();
  // pCall->SetKeepAccrued(false);
    
  size_t nNoticePeriod = 20;

  //  pCall->SetNoticePeriod(nNoticePeriod);
  //  startCallDate.AddMonths(4);
  //  endCallDate.AddMonths(-4);
          
  pCall->AddCallPeriod( CallPeriod::CreateWithStrike
                                    (startCallDate, endCallDate, 1) );//
       
  pPEPS->SetCallFixedCash(pCall);

  */

  shared_ptr<GeneralizedPEPSLikeCall> 
    pCall( new GeneralizedPEPSLikeCall(issueDate, maturityDate,
                          1.3, GeneralizedPEPSLikeCallType_FixedShare) );

  pPEPS->SetGeneralizedPEPSLikeCall(pCall);

  pPEPS->SetSessionData(pSessionData);

  return pPEPS;
}

shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  Date valuationDate(2004, Date::Apr, 30);
 
  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(120, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
  
  // Setup the rate data, and attach to the session data
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04)); 
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;

}

shared_ptr<finance::GeneralizedPEPSLike> InitAces(Date valuationDate, 
                                                 Date maturityDate,
                                                 Date avgStartDate, 
                                                 Date avgEndDate,
                                                 size_t nNbSampling,
                                                 size_t nNbSamplesUsed,
                                                 double dCurrentStockAverage)
{
  // Setup the pricing machinery
 
  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("USD"));
  shared_ptr<Equity> pEquity(new Equity(40, pCurrency));
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
  
  shared_ptr<Dividends> pDividends( new Dividends() );
  pDividends->AddCash( Date(2004, Date::Mar, 1) , .9333);
  pDividends->AddCash( Date(2004, Date::Jun, 1) , .875 );
  pDividends->AddCash( Date(2004, Date::Sep, 1) , .875 );
  pDividends->AddCash( Date(2004, Date::Dec, 1) , .875 );
 
  pDividends->AddCash( Date(2005, Date::Mar, 1) , .875);
  pDividends->AddCash( Date(2005, Date::Jun, 1) , .875 );
  pDividends->AddCash( Date(2005, Date::Sep, 1) , .875 );
  pDividends->AddCash( Date(2005, Date::Dec, 1) , .875 );

  pDividends->AddCash( Date(2006, Date::Mar, 1) , .875);
  pDividends->AddCash( Date(2006, Date::Jun, 1) , .875 );
  pDividends->AddCash( Date(2006, Date::Sep, 1) , .875 );
  pDividends->AddCash( Date(2006, Date::Dec, 1) , .875 );

  pEquity->SetDividends(pDividends);

  // Setup the rate data, and attach to the session data
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04)); 
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );
  
  Date issueDate = Date(2003, Date::Sep, 30);
  
  double dIssuePrice   = 1.;
  double dNominal      = 14.15;
  double dRecoveryRate = 0.;

  shared_ptr<BondLikeTerms> 
    pblt( new BondLikeTerms(issueDate, dIssuePrice,
                             maturityDate, dNominal, dRecoveryRate) );

 double dLowerStrike  = 48.55;
 double dHigherStrike = 60.20;
 double dDownSideConversionRatio   =  1.0299;
 double dUpsideBaseConversionRatio = .8305;

 shared_ptr<GeneralizedPEPSLike> 
   pPEPS( new GeneralizedPEPSLike(pblt, dDownSideConversionRatio, 
                                  dLowerStrike, dUpsideBaseConversionRatio, 
                                  dHigherStrike) );

  pPEPS->EnableOptionalConversion();

  if ( avgStartDate < avgEndDate )
  {
    shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
      (
           PEPSAveragingPeriod::CreateWithStock(avgStartDate, 
                                                avgEndDate, 
                                                nNbSampling)
      );

    if ( dCurrentStockAverage )
     pAveragingPeriod->SetCurrentStockAverage(dCurrentStockAverage, 
                                              nNbSamplesUsed);

    pPEPS->SetAveragingPeriod(pAveragingPeriod);
  }

  //calling part
  shared_ptr<GeneralizedPEPSLikeCall> 
    pCall( new GeneralizedPEPSLikeCall(issueDate, maturityDate,
                          1.5, GeneralizedPEPSLikeCallType_FixedShare) );

  pPEPS->SetGeneralizedPEPSLikeCall(pCall);

  pPEPS->SetSessionData(pSessionData);

  return pPEPS;

}


shared_ptr<GeneralizedPEPSLike> InitSuez()
{
  // Setup the pricing machinery
  Date valuationDate(2006, Date::May, 20);
 
  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(14.15, pCurrency));
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
  
  // Setup the rate data, and attach to the session data
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04)); 
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );
  
  Date issueDate = Date(2003, Date::May, 22);
  Date maturityDate = Date(2006, Date::May, 22);

  double dIssuePrice   = 1.;
  double dNominal      = 50.;
  double dRecoveryRate = 0.;

  shared_ptr<BondLikeTerms> 
    pblt( new BondLikeTerms(issueDate, dIssuePrice,
                            maturityDate, dNominal, dRecoveryRate)
        );

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamUniform(issueDate,
                         Date(2004,Date::May,22),
                         maturityDate,
                         4.5/100.,
                         Date::DayCountConvention_Act360,
                         Frequency_Annual)
              );

  pblt->SetCashDistribution(pInterests);

  double dLowerStrike  = 17.0;
  double dHigherStrike = 20.0;
  double dDownSideConversionRatio   = 1.;
  double dUpsideBaseConversionRatio = 1.;

  shared_ptr<GeneralizedPEPSLike> 
    pPEPS( new GeneralizedPEPSLike(pblt, dDownSideConversionRatio, 
                                   dLowerStrike, dUpsideBaseConversionRatio, 
                                   dHigherStrike) );

  pPEPS->EnableOptionalConversion();

 
  Date avgEndDate = maturityDate;
  avgEndDate.AddDays(-3);
  Date avgStartDate = avgEndDate;
  avgStartDate.AddDays( -16);
  size_t nNbSampling = 15;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
    PEPSAveragingPeriod::CreateWithConversionRatio(avgStartDate, 
                                                     avgEndDate, 
                                                     nNbSampling)
    );

  pAveragingPeriod->SetCurrentConversionRatioAverage(.5, 2);

  pPEPS->SetAveragingPeriod(pAveragingPeriod);

  //calling part
  shared_ptr<GeneralizedPEPSLikeCall>
    pCall( new GeneralizedPEPSLikeCall(issueDate, maturityDate,
                          1.2, GeneralizedPEPSLikeCallType_FixedShare) );

  pPEPS->SetGeneralizedPEPSLikeCall(pCall);

  pPEPS->SetSessionData(pSessionData);

  return pPEPS;

}


void ExtraACESTest()
{

  shared_ptr<SessionData> pSessionData = InitSessionData();

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.3;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  pModel->SetHazardRate(shared_ptr<HazardRate>(new HazardRatePower(0.01, 1.2, 100)));

  // Try changing the number of averaging days
  Date valuationDate(2005, Date::Nov, 28);
  Date maturityDate(2006, Date::Dec, 1);

  Date avgEndDate = maturityDate; avgEndDate.AddDays(-3);
  Date avgStartDate  = avgEndDate; avgEndDate.AddDays(-31);

  size_t nNbSampling = Date::DaysDiff(avgStartDate, avgEndDate);
  double dCurrentStockAverage = -1;
  size_t nNbSamplesUsed = 5;


  std::cout << "Test length of averaging period" << std::endl << std::endl;

  for ( size_t nIdx = 0; nIdx <11; nIdx++)
  {
    avgStartDate.AddDays( nIdx  );

    if ( valuationDate > avgStartDate )
      dCurrentStockAverage = 20;

    shared_ptr<finance::GeneralizedPEPSLike> pAces = InitAces(valuationDate,
      maturityDate, avgStartDate, avgEndDate, 
      nNbSampling, nNbSamplesUsed, dCurrentStockAverage );

    shared_ptr<finance::ModelOutput> output;
  
    output = pModel->Compute(*pAces);

    double dPrice = output->GetPrice();

    
    std::cout << "Days: " << Date::DaysDiff(avgEndDate, avgStartDate) 
      << " dPrice: " << dPrice << std::endl;

  } // while looping over nNbDays
  std::cout << std::endl;
/*

  // Try changing the offset 
  valuationDate = Date(2005, Date::Nov, 28);
  nNbDays = 20;
  nOffSet = 5;
  dOldPrice = 0.0;
  std::cout << "Test length of offset" << std::endl << std::endl;
  while (nOffSet < 6)
  {
    shared_ptr<finance::GeneralizedPEPSLike> pAces = InitAces(valuationDate, 
      nNbDays, nOffSet, dCurrentStockAverage);

    shared_ptr<finance::ModelOutput> output;
  
    output = pModel->Compute(*pAces);

    double dPrice = output->GetPrice();

    if (dOldPrice > 0.0)
    {
      std::cout << "Offset: " << nOffSet 
                << ", price = " << dPrice
                << ", diff = " << dOldPrice - dPrice
                << std::endl;
    }
    else
    {
      std::cout << "Offset: " << nOffSet 
                << ", price = " << dPrice
                << std::endl;
    }

    dOldPrice = dPrice;

    nOffSet -= 1;

  } // while looping over nOffSet
  std::cout << std::endl;

  // Try changing the offset 
  valuationDate = Date(2005, Date::Nov, 28);
  nNbDays = 0;
  nOffSet = 5;
  dOldPrice = 0.0;
  std::cout << "Change offset length when averaging window length is zero" 
            << std::endl << std::endl;
  while (nOffSet < 6)
  {
    shared_ptr<finance::GeneralizedPEPSLike> pAces = InitAces(valuationDate, 
      nNbDays, nOffSet, dCurrentStockAverage);

    shared_ptr<finance::ModelOutput> output;
  
    output = pModel->Compute(*pAces);

    double dPrice = output->GetPrice();

    if (dOldPrice > 0.0)
    {
      std::cout << "Offset: " << nOffSet 
                << ", price = " << dPrice
                << ", diff = " << dOldPrice - dPrice
                << std::endl;
    }
    else
    {
      std::cout << "Offset: " << nOffSet 
                << ", price = " << dPrice
                << std::endl;
    }

    dOldPrice = dPrice;

    nOffSet -= 1;

    if (nOffSet > 6)
    {
      pAces = InitAces(valuationDate, nNbDays, 10000, dCurrentStockAverage);
      output = pModel->Compute(*pAces);

      double dPrice = output->GetPrice();

      std::cout << "Price without averaging window: " << dPrice << std::endl;

    }

  } // while looping over nOffSet
  std::cout << std::endl;


  // Try changing the current average
  valuationDate = Date(2006, Date::Nov, 20);
  nNbDays = 20;
  nOffSet = 3;
  dCurrentStockAverage = 30.0;
  dOldPrice = 0.0;
  std::cout << "Test current average" << std::endl << std::endl;
  while (dCurrentStockAverage < 51.0)
  {
    shared_ptr<finance::GeneralizedPEPSLike> pAces = InitAces(valuationDate, 
      nNbDays, nOffSet, dCurrentStockAverage);

    shared_ptr<finance::ModelOutput> output;
  
    output = pModel->Compute(*pAces);

    double dPrice = output->GetPrice();

    if (dOldPrice > 0.0)
    {
      std::cout << "Current average: " << dCurrentStockAverage 
                << ", price = " << dPrice
                << ", diff = " << dOldPrice - dPrice
                << std::endl;
    }
    else
    {
      std::cout << "Current average: " << dCurrentStockAverage
                << ", price = " << dPrice
                << std::endl;
    }

    dOldPrice = dPrice;

    dCurrentStockAverage += 2.0;

  } // while looping over current average

*/
}


int main()
{

 try
 {
   bool
     bComputeRho = false,
     bComputeVega = false,
     bComputeSurface = false;
   
   shared_ptr<SessionData> pSessionData = InitSessionData();

   shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

   double dVol = 0.5;
   pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

   size_t nNbHRTimes = 2;
   static const Date pTimes[2] = { 38000, 40000};
   static const double pdValues[2] = { 0.02, .02 };
  
   pModel->SetHazardRate
           // ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );
          (shared_ptr<HazardRate>
           ( new HazardRateCombo
                 (shared_ptr<SpotComponent>(new HRSpotComponentPower(0, 100)), 
                  &pTimes[0], &pdValues[0], nNbHRTimes)) );
   
/*
   double dLambda = .02;
   pModel->SetHazardRate( new ihg::HazardRateFlat(dLambda) );
*/
   //pModel->SetDebugOutputFile("./ihg_percs.xml");

   shared_ptr<finance::ModelOutput> output;

   // percs
   shared_ptr<finance::PERCSLike> pPERCS = InitPERCS(pSessionData);

   shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
   flags->SetComputeRho(bComputeRho);
   flags->SetComputeVega(bComputeVega);
   flags->SetComputeSurface(bComputeSurface);
   flags->SetAnalysisDate( Date("2003/02/01"));

   pPERCS->SetComputationalFlags(flags);   

   output = pModel->Compute(*pPERCS);

   std::cout.precision(15);
    
   std::cout << "The percs price is: " << output->GetPrice() << std::endl;
   if (output->HasRho())
      std::cout << "The percs Rho is: " << output->GetRho() << std::endl;
   if (output->HasVega())
      std::cout << "The percs Vega is: " << output->GetVega() << std::endl;


   // peps
   shared_ptr<finance::PEPSLike> pPEPS = InitPEPS(pSessionData);
   pModel->SetDebugOutputFile("./ihg_peps.xml");
   output = pModel->Compute(*pPEPS);

   std::cout.precision(15);
    
   std::cout << "The peps price is: " << output->GetPrice() << std::endl;
   if (output->HasRho())
      std::cout << "The peps Rho is: " << output->GetRho() << std::endl;
   if (output->HasVega())
      std::cout << "The peps Vega is: " << output->GetVega() << std::endl;


   // generalized peps
   shared_ptr<finance::GeneralizedPEPSLike> 
     pNewPEPS = InitGeneralizedPEPS(pSessionData);
   pModel->SetDebugOutputFile("./ihg_newpeps.xml");
   output = pModel->Compute(*pNewPEPS);

   std::cout.precision(15);
    
   std::cout << "The generalized peps price is: " << output->GetPrice() << std::endl;
   if (output->HasRho())
      std::cout << "The generalized peps Rho is: " << output->GetRho() << std::endl;
   if (output->HasVega())
      std::cout << "The generalized peps Vega is: " << output->GetVega() << std::endl;

   // generalized peps with Averaging
   Date valuationDate(2005, Date::May, 12);
   Date maturityDate(2006, Date::Dec, 1);

  Date avgEndDate = maturityDate; avgEndDate.AddDays(-3);
  Date avgStartDate  = avgEndDate; avgEndDate.AddDays(-31);

  size_t nNbSampling = Date::DaysDiff(avgStartDate, avgEndDate);
  double dCurrentStockAverage = -1;
  size_t nNbSamplesUsed = 5;
  

   shared_ptr<finance::GeneralizedPEPSLike> pAces = InitAces(valuationDate, maturityDate,
     avgStartDate, avgEndDate, nNbSampling, nNbSamplesUsed, dCurrentStockAverage);

   
  // pModel->SetDebugOutputFile("./ihg_aces.xml");
   
   output = pModel->Compute(*pAces);

   std::cout.precision(15); 
   std::cout << "The ACES peps price  is: " << output->GetPrice() << std::endl;
   
   if (output->HasRho())
      std::cout << "The ACES peps Rho    is: " << output->GetRho() << std::endl;
   if (output->HasVega())
      std::cout << "The ACES peps Vega   is: " << output->GetVega() << std::endl;

   //ShowConvergence(*pModel, *pAces, 4);

   // generalized peps with Averaging
   shared_ptr<finance::GeneralizedPEPSLike> pSuez = InitSuez();
   //pModel->SetDebugOutputFile("./ihg_suez.xml");
   output = pModel->Compute(*pSuez);

   std::cout.precision(15); 
   std::cout << "The Suez peps price  is: " << output->GetPrice() << std::endl;
   
   if (output->HasRho())
      std::cout << "The Suez peps Rho    is: " << output->GetRho() << std::endl;
   if (output->HasVega())
      std::cout << "The Suez peps Vega   is: " << output->GetVega() << std::endl;

   //ShowConvergence(*pModel, *pSuez, 4);

  // ExtraACESTest();


   return 0;
 }
 catch ( const ::ito33::Exception& e )
 {
   std::cerr << "Exception caught:\n"
     << e.GetFullMessage() << std::endl;

   return 1;
 }
 catch ( ... )
 {
   std::cerr << "Unexpected exception caught.\n";

   return 2;
 }
}
