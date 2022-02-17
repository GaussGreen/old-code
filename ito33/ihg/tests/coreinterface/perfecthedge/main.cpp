#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/perfect_hedge_ratios.h"

#include "ito33/tests/showconvergence.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCB);
ITO33_FORCE_LINK_MODULE(IHGPriceCDS);
ITO33_FORCE_LINK_MODULE(IHGPriceOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

double AnalyticCDSPrice(const CDS& cds, double dRate, double dLambda, 
                        Date pricingDate); 

double functI(double, double, size_t i)
{
 return i / (double)100;
}

void
quickfunctI(double , const double *,  double *pdv, size_t nNb, int i)
{
 double d = functI(0, 0, i);

 for(size_t nIdx = 0; nIdx < nNb; nIdx++)
   pdv[nIdx] = d;
}



void TestConvergence(ihg::TheoreticalModel&, const finance::Derivative&);

shared_ptr<TermStructureCDS> 
MakeCDSList(const shared_ptr<SessionData>& pSessionData)
{
  shared_ptr<TermStructureCDS> tsCDS( new TermStructureCDS() );

  size_t nNbCDS = 5, nIdx;
  
  for (nIdx = 0; nIdx < nNbCDS; nIdx++)
  {
    Date IssueDate(2003, ito33::Date::Dec, 15);
    Date FirstDate(2004, ito33::Date::Jan, 1);
    
    Date MaturityDate = FirstDate;
    MaturityDate.AddMonths( (nIdx + 2) * 6);

    shared_ptr<CashFlowStreamUniform> 
      pSpreadStream( new CashFlowStreamUniform
                         (
                           IssueDate,
                           FirstDate, 
                           MaturityDate,
                           0.02,
                           Date::DayCountConvention_Act365,
                           Frequency_BiMonthly
                         )
                    );

    shared_ptr<CDS> pCDS(new CDS(0.4, pSpreadStream) );
    pCDS->SetSessionData(pSessionData);
    pCDS->SetMarketPrice(0);
    tsCDS->Add(pCDS);
  } 

  return tsCDS;
}

shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
    conversionStart = Date(2003, Date::Feb, 1),
    conversionEnd = Date(2004, Date::May, 1);

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
    pInterests( new CashFlowStreamGeneral(issueDate, paymentDates, paymentRates,
                Date::DayCountConvention_Act365, Frequency_Annual) );
  
  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  conv->AddConversionPeriod
        ( shared_ptr<ConversionPeriod>
          ( new ConversionPeriod(conversionStart, conversionEnd, 1) ) );

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);


  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
        
  pCallSchedule->AddCallPeriod
      ( shared_ptr<CallPeriod>
        ( CallPeriod::CreateWithStrike(startCallDate, endCallDate, 1) ) );
        
  pConvertibleBond->SetCallSchedule(pCallSchedule);

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

shared_ptr<finance::ConvertibleBond> 
InitCrossCurrencyCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
    conversionStart = Date(2003, Date::Feb, 1),
    conversionEnd = Date(2004, Date::May, 1);

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
    pInterests( new CashFlowStreamGeneral(issueDate, paymentDates, paymentRates,
                Date::DayCountConvention_Act365, Frequency_Annual) );
  
  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  conv->AddConversionPeriod
        ( shared_ptr<ConversionPeriod>
          ( new ConversionPeriod(conversionStart, conversionEnd, 1) ) );

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

  shared_ptr<Numeraire> pCurrencyUSD(new Numeraire("USD"));
  pConvertibleBond->SetNumeraire(pCurrencyUSD);

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);


  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
        
  pCallSchedule->AddCallPeriod
      ( shared_ptr<CallPeriod>
        ( CallPeriod::CreateWithStrike(startCallDate, endCallDate, 1) ) );
        
  pConvertibleBond->SetCallSchedule(pCallSchedule);

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01");

  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>
                           ( new YieldCurveFlat(0.0) ) );
 
  // Setup the rate data, and attach to the session data
  shared_ptr<Numeraire> pCurrencyUSD(new Numeraire("USD"));
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04));
  shared_ptr<YieldCurve> pycUSD(new YieldCurveFlat(0.02));
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);
  pRateData->SetYieldCurve(pCurrencyUSD, pycUSD);
  
  double dSpotFXRate = 0.8;
  shared_ptr<SpotFXRates> pSpotFXRates(new SpotFXRates());
  pSpotFXRates->SetFXRate(pCurrency, pCurrencyUSD, dSpotFXRate);
  pRateData->SetSpotFXRates(pSpotFXRates);

  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

shared_ptr<SessionData> InitFalseCrossCurrencySessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01");

  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>
                           ( new YieldCurveFlat(0.0) ) );
 
  // Setup the rate data, and attach to the session data
  shared_ptr<Numeraire> pCurrencyUSD(new Numeraire("USD"));
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04));
  shared_ptr<YieldCurve> pycUSD(new YieldCurveFlat(0.04));
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);
  pRateData->SetYieldCurve(pCurrencyUSD, pycUSD);
  
  double dSpotFXRate = 1.;
  shared_ptr<SpotFXRates> pSpotFXRates(new SpotFXRates());
  pSpotFXRates->SetFXRate(pCurrency, pCurrencyUSD, dSpotFXRate);
  pRateData->SetSpotFXRates(pSpotFXRates);

  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

void TestingHedge(const shared_ptr<ConvertibleBond>& pCB,
                  const shared_ptr<SessionData>& pSessionData,
                  const shared_ptr<ihg::TheoreticalModel>& pModel)
{
  std::cout.precision(15);
  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);

  std::cout << "The cb price is: " << output->GetPrice() << std::endl;
  std::cout << "The cb Delta is: " << output->GetDelta() << std::endl;
  
  if ( output->HasFXDelta() )
    std::cout << "The cb FXDelta is: " << output->GetFXDelta() << std::endl;
  
  std::cout << "The cb DefV  is: " << output->GetValueAfterDefault() 
            << std::endl;

  // hedge --------------------------------------------------------------
  Date IssueDate(2003, ito33::Date::Dec, 15);
  Date FirstDate(2004, ito33::Date::Jan, 1);
  
  Date MaturityDate(2010, ito33::Date::Jan, 1);

  shared_ptr<CashFlowStreamUniform> 
    pSpreadStream( new CashFlowStreamUniform
                       (
                         IssueDate,
                         FirstDate, 
                         MaturityDate,
                         0.02,
                         Date::DayCountConvention_Act365,
                         Frequency_SemiAnnual
                       )
                  );

  shared_ptr<CDS> pCDS(new CDS(0.4, pSpreadStream) );
  pCDS->SetSessionData(pSessionData);

  output = pModel->Compute(*pCDS);

  std::cout << "The cds price is: " << output->GetPrice() << std::endl;
  std::cout << "The cds Delta is: " << output->GetDelta() << std::endl;
  
  if ( output->HasFXDelta() )
    std::cout << "The cds FXDelta is: " << output->GetFXDelta() << std::endl;
  
  std::cout << "The cds DefV  is: " << output->GetValueAfterDefault() << std::endl;

  // option hedge --------------------------------------------------
  shared_ptr<Option> pOption(new Option(80, 
                                    Date(2004, Date::Jun, 1),
                                    Option_Call,
                                    ExerciseType_American
                                    )
                        );
  pOption->SetSessionData(pSessionData);
  output = pModel->Compute(*pOption);

  std::cout << "The option price is: " << output->GetPrice() << std::endl;
  std::cout << "The option Delta is: " << output->GetDelta() << std::endl;
  
  if ( output->HasFXDelta() )
    std::cout << "The option FXDelta is: " << output->GetFXDelta() << std::endl;
  
  std::cout << "The option DefV  is: " << output->GetValueAfterDefault() << std::endl;


  ///////////////////////////////////////////////////////////////////////////////
  
  PerfectHedgeRatios
    ratios = pModel->ComputePerfectHedgeRatios(*pCB, *pCDS);

  std::cout << "cds hedge " << ratios.GetDefaultHedgeRatio() << " "
            << ratios.GetUnderlyingHedgeRatio()  << " " 
            << ratios.GetFXHedgeRatio() << std::endl;


  ///////////////////////////////////////////////////////////////////////////////
  
  /*double dSpotTmp = pOption->GetSessionData()->GetEquity()->GetSpotSharePrice();
  pSessionData = InitSessionData();
  //pSessionData->GetEquity()->SetSpotSharePrice(dSpotTmp * 1.0001);
  pOption->SetSessionData(pSessionData);*/

  ratios = pModel->ComputePerfectHedgeRatios(*pCB, *pOption);

  std::cout << "option hedge " << ratios.GetDefaultHedgeRatio() << " "
            << ratios.GetUnderlyingHedgeRatio() << " " 
            << ratios.GetFXHedgeRatio() << std::endl;

  ///////////////////////////////////////////////////////////////////////////////
  {
    ratios = pModel->ComputePerfectHedgeRatios(*pCB, *pCB);

    std::cout << "self hedge " << ratios.GetDefaultHedgeRatio() << " "
              << ratios.GetUnderlyingHedgeRatio() << " " 
              << ratios.GetFXHedgeRatio() << std::endl;
  }

}

//*************************************************************************

int main()
{
  try {
    /*
    //--------------------------------------------------------------------------
    // Acceptance tests for call notice
    //--------------------------------------------------------------------------
    /* 
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(CallNoticeTest::suite());

    size_t result = runner.run("");
    */

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;  

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  pModel->SetDebugOutputFile("ihgcb.xml");

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>( new VolatilityFlat(dVol) ) );

  /*
  size_t nNbHRTimes = 2;
  static const Date pTimes[2] = { 38000, 40000};
  static const double pdValues[2] = { 0.02, .02 };

  pModel->SetHazardRate
          ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );
  */
   
  double dLambda = .2; 
  pModel->SetHazardRate( shared_ptr<ihg::HazardRate>
                         ( new ihg::HazardRateFlat(dLambda) ) );
  
  std::cout << std::endl << "Tests for a CB " << std::endl << std::endl;
   
  shared_ptr<SessionData> pSessionData = InitSessionData();
  
  shared_ptr<finance::ConvertibleBond> pCB = InitCB(pSessionData);

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01"));

  flags->SetComputeSurface(false);

  pCB->SetComputationalFlags(flags);

  TestingHedge(pCB, pSessionData, pModel);
  
  std::cout << std::endl << "Tests for a false cross-currency CB " << std::endl 
            << std::endl;
  
  shared_ptr<SessionData> pCC_SessionData = InitFalseCrossCurrencySessionData();
    shared_ptr<finance::ConvertibleBond> pCC_CB = InitCB(pCC_SessionData);  
  
  TestingHedge(pCC_CB, pCC_SessionData, pModel);
  
  std::cout << std::endl << "Tests for a cross-currency CB " << std::endl 
            << std::endl;
  
  pCC_SessionData = InitSessionData();
  pCC_CB = InitCrossCurrencyCB(pCC_SessionData);  
  
  TestingHedge(pCC_CB, pCC_SessionData, pModel);

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

