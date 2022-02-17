#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"
#include "ito33/dateutils.h"

#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/frequency.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/exoticoption/curran.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/blackscholes.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/tests/showconvergence.h"

ITO33_FORCE_LINK_MODULE(IHGPriceAsianOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

shared_ptr<SessionData> InitSessionData(Date valuationDate);
void BasicTest();
void ChangeValuationDate();

shared_ptr<SessionData> InitSessionData(Date valuationDate)
{ 

  double dSpot = 100.0;
  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;

  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));

  shared_ptr<Equity> pEquity( new Equity(dSpot, pCurrency) );

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
  
  pEquity->SetBorrowCurve(pyf);

  //if ( bFakeDividend )
  //{
  //  shared_ptr<Dividends> pDiv = new Dividends();
  //  pDiv->Add(Dividend::Type::Cash, Date(1900,Date::Jun,1),1.);
  //  pEquity->SetDividends(pDiv);
  //}

  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  shared_ptr<RateData> pRateData(new RateData);

  pRateData->SetYieldCurve(pCurrency, pyc);
  
  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}


void ChangeValuationDate()
{

  std::cout << "Change the valuation date" << std::endl;

  // Setup the session
  Date valuationDate(2005, Date::Jan, 1);

  // Set the asian option data
  double dStrike = 100.0;
  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);
  OptionType optionType = Option_Put;
  ExerciseType exerType = ExerciseType_European;
  Date issueDate = valuationDate;

  // setup the model
  double dVolatility = .2;
  //double dHazardRate = 0.05;
  double dHazardRate = 0.0;
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
  pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(false);
  flags->SetComputeVega(false);
  flags->SetComputeSurface(false);

  // loop over valuation dates
  Date endDate = maturityDate;
  endDate.AddDays(-1);

  std::cout << "Issue date: " << issueDate << std::endl;
  std::cout << "Maturity date: " << maturityDate << std::endl;
  std::cout << std::endl;
  
  while (valuationDate < endDate)
  {
    // Setup the session data
    shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

    // Setup the asian option
    shared_ptr<AsianOption> asianOpt(new 
    AsianOption(dStrike, maturityDate, optionType, exerType, issueDate, 252));

    double dSpot = pSessionData->GetSpotSharePrice();
    asianOpt->SetCurrentAverage(dSpot, 10);

    asianOpt->SetSessionData(pSessionData);

    asianOpt->SetComputationalFlags(flags);

    // price
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);

    std::cout << "valuation date = " << valuationDate 
              << ", price = " << output->GetPrice() << std::endl;

    // Move to next valuation date
    valuationDate.AddMonths(1);
  }

  std::cout << std::endl << std::endl;
}

void ChangeCurrentAverageAndValuation(size_t nWhich)
{
  std::cout << "Change the current average and valuation date" << std::endl;

  // Set the session data
  Date valuationDate(2005, Date::Jan, 1);

  // Set the asian option data
  double dStrike = 0.0;
  OptionType optionType;
  if (nWhich == 0)
  {
    dStrike = 100.0;
    optionType = Option_Call;
    std::cout << "Type: fixed call, strike 100" << std::endl;
  }
  else if (nWhich == 1)
  {
    dStrike = 100.0;
    optionType = Option_Put;
    std::cout << "Type: fixed put, strike 100" << std::endl;
  }
  else if (nWhich == 2)
  {
    optionType = Option_Call;
    std::cout << "Type: floating call" << std::endl;
  }
  else
  {
    optionType = Option_Put;
    std::cout << "Type: floating put" << std::endl;
  }

  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);  
  ExerciseType exerType = ExerciseType_European;
  Date issueDate = valuationDate;

  // setup the model
  double dVolatility = .2;
  //double dHazardRate = 0.05;
  double dHazardRate = 0.0;
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
  pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(false);
  flags->SetComputeVega(false);
  flags->SetComputeSurface(false);

  // loop over valuation dates
  Date endDate = maturityDate;
  endDate.AddDays(-1);

  std::cout << "Issue date: " << issueDate << std::endl;
  std::cout << "Maturity date: " << maturityDate << std::endl;
  std::cout << std::endl;
  
  while (valuationDate < endDate)
  {
    // Setup the session data
    shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

    // loop over current averages
    double dCurrentAverage = 50.0;

    for (dCurrentAverage = 50.0; dCurrentAverage < 151.0; dCurrentAverage += 25.0)
    {   
      // Setup the asian option
      shared_ptr<AsianOption> asianOpt;
      if (dStrike > 0.0)
        asianOpt = make_ptr( new AsianOption(dStrike, maturityDate, optionType, 
                                             exerType, issueDate, 252) );
      else
        asianOpt = make_ptr( new AsianOption(maturityDate, optionType, 
                                             exerType, issueDate, 252) );
      
      asianOpt->SetCurrentAverage(dCurrentAverage, 5);

      asianOpt->SetSessionData(pSessionData);

      asianOpt->SetComputationalFlags(flags);

      // price
      shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);

      std::cout << "valuation date = " << valuationDate
                << ", current average = " << dCurrentAverage
                << ", price = " << output->GetPrice() << std::endl;

    }
     
    // Move to next valuation date
    valuationDate.AddMonths(3);
    std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;
}

void BasicTestVSCurran()
{
  double dSpot           = 100.0;
  double dContinuousRate = 0.09; 
  
  double ONEWEEK       = 1./52.;
  double dMaturityTime = 72.*ONEWEEK;
  double dAvgStartTime = 20.*ONEWEEK;
  double dDividend     = 0.;

  double dVolatility = .1;
  double dStrike     = 100;


  size_t nNbAveraging = 53;

  Date valuationDate   = Date(2006, Date::Jan, 1);
  Date maturityDate    = GetDateFrom( dMaturityTime + GetDoubleFrom( valuationDate ) );
  Date avgStartDate    = GetDateFrom( dAvgStartTime + GetDoubleFrom( valuationDate )  - 1/52.);

  double dPriceCall = CurranCall(dSpot, dStrike, dMaturityTime, nNbAveraging,
                                 dAvgStartTime, dVolatility, dContinuousRate, dDividend);

  double dAnnualRate     = exp(dContinuousRate) - 1.0;

  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));

  shared_ptr<Equity> pEquity( new Equity(dSpot, pCurrency) );

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.));    
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<YieldCurve> pyc( new YieldCurveFlat( dAnnualRate ) );

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  shared_ptr<AsianOption> asianOpt( new AsianOption(dStrike, maturityDate, Option_Call, 
    ExerciseType_European, avgStartDate, nNbAveraging) );

  asianOpt->SetSessionData(pSessionData);
  

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );

  pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(0.)) );

  pModel->SetDebugOutputFile("./asian.xml");

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);

  std::cout.precision(10);
  std::cout << std::endl;
  std::cout << "Asian  (pde) price = " << output->GetPrice() << std::endl;
  std::cout << "Currant price      = " << dPriceCall << std::endl;

 // ShowConvergence(*pModel, *asianOpt, 4);
}

void BasicTest()
{
  Date valuationDate(2003, Date::Jan, 1);
  Date averageStartDate = valuationDate;
  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);

  double dSpot = 100.0;
  double dStrike = 100.0;
  double dContinuousRate = 0.09;
  double dAnnualRate = exp(dContinuousRate) - 1.0;

  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));

  shared_ptr<Equity> pEquity( new Equity(dSpot, pCurrency) );

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
  
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  OptionType optionType = Option_Call;
  ExerciseType exerType   = ExerciseType_European;
 
  shared_ptr<AsianOption> asianOpt(new 
    AsianOption(dStrike, maturityDate, optionType, exerType, averageStartDate, 53));

  asianOpt->SetSessionData(pSessionData);

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(false);
  flags->SetComputeVega(false);
  flags->SetComputeSurface(true);
  //Date testDate = maturityDate;
  //flags->SetAnalysisDate( testDate.AddDays(-1) );
  flags->SetAnalysisDate( pSessionData->GetValuationDate() );

  asianOpt->SetComputationalFlags(flags);

  // setup the model
  double dVolatility = .1;
  double dHazardRate = 0.0;
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVolatility)) );
  pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dHazardRate)) );

  pModel->SetDebugOutputFile("./asiant_ihg.xml");

  // price
  shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);

  std::cout.precision(10);
  std::cout << std::endl;
  std::cout << "Asian  (pde) price = " << output->GetPrice() << std::endl;
  std::cout << "             delta = " << output->GetDelta() << std::endl;
  std::cout << "             gamma = " << output->GetGamma() << std::endl;
  std::cout << "             theta = " << output->GetTheta() << std::endl;

  valuationDate.AddDays(1);
  averageStartDate = valuationDate;
  pSessionData = make_ptr( new SessionData(pRateData, pEquity, valuationDate) );
  shared_ptr<AsianOption> asianOpt2(new 
    AsianOption(dStrike, maturityDate, optionType, exerType, averageStartDate, 53));
  asianOpt2->SetSessionData(pSessionData);

  asianOpt2->SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output2 = pModel->Compute(*asianOpt2);
  double dPrice = output->GetPrice();
  double dPriceShifted = output2->GetPrice();

  double dThetaFD = (dPriceShifted - dPrice) / (1./365.);
  std::cout << "             theta fd = " << dThetaFD << std::endl;

  averageStartDate.AddDays(-1);
  pSessionData = make_ptr( new SessionData(pRateData, pEquity, valuationDate) );
  shared_ptr<AsianOption> asianOpt3(new 
    AsianOption(dStrike, maturityDate, optionType, exerType, averageStartDate, 252));
  asianOpt3->SetCurrentAverage(dSpot, 10);
  asianOpt3->SetSessionData(pSessionData);

  shared_ptr<finance::ModelOutput> output3 = pModel->Compute(*asianOpt3);  
  dPriceShifted = output3->GetPrice();

  dThetaFD = (dPriceShifted - dPrice) / (1./365.);
  std::cout << "      theta fd avg = " << dThetaFD << std::endl;

  std::cout << std::endl;

  if ( output->HasPriceSurface() )
  {
    std::cout << "Surface computed" << std::endl;

    SharedSurface pSurface = output->GetPriceSurface();

    shared_ptr<Domain> pDomain = pSurface->GetDomain();

    Domain::Spots pdSpots(20);
    pdSpots[0] = 0.8 * pSessionData->GetSpotSharePrice();
    double dStep = 0.1 * (pSessionData->GetSpotSharePrice() - pdSpots[0]);
    for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
      pdSpots[nIdxS] = pdSpots[nIdxS-1] + dStep;

    pDomain->SetUnderlyingSharePrices(pdSpots);

    size_t nNbDates = pDomain->GetDates().size();
    SurfaceDouble::Doubles pValues = pSurface->GetValuesAt(nNbDates-1);

    for (size_t nIdx = 0; nIdx < pValues.size(); nIdx++)
      std::cout << pdSpots[nIdx] << ", value = " << pValues[nIdx]
                << std::endl;
    
    std::cout << std::endl;
  }

  if ( output->HasPriceAtAnalysisDate() )
  {
    std::cout << "Analysis date values computed:" << std::endl;
    std::vector<double> pdSpots = output->GetSpotsAtAnalysisDate();
    std::vector<double> pdPrices = output->GetPricesAtAnalysisDate();

    for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
    {
      std::cout << pdSpots[nIdx] << " "
                << pdPrices[nIdx] << " "
                << std::endl;
    }
    std::cout << std::endl;
  }

  ShowConvergence(*pModel, *asianOpt, 4);

}

int main()
{
  try
  {
    
    //BasicTest();

    BasicTestVSCurran();

    
    //ChangeValuationDate();

    //ChangeCurrentAverageAndValuation(0);
    //ChangeCurrentAverageAndValuation(1);
    //ChangeCurrentAverageAndValuation(2);
    //ChangeCurrentAverageAndValuation(3);
  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:" << std::endl
              << e.GetFullMessage() << std::endl;

    return 1;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught." << std::endl;

    return 2;
  }

  return 0;
}
