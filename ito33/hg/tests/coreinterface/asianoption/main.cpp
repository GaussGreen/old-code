#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/frequency.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/exoticoption/curran.h"
#include "ito33/finance/domain.h"

#include "ito33/hg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"

ITO33_FORCE_LINK_MODULE(HGPriceAsianOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;


shared_ptr<SessionData> InitSessionData(Date valuationDate);
void BasicTest();
void ChangeValuationDate();

shared_ptr<SessionData> InitSessionData(Date valuationDate)
{ 
  double dSpot = 100.0;
  
  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;

  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  shared_ptr<RateData> pRateData(new RateData);

  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<Equity> pEquity( new Equity(dSpot, pCurrency) );

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
  
  pEquity->SetBorrowCurve(pyf);

  //if ( bFakeDividend )
  //{
  //  shared_ptr<Dividends> pDiv = new Dividends();
  //  pDiv->Add(Dividend::Type::Cash, Date(1900,Date::Jun,1),1.);
  //  pEquity->SetDividends(pDiv);
  //}

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<hg::TheoreticalModel> InitModel()
{
  size_t nNbRegimes = 1;

  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.2);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.0);

  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess( new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  hg::Jumps jumps;

    
  //jumps.push_back(hg::Jump(0.1, -0.2));
  //pUnderlyingProcess->SetJumps(0, 0, jumps);
   
  if (nNbRegimes > 1)
  {
      
    jumps.clear();
    jumps.push_back(hg::Jump(0.15, -0.25));
    pUnderlyingProcess->SetJumps(0, 1, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.2, -0.5));
    pUnderlyingProcess->SetJumps(1, 1, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.11, 0.35));
    pUnderlyingProcess->SetJumps(1, 0, jumps); 
      
  }

  if (nNbRegimes > 2)
  {
    jumps.clear();
    jumps.push_back(hg::Jump(0.1, -0.2));
    pUnderlyingProcess->SetJumps(0, 2, jumps); 
 
    jumps.clear();
    jumps.push_back(hg::Jump(0.1, -0.3));
    pUnderlyingProcess->SetJumps(1, 2, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.3, -0.2));
    pUnderlyingProcess->SetJumps(2, 0, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.2, -0.4));
    pUnderlyingProcess->SetJumps(2, 1, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.05, -0.25));
    pUnderlyingProcess->SetJumps(2, 2, jumps); 
  }

  shared_ptr<hg::TheoreticalModel> 
    pModel(new hg::TheoreticalModel(pUnderlyingProcess));

  return pModel;
}


void ChangeValuationDate()
{

  std::cout << "Change the valuation date" << std::endl;

  // Setup the session
  Date avgStartDate(2005, Date::Jan, 1);

  // Set the asian option data
  double dStrike        = 100.0;
  Date maturityDate     = avgStartDate;
  maturityDate.AddYears(1);
  OptionType optionType = Option_Put;
  ExerciseType exerType = ExerciseType_European;
  Date valuationDate    = avgStartDate;

  size_t nNbSampling     = 365;

  // setup the model
  shared_ptr<hg::TheoreticalModel> pModel = InitModel();

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(false);
  flags->SetComputeVega(false);
  flags->SetComputeSurface(false);

  // loop over valuation dates
  Date endDate = maturityDate;
  endDate.AddDays(-1);

  std::cout << "Average Start date: " << avgStartDate << std::endl;
  std::cout << "Maturity date: " << maturityDate << std::endl;
  std::cout << std::endl;
  size_t nNbSamplesUsed = 0;

  while (valuationDate < endDate)
  {
    // Setup the session data
    shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

    // Setup the asian option
    shared_ptr<AsianOption> asianOpt(new 
    AsianOption(dStrike, maturityDate, optionType, exerType, avgStartDate, nNbSampling) );

    double dSpot = pSessionData->GetSpotSharePrice();

    asianOpt->SetCurrentAverage(dSpot,  nNbSamplesUsed);

    asianOpt->SetSessionData(pSessionData);

    asianOpt->SetComputationalFlags(flags);

    // price
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);

    std::cout << "valuation date = " << valuationDate 
              << ", price = " << output->GetPrice() << std::endl;

    // Move to next valuation date
    valuationDate.AddMonths(1);

    nNbSamplesUsed += 30 ;
  }

  std::cout << std::endl << std::endl;
}


void BasicTest()
{
  Date valuationDate(2003, Date::Jan, 1);
  Date averageStartDate = valuationDate;
  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);

  double dSpot = 100.0;
  double dStrike = 100.0;

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  double dContinuousRate = 0.09;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<Equity> pEquity( new Equity(dSpot, pCurrency) );

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
  
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  OptionType optionType   = Option_Call;
  ExerciseType exerType   = ExerciseType_European;
 
  shared_ptr<AsianOption> asianOpt(new 
    AsianOption(dStrike, maturityDate, optionType, exerType, averageStartDate, 52));

  asianOpt->SetSessionData(pSessionData);

  // setup the model
  shared_ptr<hg::TheoreticalModel> pModel = InitModel();

  //pModel->SetDebugOutputFile("./asian_hg.xml");
  
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  //flags->SetComputeRho(false);
  //flags->SetComputeVega(false);
  //flags->SetComputeSurface(true);
  //Date testDate = maturityDate;
  //flags->SetAnalysisDate( testDate.AddDays(-1) );
  flags->SetAnalysisDate( pSessionData->GetValuationDate() );
  
  asianOpt->SetComputationalFlags(flags);

  // price
  shared_ptr<finance::ModelOutput> output = pModel->Compute(*asianOpt);

  std::cout.precision(10);
  std::cout << std::endl;
  std::cout << "Asian  (pde) price = " << output->GetPrice() << std::endl;
  std::cout << "             delta = " << output->GetDelta() << std::endl;
  std::cout << "             gamma = " << output->GetGamma() << std::endl;
  std::cout << "             theta = " << output->GetTheta() << std::endl;

  

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

  //ShowConvergence(*pModel, *asianOpt, 4);
}

int main()
{
  try
  {
    // This was the original code
    BasicTest();

   
    //ChangeValuationDate();

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
