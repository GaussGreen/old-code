#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/gammavarianceswap.h"
#include "ito33/finance/varianceswaption.h"
#include "ito33/finance/swaptype.h"
#include "ito33/finance/returntype.h"

#include "ito33/hg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceVarianceSwap);
ITO33_FORCE_LINK_MODULE(HGPriceVarianceSwaption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

using namespace std;

void RunVarianceSwapPricing();
void RunGammaVarianceSwapPricing();
void RunVarianceSwaptionPricing();

int main()
{
  try
  {   
    // Run a variance swap pricing example
    RunVarianceSwapPricing();

    // Run a gamma variance swaption pricing example
    // RunGammaVarianceSwapPricing();

    // Run a variance swaption pricing example
    // RunVarianceSwaptionPricing();

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

shared_ptr<SessionData> InitSessionData(bool bHasDividend)
{
  Date valuationDate(2002, Date::Jan, 1);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );  

  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  double dS0 = 100.0;
  shared_ptr<Equity> pEquity(new Equity(dS0, pCurrency));

  pEquity->SetPreviousSharePrice( dS0 );

  if ( bHasDividend )
  {
    shared_ptr<Dividends> pDividends(new Dividends());
    pDividends->AddCash(Date(2002, Date::Jul, 1), 2.0);
    pEquity->SetDividends(pDividends);
  }

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}


shared_ptr<hg::TheoreticalModel> InitModel()
{
  size_t nNbRegimes = 3;

  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.2);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.1);

  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess( new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  hg::Jumps jumps;

    
  jumps.push_back(hg::Jump(0.1, -0.2));
  pUnderlyingProcess->SetJumps(0, 0, jumps);
    
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

  pUnderlyingProcess->SetPostDefaultVolatility(0.10);

  shared_ptr<hg::TheoreticalModel> 
    pModel(new hg::TheoreticalModel(pUnderlyingProcess));

  return pModel;
}


void RunVarianceSwapPricing()
{

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;

  bool bHasDividend = false;

  // Create the session data
  shared_ptr<SessionData> pSessionData = InitSessionData(bHasDividend);
  Date valuationDate = pSessionData->GetValuationDate();

  // Create the model
  shared_ptr<hg::TheoreticalModel> pModel = InitModel();

  pModel->SetDebugOutputFile("hgvarswap.xml");

  // Create the variance swap
  Date maturityDate = valuationDate;
  maturityDate.AddMonths(6);
  Date startOfSamplingPeriod = valuationDate;
  double dVolatilityStrike = 0.19961;
  //double dVolatilityStrike = 0.001;
  SwapType swapType = Swap_Variance;
  ReturnType returnType = Return_Log;
  // ReturnType returnType = Return_Actual;
  size_t nNbReturns = 126;

  shared_ptr<VarianceSwapTerms> 
    pTerms( new VarianceSwapTerms(maturityDate, swapType, 
                                  startOfSamplingPeriod, nNbReturns) );
  
  pTerms->SetReturnType(returnType);

  VarianceSwap varianceSwap(pTerms, dVolatilityStrike);

  varianceSwap.SetSessionData(pSessionData);
    
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);

  Date analysisDate = valuationDate;  
  
  flags->SetAnalysisDate(analysisDate);

  //flags->SetComputeSurface(true);

  varianceSwap.SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(varianceSwap);

  std::cout.precision(10);
  std::cout << "The hg swap price is: " << output->GetPrice() << std::endl;

  if ( output->HasRho() )
    std::cout << "The swap Rho is: " << output->GetRho() << std::endl;

  if ( output->HasVega() )
    std::cout << "The swap Vega is: " << output->GetVega() << std::endl;

  
  bool bOutput = false;
  if ( bOutput && output->HasPriceAtAnalysisDate() )
  {
    std::cout.precision(12);
    std::cout << "Analysis date data" << std::endl;
    Values pdSpots = output->GetSpotsAtAnalysisDate();
    Values pdPrices = output->GetPricesAtAnalysisDate();
    Values pdDeltas = output->GetDeltasAtAnalysisDate();
    Values pdGammas = output->GetGammasAtAnalysisDate();

    for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
      std::cout << pdSpots[nIdx] 
                << " " << pdPrices[nIdx] 
                << " " << pdDeltas[nIdx] 
                << " " << pdGammas[nIdx] 
                << std::endl;
    std::cout << std::endl;
  }
 
  bOutput = false;
  if ( bOutput && output->HasPriceSurface() )
  {
    std::cout << "Surface at valuation date" << std::endl;
    Domain::Spots pdSpots(20);
    double dSpot = pSessionData->GetSpotSharePrice(); 
    pdSpots[0] =  dSpot * 0.8;
    double dIncrement = ( dSpot - pdSpots[0] ) / 10.0;
    for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
      pdSpots[nIdxS] = pdSpots[nIdxS - 1] + dIncrement;

    SharedSurface priceSurface = output->GetPriceSurface();
    SharedSurface deltaSurface = output->GetDeltaSurface();
    SharedSurface gammaSurface = output->GetGammaSurface();

    priceSurface->GetDomain()->SetUnderlyingSharePrices(pdSpots);
    Domain::Dates pdDates = priceSurface->GetDomain()->GetDates();

    for (size_t nIdxDate = 0; nIdxDate < pdDates.size(); nIdxDate++)
    {
      //std::cout << "date = " << pdDates[nIdxDate] << std::endl;
      SurfaceDouble::Doubles pdPrices = priceSurface->GetValuesAt(nIdxDate);

      SurfaceDouble::Doubles pdDeltas = deltaSurface->GetValuesAt(nIdxDate);

      SurfaceDouble::Doubles pdGammas = gammaSurface->GetValuesAt(nIdxDate);

      for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
        std::cout << nIdxDate << " " << pdSpots[nIdx] 
                  << " " << pdPrices[nIdx] 
               //   << " " << pdDeltas[nIdx] 
               //   << " " << pdGammas[nIdx] 
                  << std::endl;
       
      std::cout << std::endl;
    }
  }

  ShowConvergence(*pModel, varianceSwap, 4);

  std::cout << endl;
}

void RunGammaVarianceSwapPricing()
{

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;

  bool bHasDividend = false;

  // Create the session data
  shared_ptr<SessionData> pSessionData = InitSessionData(bHasDividend);
  Date valuationDate = pSessionData->GetValuationDate();

  // Create the model
  shared_ptr<hg::TheoreticalModel> pModel = InitModel();

  // pModel->SetDebugOutputFile("hgvarswap.xml");

  // Create the variance swap
  Date maturityDate = valuationDate;
  maturityDate.AddMonths(6);
  Date startOfSamplingPeriod = valuationDate;
  double dVolatilityStrike = 0.19961;
  //double dVolatilityStrike = 0.001;
  SwapType swapType = Swap_Variance;
  ReturnType returnType = Return_Log;
  // ReturnType returnType = Return_Actual;
  size_t nNbReturns = 126;

  shared_ptr<VarianceSwapTerms> 
    pTerms( new VarianceSwapTerms(maturityDate, swapType, 
                                  startOfSamplingPeriod, nNbReturns) );
  
  pTerms->SetReturnType(returnType);

  GammaVarianceSwap varianceSwap(pTerms, dVolatilityStrike);

  varianceSwap.SetSessionData(pSessionData);
    
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);

  Date analysisDate = valuationDate;  
  
  // flags->SetAnalysisDate(analysisDate);

  // flags->SetComputeSurface(true);

  varianceSwap.SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(varianceSwap);

  std::cout.precision(10);
  std::cout << "The hg swap price is: " << output->GetPrice() << std::endl;

  if ( output->HasRho() )
    std::cout << "The swap Rho is: " << output->GetRho() << std::endl;

  if ( output->HasVega() )
    std::cout << "The swap Vega is: " << output->GetVega() << std::endl;

  ShowConvergence(*pModel, varianceSwap, 4);

  std::cout << endl;
}

void RunVarianceSwaptionPricing()
{

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;

  bool bHasDividend = false;

  // Create the session data
  shared_ptr<SessionData> pSessionData = InitSessionData(bHasDividend);
  Date valuationDate = pSessionData->GetValuationDate();

  // Create the model
  shared_ptr<hg::TheoreticalModel> pModel = InitModel();

  pModel->SetDebugOutputFile("hg_varswaption.xml");

  // Create the variance swap terms
  Date maturityDate = valuationDate;
  maturityDate.AddMonths(6);

  Date vsMaturityDate = maturityDate;
  vsMaturityDate.AddMonths(6);

  Date startOfSamplingPeriod = maturityDate;
  SwapType swapType = Swap_Variance;
  ReturnType returnType = Return_Log;
  // ReturnType returnType = Return_Actual;
  size_t nNbReturns = 126;

  shared_ptr<VarianceSwapTerms> 
    pTerms( new VarianceSwapTerms(vsMaturityDate, swapType, 
                                  startOfSamplingPeriod, nNbReturns) );
  
  pTerms->SetReturnType(returnType);

  double dStrike = 0.2;
  VarianceSwaption varianceSwaption(pTerms, Option_Call, dStrike, maturityDate);

  varianceSwaption.SetSessionData(pSessionData);
    
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);

  Date analysisDate = valuationDate;  
  
  flags->SetAnalysisDate(analysisDate);

  //flags->SetComputeSurface(true);

  varianceSwaption.SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(varianceSwaption);

  std::cout.precision(10);
  std::cout << "The hg swaption price is: " << output->GetPrice() << std::endl;
  
  if ( output->HasRho() )
    std::cout << "The swap Rho is: " << output->GetRho() << std::endl;

  if ( output->HasVega() )
    std::cout << "The swap Vega is: " << output->GetVega() << std::endl;

  ShowConvergence(*pModel, varianceSwaption, 4);

  std::cout << endl;
}
