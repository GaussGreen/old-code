
#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/optionvarianceswap.h"
#include "ito33/finance/optiontype.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceVarianceSwap);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

void RunPricing();

int main()
{
  try
  {   

    // Run a pricing example
    RunPricing();

    return 0;
  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n" << e.GetFullMessage() << std::endl;

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

  pEquity->SetPreviousSharePrice( .96 * dS0 );

  if ( bHasDividend )
  {
    shared_ptr<Dividends> pDividends(new Dividends());
    pDividends->AddCash(valuationDate, .0000000000001);
    pEquity->SetDividends(pDividends);
  }

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    
  pEquity->SetBorrowCurve(pyf);
 
  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}


shared_ptr<ihg::TheoreticalModel> InitModel()
{
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.2;
  shared_ptr<Volatility> pVolatility(new VolatilityFlat(dVol));
  pModel->SetVolatility( pVolatility );

  //double dAlpha = 0.2;
  //double dBeta  = 0.5;
  //shared_ptr<Volatility> pVolatility(new VolatilityPower(dAlpha, dBeta, 100.0));
  //pModel->SetVolatility( pVolatility );

  double dLambda = 0.0;  
  shared_ptr<HazardRate> pHazardRate(new ihg::HazardRateFlat(dLambda));
  pModel->SetHazardRate( pHazardRate );

  return pModel;
}


void RunPricing()
{

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;

  bool bHasDividend = true;

  // Create the session data
  shared_ptr<SessionData> pSessionData = InitSessionData(bHasDividend);
  Date valuationDate = pSessionData->GetValuationDate();

  // Create the model
  shared_ptr<ihg::TheoreticalModel> pModel = InitModel();

  pModel->SetDebugOutputFile("./ihgvarswap.xml");

  // Create the variance swap
  Date maturityDate = valuationDate;
  maturityDate.AddMonths(6);
  Date startOfSamplingPeriod = valuationDate;
  //valuationDate.AddDays(-1);
  pSessionData->SetValuationDate( valuationDate );
  //startOfSamplingPeriod.AddDays(-2);
  //double dVolatilityStrike = 0.19961;
  double dVolatilityStrike = 0.2;
  SwapType swapType = Swap_Volatility;
  size_t nNbReturns = 252;

  shared_ptr<VarianceSwapTerms> 
    pTerms( new VarianceSwapTerms(maturityDate, swapType, 
                                  startOfSamplingPeriod, nNbReturns) );
  
  ReturnType returnType = Return_Actual;
  pTerms->SetReturnType(returnType);

  //pTerms->SetCapMultiplier( 1.2 );
  
  //VarianceSwap varianceSwap(pTerms, dVolatilityStrike);
  OptionVarianceSwap varianceSwap(pTerms, dVolatilityStrike, Option_Call);

  //pVarianceSwap->SetCurrentValues(dVolatilityStrike, 2);

  //pVarianceSwap->SetCapMultiplier( 5. );

  varianceSwap.SetSessionData(pSessionData);
   
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);

  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( valuationDate );

  Date analysisDate = valuationDate;   
  flags->SetAnalysisDate(analysisDate);

  varianceSwap.SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(varianceSwap);

  std::cout.precision(10);
  std::cout << "The ihg swap price is: " << output->GetPrice() << std::endl;
  
  if ( output->HasRho() )
    std::cout << "The swap Rho is: " << output->GetRho() << std::endl;

  if ( output->HasVega() )
    std::cout << "The swap Vega is: " << output->GetVega() << std::endl;

  bool bOutput = true;
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

  //ShowConvergence(*pModel, *pVarianceSwap, 4);

  std::cout << endl;
}
