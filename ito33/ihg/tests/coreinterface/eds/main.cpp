#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h" 
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/domain.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratepower.h"

#ifdef ITO33_TEST_MODENV
#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#endif

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceEDS);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

void RunPricing();
void RunMeshTests();
extern void RunImpliedSpreads();


int main()
{
  try
  {   

    // Run a pricing example
    RunPricing();

    // Run an implied spread calibration example
//    RunImpliedSpreads();

    // Run some mesh tests
//    RunMeshTests();

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


void RunPricing()
{

  bool
    bComputeRho = true,
    bComputeVega = true,
    bComputeSurface = false;

  Date valuationDate(2002, Date::Jan, 1);

  double dS0 = 100.0;
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dS0, pCurrency));

  shared_ptr<Dividends> pDividends(new Dividends());
  pDividends->AddCash(Date(2002, Date::Jul, 1), 2.0);
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    
  pEquity->SetBorrowCurve(pyf);
    
  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  //dAnnualRate = dContinuousRate;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  pModel->SetDebugOutputFile("ihgcds.xml");

  double dVol = 0.2;
  //double dVol = 0.2369;
  shared_ptr<Volatility> pVolatility(new VolatilityFlat(dVol));
  pModel->SetVolatility( pVolatility );

  double dLambda = .02;
  double dBeta = 1.2;
  //dLambda = 0.0;  
  //shared_ptr<HazardRate> pHazardRate(new ihg::HazardRateFlat(dLambda));
  shared_ptr<HazardRate> pHazardRate(
    new ihg::HazardRatePower(dLambda, dBeta, dS0));
  pModel->SetHazardRate( pHazardRate );



  // the actual recovery is 1.0 - recoveryrate
  double dRecoveryRate = 0.0;

  // multiply by 2 since payments are semi-annual
  double dSpread  = 0.009412 * 2.0;
  Date firstSpreadDate(2002, Date::Jul, 1);
  Date lastSpreadDate(2007, Date::Jan, 1);

  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( new CashFlowStreamUniform
                      (
                        valuationDate,
                        firstSpreadDate,
                        lastSpreadDate,
                        dSpread,
                        Date::DayCountConvention_Act365,
                        Frequency_SemiAnnual
                      )
                );

  std::cout << "Maturity = " << pSpreadStream->GetLastPaymentDate() << std::endl;
  std::cout << "Spread = " << dSpread << std::endl;
  std::cout << "Recovery = " << dRecoveryRate << std::endl;

  /*
  std::cout << "The coupon payments:" << std::endl;
  CashFlowStream::const_iterator iter;
  for (iter = pSpreadStream->begin();
      iter != pSpreadStream->end();
      ++iter)
  {
    std::cout << "Date = " << iter->first
              << ", amount = " << iter->second 
              << std::endl;
  }
  */

  std::cout << endl;
  
  double dBarrier = 30.0;
  shared_ptr<finance::EDS>
    pEDS( new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);

  Date analysisDate = valuationDate;
  analysisDate.AddMonths(2);

  //Date tmp(2007, Date::Jan, 1);
  //analysisDate = tmp;
  //analysisDate.AddDays(-1);

  flags->SetAnalysisDate(analysisDate);

  flags->SetComputeSurface(true);

  pEDS->SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pEDS);

  std::cout.precision(10);
  std::cout << "The EDS price is: " << output->GetPrice() << std::endl;
  
  if ( output->HasRho() )
    std::cout << "The EDS Rho is: " << output->GetRho() << std::endl;

  if ( output->HasVega() )
    std::cout << "The EDS Vega is: " << output->GetVega() << std::endl;

  bool bTestExtraPoint = true;
  if ( bTestExtraPoint && output->HasPriceAtAnalysisDate() )
  {
    std::cout.precision(12);
    std::cout << "Analysis date data" << std::endl;
    finance::Values pdSpots = output->GetSpotsAtAnalysisDate();
    finance::Values pdPrices = output->GetPricesAtAnalysisDate();
    finance::Values pdDeltas = output->GetDeltasAtAnalysisDate();
    finance::Values pdGammas = output->GetGammasAtAnalysisDate();

    for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
      std::cout << pdSpots[nIdx] 
                << " " << pdPrices[nIdx] 
                << " " << pdDeltas[nIdx] 
                << " " << pdGammas[nIdx] 
                << std::endl;
    std::cout << std::endl;
  }
   

  if ( bTestExtraPoint && output->HasPriceSurface() )
  {
    std::cout << "Surface at valuation date" << std::endl;
    finance::Domain::Spots pdSpots(20);
    pdSpots[0] = dBarrier - 1.0;
    for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
      pdSpots[nIdxS] = pdSpots[nIdxS - 1] + 0.5;

    finance::SharedSurface priceSurface = output->GetPriceSurface();
    finance::SharedSurface deltaSurface = output->GetDeltaSurface();
    finance::SharedSurface gammaSurface = output->GetGammaSurface();

    priceSurface->GetDomain()->SetUnderlyingSharePrices(pdSpots);
    finance::Domain::Dates pdDates = priceSurface->GetDomain()->GetDates();

    for (size_t nIdxDate = 0; nIdxDate < pdDates.size(); nIdxDate++)
    {
      //std::cout << "date = " << pdDates[nIdxDate] << std::endl;
      finance::SurfaceDouble::Doubles pdPrices 
        = priceSurface->GetValuesAt(nIdxDate);

      finance::SurfaceDouble::Doubles pdDeltas 
        = deltaSurface->GetValuesAt(nIdxDate);

      finance::SurfaceDouble::Doubles pdGammas 
        = gammaSurface->GetValuesAt(nIdxDate);

      for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
        std::cout << nIdxDate << " " << pdSpots[nIdx] 
               //   << " " << pdPrices[nIdx] 
                  << " " << pdDeltas[nIdx] 
               //   << " " << pdGammas[nIdx] 
                  << std::endl;
       
      std::cout << std::endl;
    }
  }

  //ShowConvergence(*pModel, *pEDS, 5);

  std::cout << endl;
}


void RunMeshTests()
{

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;

  Date valuationDate(2002, Date::Jan, 1);

  double dS0 = 100.0;
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dS0, pCurrency));

  shared_ptr<Dividends> pDividends(new Dividends());
  pDividends->AddCash(Date(2002, Date::Jul, 1), 2.0);
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    
  pEquity->SetBorrowCurve(pyf);
    
  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  //dAnnualRate = dContinuousRate;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));  
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.2;
  //double dVol = 0.2369;
  shared_ptr<Volatility> pVolatility(new VolatilityFlat(dVol));
  pModel->SetVolatility( pVolatility );

  double dLambda = .02;
  //double dBeta = 1.2;
  //dLambda = 0.0;  
  shared_ptr<HazardRate> pHazardRate(new ihg::HazardRateFlat(dLambda));
  //shared_ptr<HazardRate> pHazardRate(
  //  new ihg::HazardRatePower(dLambda, dBeta, dS0));
  pModel->SetHazardRate( pHazardRate );

  // the actual recovery is 1.0 - recoveryrate
  double dRecoveryRate = 0.0;

  // multiply by 2 since payments are semi-annual
  //double dSpread  = 0.009412 * 2.0;
  double dSpread  = 0.02;
  Date firstSpreadDate(2002, Date::Jul, 1);
  Date lastSpreadDate(2007, Date::Jan, 1);

  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( new CashFlowStreamUniform
                      (
                        valuationDate,
                        firstSpreadDate,
                        lastSpreadDate,
                        dSpread,
                        Date::DayCountConvention_Act365,
                        Frequency_SemiAnnual
                      )
                );

  std::cout << "Maturity = " << pSpreadStream->GetLastPaymentDate() << std::endl;
  std::cout << "Spread = " << dSpread << std::endl;
  std::cout << "Recovery = " << dRecoveryRate << std::endl;

  /*
  std::cout << "The coupon payments:" << std::endl;
  CashFlowStream::const_iterator iter;
  for (iter = pSpreadStream->begin();
      iter != pSpreadStream->end();
      ++iter)
  {
    std::cout << "Date = " << iter->first
              << ", amount = " << iter->second 
              << std::endl;
  }
  */

  std::cout << endl;
  /*
  double dBarrier = 3.0;  
  shared_ptr<finance::EDS>
    pEDS( new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);

  pEquity->SetSpotSharePrice(4.0);
  shared_ptr<ihg::ModelOutput> output = pModel->Compute(*pEDS);
  std::cout.precision(10);
  std::cout << "The EDS price is: " << output->GetPrice() << std::endl;
  */

  
  double dBarrier = 30.0;
  shared_ptr<finance::EDS>
    pEDS( new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);
  
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);

  pEDS->SetComputationalFlags(flags);

  std::cout << "Barrier at " << dBarrier << std::endl << std::endl;
  std::cout << " spot = " << 29.0 << std::endl;
  pEquity->SetSpotSharePrice(29.0);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;

  std::cout << " spot = " << 31.0 << std::endl;
  pEquity->SetSpotSharePrice(31.0);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;

  std::cout << " spot = " << 100.0 << std::endl;
  pEquity->SetSpotSharePrice(100.0);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;


  dBarrier = 3.0;
  pEDS = make_ptr( new EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);
  std::cout << "Barrier at " << dBarrier << std::endl << std::endl;
  std::cout << " spot = " << 4.0 << std::endl;
  pEquity->SetSpotSharePrice(4.0);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;

  std::cout << " spot = " << 31.0 << std::endl;
  pEquity->SetSpotSharePrice(31.0);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;

  std::cout << " spot = " << 100.0 << std::endl;
  pEquity->SetSpotSharePrice(100.0);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;


  dBarrier = 300000;
  pEDS = make_ptr( new EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);
  std::cout << "Barrier at " << dBarrier << std::endl << std::endl;
  std::cout << " spot = " << 310000 << std::endl;
  pEquity->SetSpotSharePrice(310000);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;

  std::cout << " spot = " << 1000000 << std::endl;
  pEquity->SetSpotSharePrice(1000000);
  ShowConvergence(*pModel, *pEDS, 5);
  std::cout << std::endl;

}

