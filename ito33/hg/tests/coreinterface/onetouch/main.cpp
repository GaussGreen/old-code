#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/domain.h"

#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/jumps.h"

#include "hg/computesensitivity.h"
#include "hg/numoutput.h"

#ifdef ITO33_TEST_MODENV
#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#endif

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceOneTouch);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

shared_ptr<SessionData> InitSessionData()
{
  Date valuationDate(2002, Date::Jul, 1);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  double dSpot = 30.0;
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));

  shared_ptr<Dividends> pDividends( new Dividends() );
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.00) );
     
  pEquity->SetBorrowCurve(pyf);
     
  double dRate = 0.05;
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<hg::UnderlyingProcess> InitProcess()
{
  size_t nNbRegimes = 2;
  std::vector<double> pdVols(nNbRegimes);
  std::vector<double> pdDefaultIntensities(nNbRegimes);
  
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbRegimes; nIdx++)
  {
    pdVols[nIdx] = 0.2 + nIdx * 0.1;
    pdDefaultIntensities[nIdx] = 0.05;
  }

  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess(new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities));

  if (nNbRegimes > 1)
  {
    std::vector<double> pdIntensity(1);
    std::vector<double> pdAmplitude(1);

    size_t nIdx1, nIdx2;
    for (nIdx1 = 0; nIdx1 < nNbRegimes; nIdx1++)
    {
      for (nIdx2 = 0; nIdx2 < nNbRegimes; nIdx2++)
      {
        pdIntensity[0] = 0.1;
        pdAmplitude[0] = 0.0;

        if (nIdx1 != nIdx2)
          pUnderlyingProcess->SetJumps(nIdx1, nIdx2, pdIntensity, pdAmplitude);
      }
    }
  } // if more than one regime

  return pUnderlyingProcess;
}

void OutputResults(shared_ptr<finance::ModelOutput> pOutput,
                   shared_ptr<hg::TheoreticalModel> pModel,
                   shared_ptr<finance::OneTouch> pOneTouch)
{
  std::cout.precision(10);
  std::cout << "The price is: " << pOutput->GetPrice() << std::endl;
  std::cout << "    delta is: " << pOutput->GetDelta() << std::endl;
  std::cout << "    gamma is: " << pOutput->GetGamma() << std::endl;
  std::cout << "    theta is: " << pOutput->GetTheta() << std::endl;
  std::cout << std::endl;
  std::cout << "The spot is: " << pOneTouch->GetSessionData()->GetSpotSharePrice() 
            << std::endl;
  std::cout << "The barrier is: " << pOneTouch->GetBarrier() << std::endl;
  
  // Output sensitivities, if computed
  shared_ptr<hg::NumOutput>
    pNumOutput( dynamic_pointer_cast<hg::NumOutput>( pOutput->GetNumOutput() ) );

  // Output sensitivities, if computed
  if ( pNumOutput->HasSensitivities() )
  {
    std::cout << std::endl;
    std::cout << "Sensitivities:" << std::endl;
    std::vector<double> pdSensitivities = pNumOutput->GetSensitivities();
    for (size_t nIdx = 0; nIdx < pdSensitivities.size(); nIdx++)
      std::cout << "index = " << nIdx << ", value = " << pdSensitivities[nIdx] << std::endl;
  }
  std::cout << std::endl;

  // Compute derivatives by finite differences
  if ( pNumOutput->HasSensitivities() )
  {
    double dShift = 1.e-8;
    
    std::vector<double> 
      sensitivityByFD( hg::ComputeSensitivity(*pModel, *pOneTouch, dShift) );

    std::cout << "sensitivity by finite difference" << std::endl;

    for (size_t nIdx = 0; nIdx < sensitivityByFD.size(); nIdx++)
      std::cout << "index = " << nIdx << ", value = "
                << sensitivityByFD[nIdx] << std::endl;
  }

  bool bTestExtraPoint = false;
  if ( bTestExtraPoint && pOutput->HasPriceAtAnalysisDate() )
  {
    std::cout.precision(12);
    std::cout << "Analysis date data" << std::endl;
    finance::Values pdSpots = pOutput->GetSpotsAtAnalysisDate();
    finance::Values pdPrices = pOutput->GetPricesAtAnalysisDate();
    finance::Values pdDeltas = pOutput->GetDeltasAtAnalysisDate();
    finance::Values pdGammas = pOutput->GetGammasAtAnalysisDate();

    for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
      std::cout << pdSpots[nIdx] 
                << " " << pdPrices[nIdx] 
                << " " << pdDeltas[nIdx] 
                << " " << pdGammas[nIdx] 
                << std::endl;
    std::cout << std::endl;
  }
    
  if ( bTestExtraPoint && pOutput->HasPriceSurface() )
  {
    std::cout << "Surface at valuation date" << std::endl;
    finance::Domain::Spots pdSpots(20);
    pdSpots[19] = pOneTouch->GetBarrier() + 1.0;
    for (size_t nIdxS = 18; nIdxS < 20; nIdxS--)
      pdSpots[nIdxS] = pdSpots[nIdxS + 1] - 0.5;

    finance::SharedSurface priceSurface = pOutput->GetPriceSurface();
    finance::SharedSurface deltaSurface = pOutput->GetDeltaSurface();
    finance::SharedSurface gammaSurface = pOutput->GetGammaSurface();

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
                  << " " << pdPrices[nIdx] 
                  << " " << pdDeltas[nIdx] 
                  << " " << pdGammas[nIdx] 
                  << std::endl;
        
      std::cout << std::endl;
    }
  }
}

void TestOneTouch()
{
  std::cout << "Testing OneTouch" << std::endl;
  std::cout << std::endl;

  shared_ptr<SessionData> pSessionData = InitSessionData();

  shared_ptr<hg::UnderlyingProcess> pUnderlyingProcess = InitProcess();

  shared_ptr<hg::TheoreticalModel> 
    pModel( new hg::TheoreticalModel(pUnderlyingProcess) );

  pModel->SetDebugOutputFile("hg_onetouch.xml");    

  Date maturityDate(2003, Date::May, 1);

  double dUpBarrier = pSessionData->GetSpotSharePrice() + 10.0;
  shared_ptr<OneTouch>
    pOneTouch(new OneTouch(maturityDate, 
                           dUpBarrier, Barrier_UpAndOut, Rebate_Immediate) );

  pOneTouch->SetSessionData(pSessionData);
  
  shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);

  pFlags->ActivateAllSensitivities(true);

  Date analysisDate = pSessionData->GetValuationDate();
  analysisDate.AddMonths(2);

  pFlags->SetAnalysisDate(analysisDate);

  pFlags->SetComputeSurface(false);

  pOneTouch->SetComputationalFlags(pFlags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pOneTouch);

  OutputResults(output, pModel, pOneTouch);

  // ShowConvergence(*pModel, *pOneTouch, 3);

  std::cout << std::endl;
}

void TestFXOneTouch()
{
  std::cout << "Testing FXOneTouch" << std::endl;
  std::cout << std::endl;
  
  shared_ptr<SessionData> pSessionData = InitSessionData();

  shared_ptr<hg::UnderlyingProcess> pUnderlyingProcess = InitProcess();

  shared_ptr<hg::TheoreticalModel> 
    pModel( new hg::TheoreticalModel(pUnderlyingProcess) );

  pModel->SetDebugOutputFile("hg_fxonetouch.xml");

  Date maturityDate(2003, Date::May, 1);

  double dBSBarrierPrice = 0.18;
  double dRefVol = 0.2;
  double dQuote = 0.1;
  shared_ptr<FXOneTouch>
    pFXOneTouch(new FXOneTouch(maturityDate, dBSBarrierPrice,
                               Barrier_UpAndOut, dRefVol) );

  // Exception testing.  These should fail before the quote is set
  //double dP = pFXOneTouch->GetMarketPrice();
  //double dQ = pFXOneTouch->GetMarketQuote();
  //pFXOneTouch->SetMarketPrice(0.3);
  
  pFXOneTouch->SetMarketQuote(dQuote);
  
  // These should work after the quote is set.
  //double dP = pFXOneTouch->GetMarketPrice();
  //double dQ = pFXOneTouch->GetMarketQuote();
  
  pFXOneTouch->SetSessionData(pSessionData);

  shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);

  pFlags->ActivateAllSensitivities(true);

  Date analysisDate = pSessionData->GetValuationDate();
  analysisDate.AddMonths(2);

  pFlags->SetAnalysisDate(analysisDate);

  pFlags->SetComputeSurface(false);

  pFXOneTouch->SetComputationalFlags(pFlags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pFXOneTouch);

  OutputResults(output, pModel, pFXOneTouch);

  // ShowConvergence(*pModel, *pFXOneTouch, 3);

  std::cout << std::endl;
}



int main()
{
  try
  {   
    TestOneTouch();

    TestFXOneTouch();

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
