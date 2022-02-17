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
#include "ito33/finance/computationalflags.h"

#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOneTouch);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

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

shared_ptr<ihg::UnderlyingProcess> InitProcess()
{
  double dVol = 0.2;
  shared_ptr<Volatility> pVolatility( new VolatilityFlat(dVol) );

  double dLambda = .02;
  shared_ptr<HazardRate> pHazardRate(new ihg::HazardRateFlat(dLambda));

  return make_ptr( new ihg::UnderlyingProcess(pVolatility, pHazardRate) );
}

void OutputResults(shared_ptr<finance::ModelOutput> pOutput,
                   shared_ptr<ihg::TheoreticalModel> pModel,
                   shared_ptr<finance::OneTouch> pOneTouch)
{
  std::cout.precision(10);
  std::cout << "The price is: " << pOutput->GetPrice() << std::endl;
  std::cout << "    delta is: " << pOutput->GetDelta() << std::endl;
  std::cout << "    gamma is: " << pOutput->GetGamma() << std::endl;
  std::cout << "    theta is: " << pOutput->GetTheta() << std::endl;
  std::cout << std::endl;
  std::cout << "The spot is: " 
            << pOneTouch->GetSessionData()->GetSpotSharePrice() << std::endl;
  std::cout << "The barrier is: " << pOneTouch->GetBarrier() << std::endl;
}

void TestOneTouch()
{
  std::cout << "Testing OneTouch" << std::endl;
  std::cout << std::endl;

  shared_ptr<SessionData> pSessionData = InitSessionData();

  shared_ptr<ihg::UnderlyingProcess> pUnderlyingProcess = InitProcess();

  shared_ptr<ihg::TheoreticalModel> 
    pModel( new ihg::TheoreticalModel(pUnderlyingProcess) );

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

  shared_ptr<ihg::UnderlyingProcess> pUnderlyingProcess = InitProcess();

  shared_ptr<ihg::TheoreticalModel> 
    pModel( new ihg::TheoreticalModel(pUnderlyingProcess) );

  // pModel->SetDebugOutputFile("ihg_fxonetouch.xml");

  Date maturityDate(2003, Date::May, 1);

  double dBSBarrierPrice = 0.18;
  double dRefVol = 0.2;
  double dQuote = 0.1;
  shared_ptr<FXOneTouch>
    pFXOneTouch(new FXOneTouch(maturityDate, dBSBarrierPrice,
                               Barrier_UpAndOut, dRefVol) );

  pFXOneTouch->SetMarketQuote(dQuote);
  
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
