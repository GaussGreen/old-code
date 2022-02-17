#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/impliededsspreads.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceEDS);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

void RunImpliedSpreads()
{
 
  Date valuationDate(2002, Date::Jan, 1);

  double dS0 = 100.0;
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dS0, pCurrency));

  shared_ptr<Dividends> pDividends(new Dividends());
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.00));
    
  pEquity->SetBorrowCurve(pyf);
    
  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));  
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.2;
  shared_ptr<Volatility> pVolatility(new VolatilityFlat(dVol));
  pModel->SetVolatility( pVolatility );

  /*
  size_t nNbHRTimes = 2;
  static const Date pTimes[2] = { 38000, 40000};
  static const double pdValues[2] = { 0.02, .02 };

  pModel->SetHazardRate
          ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );
  */

  double dAlpha = 0.02;
  double dBeta = 1.2;
  shared_ptr<HazardRate> pHR(new ihg::HazardRatePower(dAlpha, dBeta, dS0));
  pModel->SetHazardRate( pHR );

  double dRecoveryRate = 0.0;
  double dBarrier = 30;

  Date issueDate(2002, Date::Jan, 1);
  Date firstSpreadDate(2002, Date::Jul, 1);
  Frequency freq = Frequency_SemiAnnual;
  
  Date lastSpreadDate(2003, Date::Jan, 1);
  size_t nNbSpreads = 4;
  std::vector<Date> maturityDates(nNbSpreads);
  for (size_t nIdx = 0; nIdx < nNbSpreads; nIdx++)
    maturityDates[nIdx] = lastSpreadDate.AddMonths(12);
    //maturityDates[nIdx] = lastSpreadDate.AddMonths(12 * nIdx);

  shared_ptr<ImpliedEDSSpreads>
    pImpliedSpreads( new ImpliedEDSSpreads
                          (
                            issueDate,
                            firstSpreadDate,
                            Date::DayCountConvention_Act365,
                            freq,
                            dBarrier,
                            dRecoveryRate
                          )
                    );

  std::vector<double>
    spreads = pImpliedSpreads->Compute(pModel, pSessionData, maturityDates);

  std::cout << "issue date: " << issueDate << std::endl;
  std::cout << "first spread: " << firstSpreadDate << std::endl;
  std::cout << "payment frequency: " << freq << std::endl;
  std::cout << "barrier: " << dBarrier << std::endl;
  std::cout << "recovery rate: " << dRecoveryRate << std::endl;
  std::cout << std::endl;

  for (size_t nIdx = 0; nIdx < nNbSpreads; nIdx++)
  {
    shared_ptr<CashFlowStreamUniform>
      pSpreadStream( new CashFlowStreamUniform
                          (
                            issueDate,
                            firstSpreadDate,
                            maturityDates[nIdx],
                            spreads[nIdx],
                            Date::DayCountConvention_Act365,
                            Frequency_SemiAnnual
                          )
                    );

    shared_ptr<EDS>
      pEDS( new EDS(dRecoveryRate, pSpreadStream, dBarrier) );

    pEDS->SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pEDS);

    std::cout.precision(10);
    std::cout << "maturity = " << maturityDates[nIdx]
              << ", price = " << output->GetPrice() 
              << ", spread = " << spreads[nIdx] 
              << std::endl;
  }

  std::cout << endl;
}
