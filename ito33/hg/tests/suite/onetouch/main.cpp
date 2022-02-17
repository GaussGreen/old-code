#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"

#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/jumps.h"

#include "hg/computesensitivity.h"
#include "hg/tests/onetouch.h"

#ifdef ITO33_TEST_MODENV
#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#endif

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceOneTouch);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

shared_ptr<SessionData> InitSessionData()
{
  Date valuationDate(2002, Date::Jul, 1);

  double dRate = 0.05;
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  shared_ptr<RateData> pRateData(new RateData );
  pRateData->SetYieldCurve(pCurrency, pyc );

  double dSpot = 30.0;
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));

  shared_ptr<Dividends> pDividends( new Dividends() );
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.02) );
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<SessionData> pSessionData(
    new SessionData(pRateData, pEquity, valuationDate));

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
        pdIntensity[0] = 0.1 * (nIdx1 + 1);
        pdAmplitude[0] = 0.02 * nIdx2 - 0.012;

        if (nIdx1 != nIdx2)
          pUnderlyingProcess->SetJumps(nIdx1, nIdx2, pdIntensity, pdAmplitude);
      }
    }
  } // if more than one regime

  return pUnderlyingProcess;
}

void TestOneTouch()
{
  std::cout << "Testing OneTouch" << std::endl;
  std::cout << std::endl;

  shared_ptr<SessionData> pSessionData = InitSessionData();

  shared_ptr<hg::UnderlyingProcess> pUnderlyingProcess = InitProcess();

  shared_ptr<hg::TheoreticalModel> 
    pModel( new hg::TheoreticalModel(pUnderlyingProcess) );

  Date maturityDate(2003, Date::May, 1);

  double dUpBarrier = pSessionData->GetSpotSharePrice() + 10.0;
  shared_ptr<OneTouch>
    pUpAndOut(new OneTouch(maturityDate, 
                           dUpBarrier, Barrier_UpAndOut, Rebate_Immediate) );

  pUpAndOut->SetSessionData(pSessionData);

  OneTouchTest::Setup(pModel, pUpAndOut);

  CppUnit::TextUi::TestRunner runner;

  runner.addTest(OneTouchTest::OneTouchTest::suite());

  runner.run("");

  double dDownBarrier = pSessionData->GetSpotSharePrice() - 10.0;
  shared_ptr<OneTouch>
    pDownAndOut(new OneTouch
        (maturityDate, dDownBarrier, Barrier_DownAndOut, Rebate_Immediate) );

  pDownAndOut->SetSessionData(pSessionData);

  OneTouchTest::Setup(pModel, pDownAndOut);

  CppUnit::TextUi::TestRunner runner2;

  runner2.addTest(OneTouchTest::OneTouchTest::suite());

  runner2.run("");
}

int main()
{
  try
  {   
    TestOneTouch();

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
