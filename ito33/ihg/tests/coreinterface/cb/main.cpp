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
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/frequency.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_MODULE(IHGPriceCB);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

static Date valuationDate;
static shared_ptr<SpotFXRates> fx;
shared_ptr<SessionData> pSessionData;

shared_ptr<finance::ConvertibleBond> InitCB()
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
    conversionStart = Date(2003, Date::Apr, 1),
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
    Date::DayCountConvention_Act365,Frequency_SemiAnnual) );
  
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
  
  size_t nNoticePeriod = 20;

//    pCallSchedule->SetNoticePeriod(nNoticePeriod);
//    startCallDate.AddMonths(4);
//    endCallDate.AddMonths(-4);
        
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithYield
                             (startCallDate, endCallDate, .2) );

  pCallSchedule->AddCallPeriod( pCallPeriod);
        
  // pConvertibleBond->SetCallSchedule(pCallSchedule);

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  valuationDate = Date("2003/02/01");

  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.0) ) );
 
  // Setup the rate data and attach to the session data
  shared_ptr<RateData> pRateData( new RateData);
  
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.04) );  
 
  pRateData->SetYieldCurve(pCurrency, pyc);
 
  pSessionData = shared_ptr<SessionData>(new SessionData(pRateData, pEquity, 
                                                        valuationDate));

  return pSessionData;
}

//*************************************************************************

int main()
{
  try {
   

  bool
    bComputeRho = true,
    bComputeVega = true,
    bComputeSurface = true;
   
  InitSessionData();  

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  pModel->SetDebugOutputFile("./ihgcb.xml");

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
  pModel->SetHazardRate( make_ptr(new ihg::HazardRateFlat(dLambda) ) );

  shared_ptr<finance::ConvertibleBond> pCB = InitCB();
  
  shared_ptr<finance::ComputationalFlags> flags(new finance::ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01"));

  pCB->SetComputationalFlags(flags); 

  shared_ptr<finance::BondLikeOutput> 
    output( static_pointer_cast<finance::BondLikeOutput>
            ( pModel->Compute(*pCB) ) );

  std::cout.precision(15);
  std::cout << "The cb price is: " << output->GetPrice() << std::endl;
  if (output->HasRho())
    std::cout << "The cb Rho is: " << output->GetRho() << std::endl;
  if (output->HasVega())
    std::cout << "The cb Vega is: " << output->GetVega() << std::endl;
  
  //ShowConvergence(*pModel, *pCB);

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

