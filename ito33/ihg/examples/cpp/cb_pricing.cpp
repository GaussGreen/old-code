#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

shared_ptr<ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
    conversionStart = issueDate,
    conversionEnd = maturityDate,
    startCallDate(issueDate),
    endCallDate(maturityDate);

  double
    dIssuePrice = 1,
    dParValue = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  /////////////////////////////////////////////////////////////////////////////
  /// Bond Terms //////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                      maturityDate, dParValue, dRedemptionRate,
                      dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral
                  (issueDate, paymentDates, paymentRates,
                   Date::DayCountConvention_Act365,
                   Frequency_Annual) );
  
  bc->SetCashDistribution(pInterests);

  /////////////////////////////////////////////////////////////////////////////
  /// Conversion Schedule /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  conv->AddConversionPeriod
        ( shared_ptr<ConversionPeriod>
          (new ConversionPeriod(conversionStart, conversionEnd, 1)) );

  /////////////////////////////////////////////////////////////////////////////
  /// create a convertible bond by bond terms and conversion schedule
  shared_ptr<ConvertibleBond> 
    pConvertibleBond( new ConvertibleBond(bc, conv) );


  /////////////////////////////////////////////////////////////////////////////
  //// add call schedule if any
  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
          
  pCallSchedule->AddCallPeriod
                 ( shared_ptr<CallPeriod>
                   ( CallPeriod::CreateWithStrike 
                     (startCallDate, endCallDate, 1) ) );
        
  pConvertibleBond->SetCallSchedule(pCallSchedule);

  
  /////////////////////////////////////////////////////////////////////////////
  /// set exchangeable property if any
  const bool bIsExchangeable = true;
  if ( bIsExchangeable )
  {
    shared_ptr<Issuer> pIssuer(new Issuer);
    std::vector<Date> pDates;
    std::vector<double> pdValues;

    pDates.push_back(maturityDate);
    pdValues.push_back(0.01);

    pIssuer->SetDefaultIntensity(pDates, pdValues);

    pConvertibleBond->SetIssuer(pIssuer);

    pConvertibleBond->SetExchangeable(false);
  }

  /////////////////////////////////////////////////////////////////////////////
  /// attach session data
  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}


shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01");

  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pNumeraire( new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pNumeraire));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.0) ) );
 
  // Setup the rate data, and attach to the session
  shared_ptr<RateData> pRateData( new RateData );
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.04) );
  pRateData->SetYieldCurve(pNumeraire, pyc);
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

//*************************************************************************

int CBPricing()
{
  try {

  bool
    bComputeRho = true,
    bComputeVega = true,
    bComputeSurface = false;
   
  shared_ptr<SessionData> pSessionData = InitSessionData();  

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  pModel->SetDebugOutputFile("ihg_cb.xml");

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<ihg::Volatility>
                         ( new VolatilityFlat(dVol) ) );

  size_t nNbHRTimes = 2;
  static const Date pDates[2] = { 38000, 40000};
  static const double pdValues[2] = { 0.02, .02 };

  pModel->SetHazardRate
          ( shared_ptr<HazardRate>
            (new HazardRateTimeOnly(&pDates[0], &pdValues[0], nNbHRTimes) ) );

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01"));

  shared_ptr<ConvertibleBond> pCB( InitCB(pSessionData) );
  
  pCB->SetComputationalFlags(flags);

  /*
    TheoreticalModel::Compute() return shared_ptr<finance::ModelOutput>
    which is indeed BondLikeOutput. We should do cast here to have acces 
    to more output data.
  */
  shared_ptr<BondLikeOutput> 
    output( static_pointer_cast<BondLikeOutput>(pModel->Compute(*pCB)) );

  std::cout.precision(15);
  std::cout << "\n";
  std::cout << "CB\n";
  std::cout << "\tprice is: " << output->GetPrice() << std::endl;
  if (output->HasRho())
    std::cout << "\tRho is: " << output->GetRho() << std::endl;
  if (output->HasVega())
    std::cout << "\tVega is: " << output->GetVega() << std::endl;
  
  std::cout << "\tBondFloor is: " << output->GetBondFloor() << std::endl;
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

