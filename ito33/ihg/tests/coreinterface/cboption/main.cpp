/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/coreinterface/cboption/main.cpp
// Purpose:     To test the financial interface of the cb option
// Author:      Nabil Ouachani
// Created:     2005/07/06
// RCS-ID:      $Id: main.cpp,v 1.10 2006/08/20 09:49:26 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

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
#include "ito33/finance/floatingrates.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/cboption.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/ihg/bondlikeoutput.h"
#include "ito33/ihg/cboptionoutput.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCBOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;

shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2000, Date::Feb, 1),
    firstCouponDate = Date(2000, Date::Nov, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
    //conversionStart = Date(2003, Date::Apr, 1),
    conversionStart = Date(2003, Date::Feb, 1),
    conversionEnd = Date(2004, Date::May, 1);

  double
    dIssuePrice = 1,
    dParValue = 110,
    dRedemptionRate = 1.0,
    dRecoveryRate = 0.;

  shared_ptr<BondTerms> bc ( new BondTerms(issueDate, dIssuePrice,
                                maturityDate, dParValue, dRedemptionRate,
                                dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamUniform(issueDate, firstCouponDate, 
    maturityDate, 0.02, Date::DayCountConvention_Act365, Frequency_SemiAnnual) );
  
  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> conv ( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  conv->AddConversionPeriod(
    shared_ptr<ConversionPeriod> 
    (new ConversionPeriod(conversionStart, conversionEnd, 1)) );

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond ( new finance::ConvertibleBond(bc, conv) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  shared_ptr<CallSchedule> pCallSchedule ( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(true);
  
//  size_t nNoticePeriod = 20;

//    pCallSchedule->SetNoticePeriod(nNoticePeriod);
//    startCallDate.AddMonths(4);
//    endCallDate.AddMonths(-4);
        
  shared_ptr<CallPeriod> 
    pCallPeriod ( CallPeriod::CreateWithStrike
                              (startCallDate, endCallDate, 1.05) );

  pCallSchedule->AddCallPeriod( pCallPeriod);

        
  pConvertibleBond->SetCallSchedule(pCallSchedule);

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}


shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/10");

  // Setup the equity, and attach to session
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve> (new YieldCurveFlat(0.0)) );
 
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04));

  pRateData->SetYieldCurve(pCurrency, pyc);  
  
  shared_ptr<SessionData> 
    pSessionData (new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<finance::CBOption> 
InitCBOption(const shared_ptr<ConvertibleBond>& pCB, 
             const shared_ptr<SessionData>& pSessionData)
{
  double 
    dRecallSpread = 0.01;

  Date 
    maturityDate("2004/05/01"),
    startOfAccruedDate("2003/02/01"), 
    firstUnknownCouponDate("2003/08/01"),
    lastButOneUnknownCouponDate("2004/02/01");
  
  Date::DayCountConvention dcc = Date::DayCountConvention_Act365;
  
  Frequency frequency = Frequency_Quarterly;

  shared_ptr<finance::FloatingRates>
    pFloatingRates ( new finance::FloatingRates
                                  (
                                    dRecallSpread, 
                                    startOfAccruedDate, 
                                    firstUnknownCouponDate, 
                                    maturityDate, frequency
                                  ) 
                   );
  
  pFloatingRates->SetDayCountConvention(dcc);

  vector<Date> pknownPaymentDates;
  vector<double> pknownPaymentAmounts;
  pknownPaymentDates.push_back("2003/05/01");
  pknownPaymentAmounts.push_back(0.01);

  pFloatingRates->SetKnownPaymentStream(pknownPaymentDates, 
    pknownPaymentAmounts);

  shared_ptr<finance::CBOption> 
    pCBOption ( new finance::CBOption(pCB,
                                      pFloatingRates,
                                      maturityDate) 
              );
  
  pCBOption->SetSessionData(pSessionData);
  
  return pCBOption;
}

//*************************************************************************

int main()
{
  try {
   

  bool
    bComputeRho = true,
    bComputeVega = true,
    bComputeSurface = false;
   
  shared_ptr<SessionData> pSessionData = InitSessionData();  

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  double dVol = 0.3;
  pModel->SetVolatility( shared_ptr<Volatility> (new VolatilityFlat(dVol)) );

  /*
  size_t nNbHRTimes = 2;
  static const Date pTimes[2] = { 38000, 40000};
  static const double pdValues[2] = { 0.02, .02 };

  pModel->SetHazardRate
          ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );
  */
   
  double dLambda = 0.; 
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  shared_ptr<finance::ConvertibleBond> pCB = InitCB(pSessionData);
  
  shared_ptr<finance::ComputationalFlags> flags(new finance::ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01"));

  shared_ptr<finance::CBOption> pCBOption = InitCBOption(pCB, pSessionData);

  //shared_ptr<ihg::ModelOutput> cboutput = pModel->Compute(*pCB);

  pCBOption->SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCBOption);
  
  shared_ptr<finance::CBOptionOutput>
    cboptionoutput( dynamic_pointer_cast<finance::CBOptionOutput>(output) );

  shared_ptr<finance::BondLikeOutput>
    cboutput = cboptionoutput->GetCBOutput();

  std::cout.precision(15);

  //CB Option output
  std::cout << "CB option output " << std::endl;
  std::cout << "The cb option price is: " 
            << cboptionoutput->GetPrice() << std::endl;
  std::cout << "The cb option delta is: " 
            << cboptionoutput->GetDelta() << std::endl;
  std::cout << "The cb option gamma is: " 
            << cboptionoutput->GetGamma() << std::endl;
  std::cout << "The cb option theta is: " 
            << cboptionoutput->GetTheta() << std::endl;
  
  if ( cboptionoutput->HasRho() )
    std::cout << "The cb option rho is: " 
              << cboptionoutput->GetRho() << std::endl;
  if ( cboptionoutput->HasVega() )
    std::cout << "The cb option Vega is: " 
              << cboptionoutput->GetVega() << std::endl;

  std::cout << std::endl << "Strike is : " 
            << pCBOption->ComputeStrike() << "\n";
  std::cout <<"________________________________________\n\n";

  // CB output
  std::cout << "CB output " << std::endl;
  std::cout << "The cb price is: " << cboutput->GetPrice() << std::endl;
  std::cout << "The cb delta is: " << cboutput->GetDelta() << std::endl;
  std::cout << "The cb gamma is: " << cboutput->GetGamma() << std::endl;
  std::cout << "The cb theta is: " << cboutput->GetTheta() << std::endl;
  std::cout << "The cb bond floor is: " << cboutput->GetBondFloor() 
            << std::endl;
  if ( cboutput->HasRho() )
    std::cout << "The cb Rho is: " << cboutput->GetRho() << std::endl;
  if ( cboutput->HasVega() )
    std::cout << "The cb Vega is: " << cboutput->GetVega() << std::endl;
  
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
