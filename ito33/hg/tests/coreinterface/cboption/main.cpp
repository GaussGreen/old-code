/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/coreinterface/cboption/main.cpp
// Purpose:     To test the financial interface of the cb option
// Created:     2006/01/19
// RCS-ID:      $Id: main.cpp,v 1.9 2006/08/19 23:47:17 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/termstructurecds.h"
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

#include "ito33/hg/bondlikeoutput.h"
#include "ito33/hg/cboptionoutput.h"
#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceCBOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

using namespace std;

shared_ptr<ConvertibleBond> InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2000, Date::Feb, 1),
    firstCouponDate = Date(2000, Date::Nov, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
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
                                          maturityDate, 0.02,
                                          Date::DayCountConvention_Act365, 
                                          Frequency_SemiAnnual) );
  
  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> conv ( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  conv->AddConversionPeriod(
    shared_ptr<ConversionPeriod> 
    (new ConversionPeriod(conversionStart, conversionEnd, 1)) );

  shared_ptr<ConvertibleBond> 
    pConvertibleBond ( new ConvertibleBond(bc, conv) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  shared_ptr<CallSchedule> pCallSchedule ( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(true);
  
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

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  // Setup the equity, and attach to session
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve> (new YieldCurveFlat(0.0)) );
 
  // Setup the issuer, and attach to the session
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, 
    shared_ptr<YieldCurve> (new YieldCurveFlat(0.04) ) );
  
  shared_ptr<SessionData> 
    pSessionData (new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<CBOption> 
InitCBOption(const shared_ptr<ConvertibleBond>& pCB, 
             const shared_ptr<SessionData>& pSessionData)
{
  double 
    dRecallSpread = 0.01;

  Date 
    maturityDate("2004/05/01"),
    startOfAccruedDate("2003/02/21"), 
    firstUnknownCouponDate("2003/08/01"),
    lastButOneUnknownCouponDate("2004/02/01");
  
  Date::DayCountConvention dcc = Date::DayCountConvention_Act365;
  
  Frequency frequency = Frequency_Quarterly;

  shared_ptr<FloatingRates>
    pFloatingRates( new FloatingRates
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

  pFloatingRates->SetKnownPaymentStream(pknownPaymentDates, pknownPaymentAmounts);

  shared_ptr<CBOption> 
    pCBOption( new CBOption(pCB, pFloatingRates, maturityDate) );
  
  pCBOption->SetSessionData(pSessionData);
  
  return pCBOption;
}

//*************************************************************************

int main()
{
  try {
   
  bool
    bComputeRho = true,
    bComputeSurface = false;
   
  shared_ptr<SessionData> pSessionData = InitSessionData();  

  size_t nNbRegimes = 2;

  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.1);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.02);

  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess( new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  Jumps jumps;

  jumps.push_back( Jump(0.1, -0.2) );
  pUnderlyingProcess->SetJumps(0, 0, jumps);
  
  jumps.clear();
  jumps.push_back( Jump(0.15, -0.25) );
  pUnderlyingProcess->SetJumps(0, 1, jumps); 

  jumps.clear();
  jumps.push_back( Jump(0.2, -0.5) );
  pUnderlyingProcess->SetJumps(1, 1, jumps); 

  jumps.clear();
  jumps.push_back( Jump(0.11, 0.35) );
  pUnderlyingProcess->SetJumps(1, 0, jumps); 
  /**/

  shared_ptr<hg::TheoreticalModel> 
    pModel( new hg::TheoreticalModel(pUnderlyingProcess) );
  
  // pModel->EnableDebugOutput();

  shared_ptr<ConvertibleBond> pCB = InitCB(pSessionData);
  
  shared_ptr<CBOption> pCBOption = InitCBOption(pCB, pSessionData);

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01") );

  pCBOption->SetComputationalFlags(flags);

  std::cout.precision(15);
  
  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCBOption);
  
  shared_ptr<hg::CBOptionOutput>
    cboptionoutput( dynamic_pointer_cast<hg::CBOptionOutput>(output) );

  shared_ptr<hg::BondLikeOutput>
    cboutput( dynamic_pointer_cast<hg::BondLikeOutput> 
              ( cboptionoutput->GetCBOutput() ) );

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

  std::cout << std::endl << "Strike is : " 
            << pCBOption->ComputeStrike() << "\n";
  std::cout <<"________________________________________\n\n";

  // CB output
  std::cout << "CB output " << std::endl;
  std::cout << "The cb price is: " << cboutput->GetPrice() << std::endl;
  std::cout << "The cb delta is: " << cboutput->GetDelta() << std::endl;
  std::cout << "The cb gamma is: " << cboutput->GetGamma() << std::endl;
  std::cout << "The cb theta is: " << cboutput->GetTheta() << std::endl;
  std::cout << "The cb bond floor is: " << cboutput->GetBondFloor() << std::endl;
  if ( cboutput->HasRho() )
    std::cout << "The cb Rho is: " << cboutput->GetRho() << std::endl;

  // ShowConvergence(*pModel, *pCBOption);

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
