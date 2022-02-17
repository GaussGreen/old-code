

//                          Carnival 1.132%2033 
// using a swap curve

#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_swap.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/option.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
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

shared_ptr<YieldCurveSwap> InitSwapCurve(Date valuationDate)
{
  shared_ptr<YieldCurveSwap> pyc(new YieldCurveSwap(valuationDate));

  CashRates cashRates;
  cashRates.SetBasis(Date::DayCountConvention_30360);
  cashRates.AddLeg(0.05, 1, TimeUnit_Day);
  cashRates.AddLeg(0.0502, 7, TimeUnit_Day);
  cashRates.AddLeg(0.0513, 1, TimeUnit_Month);
  cashRates.AddLeg(0.0523, 3, TimeUnit_Month);
  cashRates.AddLeg(0.054, 6, TimeUnit_Month);
  cashRates.AddLeg(0.0562, 1, TimeUnit_Year);
  
  SwapRates swapRates;
  swapRates.SetBasis(Date::DayCountConvention_Act365);
  swapRates.AddLeg(0.0571, 2, TimeUnit_Year, Frequency_SemiAnnual);
  swapRates.AddLeg(0.058, 3, TimeUnit_Year, Frequency_SemiAnnual);
  swapRates.AddLeg(0.061, 4, TimeUnit_Year, Frequency_SemiAnnual);

  pyc->SetSwapRates(swapRates);
  pyc->SetCashRates(cashRates);

  return pyc;

}


shared_ptr<AttachedWarrantConvertibleBond> 
  InitCarnivalTest()
{

  // The session data variables
  Date valuationDate("2004/01/01");

  double dSpot = 646.88 / 12.18;
  shared_ptr<Numeraire> pNumeraire( new Numeraire("EUR") );
  shared_ptr<Equity> pEquity(new Equity(dSpot, pNumeraire));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.0) ) );

  shared_ptr<RateData> pRateData( new RateData );
  shared_ptr<YieldCurve> pyc = InitSwapCurve(valuationDate);
  pRateData->SetYieldCurve(pNumeraire, pyc);
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );


  // Setup the bond terms
  Date issueDate = Date(2003, Date::Apr, 29);
  Date maturityDate = Date(2033, Date::Apr, 29);

  double dIssuePrice     = 646.88/1000; //issue price percentage of nominal
  double dParValue       = 1000;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms>
    bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
                      dParValue, dRedemptionRate, dRecoveryRate) );

  Date firstCashDate = Date(2003, Date::Oct, 29);
  Date lastCashDate = Date(2008, Date::Apr, 29);

  double dCashYield = 1.132/100.;
  shared_ptr<CashFlowStreamUniform> 
      pInterests( new CashFlowStreamUniform
                         (
                           issueDate,
                           firstCashDate, 
                           lastCashDate,
                           dCashYield,
                           Date::DayCountConvention_30360,
                           Frequency_SemiAnnual
                         )
                   );

  bc->SetCashDistribution(pInterests);

  double dOIDYield = 1.75/100.;
  bc->SetCashPayToZero(dOIDYield);
  bc->SetYieldCompoundingFrequency(Frequency_SemiAnnual);

  // Conversion provisions
  Date conversionStart = Date(2003,Date::Aug,31);
  Date conversionEnd = maturityDate;
  Date resetDate = Date(2008, Date::Apr, 29);

  double dInitialRatio = 12.18;
  double dShareFactor  = 11.3258;
  double dStrike = dIssuePrice * dParValue / dInitialRatio;
  
  shared_ptr<ShareDependentConversion> 
    pConv( new ShareDependentConversion
               (conversionStart, conversionEnd, dInitialRatio, dShareFactor) );

  pConv->SetResetDate(resetDate);
  pConv->SetFixedStrike(dStrike);

  pConv->SetKeepAccrued(false);
  pConv->SetForfeitCoupon(false);
   
  double dTriggerRate = 120./100.;
  double dChangeRate  = 0.0;
  double dExtremeTriggerRate = dTriggerRate;
  bool bIsCurrentlyActive = true;
  pConv->SetCoCo(dTriggerRate, 
                 CoCoType_CheckQuarterlyAndConvertDuringNextQuarter,
                 dChangeRate, dExtremeTriggerRate, bIsCurrentlyActive);

  // Call schedule
  Date callStart = Date(2008, Date::Apr, 29);
  Date callEnd = maturityDate;

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
    
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike(callStart, callEnd, 1.) );

  pCallSchedule->AddCallPeriod( pCallPeriod );

  // Put provisions
  shared_ptr<PutSchedule> pPutSchedule(new PutSchedule);
  pPutSchedule->AddPutWithStrike(Date(2008, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2013, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2018, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2023, Date::Apr, 29), 1.0);
  pPutSchedule->AddPutWithStrike(Date(2028, Date::Apr, 29), 1.0);

  // Construct the attached warrant cb
  shared_ptr<AttachedWarrantConvertibleBond> 
    pAttachedWarrantCB( new AttachedWarrantConvertibleBond(bc,pConv) );

  // Finalize the contract (call schedule, sessiondata, etc)
  pAttachedWarrantCB->SetCallSchedule(pCallSchedule);
  pAttachedWarrantCB->SetPutSchedule(pPutSchedule);
  pAttachedWarrantCB->SetConversionTriggerAsPercentageOf
                      ( TriggerAsPercentageOf_IssuePrice );
  pAttachedWarrantCB->SetSessionData(pSessionData);

  return pAttachedWarrantCB;
}



//*************************************************************************

int AttachedWarrantCBPricing()
{
  try {

  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;
   
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  pModel->SetDebugOutputFile("ihg_attachedwarrantcb.xml");

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol) ) );

  size_t nNbHRTimes = 2;
  static const Date pDates[2] = { 38000, 40000};
  static const double pdValues[2] = { 0.02, .02 };

  pModel->SetHazardRate
          ( shared_ptr<HazardRate>
            ( new HazardRateTimeOnly(&pDates[0], &pdValues[0], nNbHRTimes) ) );

  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01"));

  shared_ptr<AttachedWarrantConvertibleBond>
    pCB = InitCarnivalTest();
   
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
  std::cout << "Attached Warrant Convertible Bond: (Carnival 1.132% 2003)\n";
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

