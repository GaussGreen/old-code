
#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/frequency.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/conversionpricereset.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"


#include "ito33/tests/showconvergence.h"

using namespace ito33;
using namespace ito33::finance;

//Based on the prospectus called: Transmile 1%25 2010 OID.pdf
void PriceTransmile()
{
  // stock underlying  
  shared_ptr<Numeraire> pforeignCurrency(new Numeraire("RM"));
  // Bond  
  shared_ptr<Numeraire> pbaseCurrency(new Numeraire("USD"));
  
  double dFixedFXRate      = 3.80; //3.80RM = 1 US

  //today
  Date valuationDate(2006, Date::Feb, 2);
  
  //page 1 prospectus
  Date issueDate(2005, Date::May, 17);

  //page 1 prospectus
  Date maturityDate(2010, Date::May, 17);

  double dSpotPrice = 9.40; //in RM first page prospectus

  shared_ptr<Equity> pEquity( new Equity(dSpotPrice, pforeignCurrency) );

  //make a yield curve
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.01)) );
  
  //Maybe not dividends after all
  shared_ptr<finance::Dividends> pDiv( new Dividends() );
  pDiv->Add(finance::Dividend::Type::Cash, Date(2004, Date::Dec, 31), .03);
  pDiv->Add(finance::Dividend::Type::Cash, Date(2005, Date::Dec, 31), .033);
  pDiv->Add(finance::Dividend::Type::Cash, Date(2006, Date::Dec, 31), .036);
  pDiv->Add(finance::Dividend::Type::Cash, Date(2007, Date::Dec, 31), .040);
  pDiv->Add(finance::Dividend::Type::Cash, Date(2008, Date::Dec, 31), .044);
  pDiv->Add(finance::Dividend::Type::Cash, Date(2009, Date::Dec, 31), .048);

  //pEquity->SetDividends(pDiv);
 
  // Setup the issuer, and attach to the session
  shared_ptr<YieldCurve> pYC( new YieldCurveFlat(.03) );
  shared_ptr<YieldCurve> pBondYC( new YieldCurveFlat(.02) );  
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pforeignCurrency, pYC);
  pRateData->SetYieldCurve(pbaseCurrency, pBondYC);

  shared_ptr<SpotFXRates> pFX(new SpotFXRates() );
  pFX->SetFXRate( pforeignCurrency, pbaseCurrency, 1./dFixedFXRate);
  pRateData->SetSpotFXRates(pFX);

  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );
  
  double dIssuePrice      = 1.0;
  double dParValue        = 100; //in effect but scaled down 100,000.0;

  //dRedemptionPrice the redemption price, as a percentage of nominal
  double dRedemptionRate = 1.2241;
  double dRecoveryRate   = 0.; 
  double dRate           = 1./100.;
  Date firstDate(2005, Date::Nov, 17);

  shared_ptr<BondTerms> bc
    ( 
      new BondTerms(issueDate, dIssuePrice, maturityDate, 
                    dParValue, dRedemptionRate, dRecoveryRate) 
     );

  shared_ptr<CashFlowStream> pInterests
    ( 
     new CashFlowStreamUniform( issueDate, 
                                firstDate,
                                maturityDate,
                                dRate,
                                Date::DayCountConvention_Act365,
                                finance::Frequency_SemiAnnual)
     );

  bc->SetCashDistribution(pInterests);

  // bc->SetAccretingBond(.05); //page 79
  
  
  //Call Period pp77 Sec 8.
  double dStrike = 1.0;
  double dTrigger = 1.3;
  Date startCallDate(2008, Date::May, 17);
  Date endCallDate = maturityDate; endCallDate.AddDays(-20);
  
  shared_ptr<CallPeriod> pCallPeriod( 
    CallPeriod::CreateWithStrike(startCallDate, endCallDate, dStrike) );
  pCallPeriod->SetTrigger( dTrigger);

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->AddCallPeriod(  pCallPeriod );

  //Put Period section 8.6
  Date putDate(2008, Date::May, 17);
  double dPutStrike = 1.;

  shared_ptr<PutSchedule> pPutSchedule( new PutSchedule() );
  pPutSchedule->AddPutWithStrike(putDate, dPutStrike);


  //Reset conversion schedule
  Date resetDate(2006, Date::May, 17);
  double dInitialConvPrice = 10.81;
  double dCurrentConvPrice = dInitialConvPrice;


  shared_ptr<ConversionPriceReset> pConvPriceReset
    ( new ConversionPriceReset(resetDate, .8) );

  Date convStartDate = Date(2005, Date::Jun, 27 );
  Date convEndDate   = Date(2010, Date::May, 3 );

  shared_ptr<ResetConversionSchedule> pResetConvSchedule
    ( 
     new ResetConversionSchedule(convStartDate, convEndDate, dInitialConvPrice, 
                      dCurrentConvPrice, ResetFlooredBy_InitialConversionPrice) 
    );

  pResetConvSchedule->AddConversionPriceReset( pConvPriceReset );

  //Create the reset
  shared_ptr<Reset> pReset( new Reset(bc, pResetConvSchedule) );
  
  pReset->SetNumeraire(pbaseCurrency);

  pReset->SetSessionData(pSessionData);
    
  pReset->SetCallSchedule( pCallSchedule );

  pReset->SetPutSchedule( pPutSchedule );

  //Note it is not clear that this is the actual case
  //It is simply a test.
  pReset->SetTriggerInCurrencyOf(TriggerInCurrencyOf_Underlying);
  pReset->SetFixedFXRate( 1./dFixedFXRate );

  double dVol = 0.5;
  shared_ptr<ihg::Volatility> pVolatility(new ihg::VolatilityFlat(dVol) );
  double dLambda = .02;
  shared_ptr<ihg::HazardRate> pHazardRate(new ihg::HazardRateFlat(dLambda) );

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
  pModel->SetVolatility( pVolatility );
  pModel->SetHazardRate( pHazardRate );

  // Actually price
  pModel->SetDebugOutputFile("./ihg_transmile.xml");
  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pReset);

  std::cout.precision(12);
  std::cout << "Price: " << output->GetPrice() << std::endl;
  std::cout << "Delta: " << output->GetDelta() << std::endl;
  std::cout << "Gamma: " << output->GetGamma() << std::endl;

 // ShowConvergence(*pModel, *pReset);
}