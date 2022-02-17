

#include "utils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/qualitycontrol.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"

using namespace ito33;
using namespace ito33::finance;


shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData,
       const Date issueDate, const Date maturityDate,
       const shared_ptr<ConversionSchedule>& pConversionSchedule )
{
  double dIssuePrice     = 1;
  double dParValue       = 110;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> bc( new BondTerms(issueDate, dIssuePrice, 
    maturityDate, dParValue, dRedemptionRate, dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral(issueDate, 
                                          paymentDates, paymentRates,
                                          Date::DayCountConvention_Act365, 
                                          Frequency_Annual) );

  bc->SetCashDistribution(pInterests);

  shared_ptr<finance::ConvertibleBond>  
    pConvertibleBond ( new ConvertibleBond(bc, pConversionSchedule) );

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData,
       const Date issueDate,const Date maturityDate)
{
  double dIssuePrice     = 1;
  double dParValue       = 110;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> bc( new BondTerms(issueDate, dIssuePrice, maturityDate,
    dParValue, dRedemptionRate, dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamGeneral(issueDate, 
                                          paymentDates, paymentRates, 
                                          Date::DayCountConvention_Act365,
                                          Frequency_Annual)
              );

  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> pConvSchedule( new ConversionSchedule() );
  
  shared_ptr<ConversionPeriod> 
    pConvPeriod(new ConversionPeriod(issueDate, maturityDate, 1.e-10) );

  pConvSchedule->AddConversionPeriod(pConvPeriod);

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, pConvSchedule) );

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

shared_ptr<SessionData> InitSessionData(const Date issueDate)
{ 
  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pNumeraire));
  shared_ptr<YieldCurve> pYF( new YieldCurveFlat(0.0) ); 
  pEquity->SetBorrowCurve(  pYF );
 
  // Setup the rate data, and attach to the session  
  shared_ptr<YieldCurve> pYC( new YieldCurveFlat(0.04) );
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pNumeraire, pYC);
 
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, issueDate) );

  return pSessionData;
}
