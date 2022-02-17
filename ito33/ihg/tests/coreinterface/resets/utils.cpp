

#include "utils.h"

#include "ito33/beforestd.h"
#include <cmath>
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/frequency.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/ResetConversionSchedule.h"


#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/finance/computationalflags.h"


#include "ihg/model.h"
#include "ito33/pricing/model.h"

#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/extrapolationmode.h"
#include "ito33/numeric/interpolation.h"


using namespace ito33;
using namespace ito33::numeric;
using namespace ito33::pricing;
using namespace ito33::finance;
using namespace ito33::ihg;

namespace ito33
{

//--------------------------------------------------------------------

shared_ptr<finance::ConvertibleBond> 
InitCB
  (
    const shared_ptr<SessionData>& pSessionData, 
    Date issueDate, Date maturityDate,
    const shared_ptr<ConversionSchedule>& pConversionSchedule, double dParValue
  )
{
  double
    dIssuePrice = 1,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

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
                    Date::DayCountConvention_Act365,finance::Frequency_Annual)
              );

  bc->SetCashDistribution(pInterests);

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, pConversionSchedule) );

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}

shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData, 
       Date issueDate, Date maturityDate, 
       double dParValue)
{
  double
    dIssuePrice = 1,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

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
                    Date::DayCountConvention_Act365,finance::Frequency_Annual)
              );

  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> pConvSchedule( new ConversionSchedule() );
  pConvSchedule->AddConversionPeriod(shared_ptr<ConversionPeriod>(
    new ConversionPeriod(issueDate, maturityDate, 1.0) ) );

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, pConvSchedule) );

  pConvertibleBond->SetSessionData(pSessionData);

  return pConvertibleBond;
}


shared_ptr<finance::Reset> 
InitReset(const shared_ptr<SessionData>& pSessionData,
               Date issueDate,
               Date maturityDate,
               const shared_ptr<ResetConversionSchedule>& pResetConversionSchedule,
               double dParValue)
{
  double
    dIssuePrice = 1,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;  

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
                   Date::DayCountConvention_Act365,finance::Frequency_Annual)
             );

  bc->SetCashDistribution(pInterests);

  shared_ptr<finance::Reset> 
    pReset( new finance::Reset(bc, pResetConversionSchedule) );

  pReset->SetSessionData(pSessionData);

  return pReset;
}

shared_ptr<SessionData> InitSessionData(Date issueDate,
                                   bool bFakeDividend,Date valuationDate)
{
  // Setup the equity, and attach to session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
  
  if ( bFakeDividend )
  {
    shared_ptr<finance::Dividends> pDiv( new Dividends() );
    pDiv->Add(finance::Dividend::Type::Cash, Date(1900,Date::Jun,1),1.);
    pEquity->SetDividends(pDiv);
  }

  // Setup the rate data, and attach to the session data
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(0.04));
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity,
                                 valuationDate.IsValid() ? valuationDate : 
                                 issueDate) );
  return pSessionData;
}

}
