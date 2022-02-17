#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessionData.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/yieldcurve_flat.h"
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

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/showconvergence.h"

#include "ihg/tests/testutils.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCB);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

using namespace std;


shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  double
    dIssuePrice = 1,
    dParValue = 110,
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
                     Date::DayCountConvention_Act365,
                     Frequency_Annual)
              );
    
  bc->SetCashDistribution(pInterests);
  
  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(false);
  
  conv->AddConversionPeriod(
    shared_ptr<ConversionPeriod>(new ConversionPeriod(issueDate, maturityDate, 1)));

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
 
  Date startCallDate(issueDate);
  Date endCallDate(maturityDate);  
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike
                             (startCallDate, endCallDate, 1) );
  pCallSchedule->AddCallPeriod( pCallPeriod );    
    
  pConvertibleBond->SetCallSchedule(pCallSchedule);

  pConvertibleBond->SetSessionData(pSessionData);

  shared_ptr<Numeraire> pCurrency(new Numeraire("USD"));
  pConvertibleBond->SetNumeraire(pCurrency);

  return pConvertibleBond;
}

//*************************************************************************

int main()
{
  try {
  
  bool
    bComputeRho = false,
    bComputeVega = false,
    bComputeSurface = false;
  
  shared_ptr<SessionData> pSessionData = MakeSessionData(100.0);

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  // pModel->SetDebugOutputFile("ihgcrosscurrency.xml");

  double dVol = 0.5;
  pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

  /*
  size_t nNbHRTimes = 2;
  static const Date pTimes[2] = { 38000, 40000};
  static const double pdValues[2] = { 0.02, .02 };

  pModel->SetHazardRate
          ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );
  */

  double dLambda = .2;
  pModel->SetHazardRate( 
    shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda)) );

  shared_ptr<finance::ConvertibleBond> pCB = InitCB(pSessionData);
  
  shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
  flags->SetComputeRho(bComputeRho);
  flags->SetComputeVega(bComputeVega);
  flags->SetComputeSurface(bComputeSurface);
  flags->SetAnalysisDate( Date("2003/02/01"));

  pCB->SetComputationalFlags(flags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);

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

